#ifndef SEGBYCUT_BOUNDEDFLOW_H
#define SEGBYCUT_BOUNDEDFLOW_H

#include "../../FlowGraph/ImageFlowData.h"
#include "../../FlowGraph/FlowGraphBuilder.h"
#include "../../FlowGraph/FlowGraphDebug.h"
#include "BoundingBox.h"
#include "InnerProjection.h"
#include "../RefundFlow/RefundFlow.h"

class BoundedFlow
{
public:
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
    typedef FlowGraphBuilder::LinelWeightMap LinelWeightMap;

public:
    BoundedFlow(int gluedCurveLength):gluedCurveLength(gluedCurveLength)
    {}


    void operator()(Image2D imageIn, Image2D& imageOut)
    {
        Curve optimizationBoundCurve;
//        boundingBox(boundingBoxCurve,imageIn);
        enveloppe(optimizationBoundCurve,imageIn);

        imageOut=imageIn;
        while( dilationStepRefundFlow(imageOut,optimizationBoundCurve) )
        {
            ImageProc::erode(imageOut,imageOut,1,cv::MORPH_CROSS);
        }
    }

    double cutValue(){ return _cutValue; }

private:

    bool dilationStepRefundFlow(Image2D& image, Curve& optimizationBound)
    {
        LinelWeightMap exampleWeight;

        Curve projectedCurve;
        innerProjectedCurve(projectedCurve, image,optimizationBound);

        {
            ImageFlowData temp(image);
            temp.init(ImageFlowData::DilationOnly,5);

            Board2D board;
            board << DGtal::CustomStyle(projectedCurve.begin()->className(),
                                        new DGtal::CustomColors(DGtal::Color::Red, DGtal::Color::Green));
            board << temp.getMostInnerCurve();
            board.saveEPS("innerCurve.eps");

            board << DGtal::CustomStyle(projectedCurve.begin()->className(),
                                        new DGtal::CustomColors(DGtal::Color::Green, DGtal::Color::Green));
            board << projectedCurve;
            board.saveEPS("projectedCurvePlusInner.eps");
        }

        ImageFlowData imf(image);
        imf.init(ImageFlowData::DilationOnly, gluedCurveLength,projectedCurve);

        {
            //Graph Draw

            FlowGraph fg;
            FlowGraphBuilder fgb(fg, imf, exampleWeight);

            FlowGraphDebug fgd(fg);
            fgd.drawFlowGraph("../cmake-build-debug", "test.eps");
        }

        RefundFlow refundFlow(imf,1);

        _cutValue = refundFlow.run(0);
        image = refundFlow.outputImage();

        std::cout << _cutValue << std::endl;

        return false;
    }

    bool dilationStep(Image2D& image, Curve& optimizationBound)
    {
        double previousCutValue=-10;
        _cutValue=-12;

        Image2D previousImage = image;
        int k = 0;
        LinelWeightMap exampleWeight;
        while( fabs(previousCutValue-_cutValue)>0.0001 )
        {
            previousCutValue = _cutValue;

            Curve projectedCurve;
            innerProjectedCurve(projectedCurve, image,optimizationBound);

            ImageFlowData temp(image);
            temp.init(ImageFlowData::DilationOnly,5);

            Board2D board;
            board << DGtal::CustomStyle(projectedCurve.begin()->className(),
                                        new DGtal::CustomColors( DGtal::Color::Red, DGtal::Color::Green ) );
            board << temp.getMostInnerCurve();
            board.saveEPS("innerCurve.eps");

            board << DGtal::CustomStyle(projectedCurve.begin()->className(),
                                        new DGtal::CustomColors( DGtal::Color::Green, DGtal::Color::Green ) );
            board << projectedCurve;
            board.saveEPS("projectedCurvePlusInner.eps");

            ImageFlowData imf(image);
            imf.init(ImageFlowData::DilationOnly, gluedCurveLength,projectedCurve);

            {
                //Graph Draw

                FlowGraph fg;
                FlowGraphBuilder fgb(fg, imf, exampleWeight);

                FlowGraphDebug fgd(fg);
                fgd.drawFlowGraph("../cmake-build-debug", "test.eps");
            }

            LinelWeightMap weightMap;
            setArcsWeight(imf, weightMap);
            exampleWeight = weightMap;

            FlowGraph fg;
            FlowGraphBuilder fgb(fg, imf, weightMap);

            _cutValue = FlowGraphQuery::cutValue(fg);

            ImageUpdater imageUpdater(imf.getKSpace());
            imageUpdater(image,fg,previousImage);
            previousImage=image;

            ++k;

            std::cout << _cutValue << std::endl;
        }

        return k>=2;
    }


    void boundingBox(Curve& boundingBoxCurve, Image2D& imageIn)
    {
        ImageFlowData imageFlowData(imageIn);
        imageFlowData.init(ImageFlowData::DilationOnly,gluedCurveLength);

        Point lowerLeft,upperRight;
        KSpace& KImage = imageFlowData.getKSpace();
        Curve& outerCurve = imageFlowData.getMostOuterCurve();

        BoundingBox BB(KImage);
        BB(lowerLeft,upperRight,boundingBoxCurve,outerCurve);
    }

    void enveloppe(Curve& boundingBoxCurve, Image2D& imageIn)
    {
        Image2D imageOut = imageIn;
        ImageProc::dilate(imageOut,imageIn,2);
        ImageProc::computeBoundaryCurve(imageOut,boundingBoxCurve,100);
    }

    void innerProjectedCurve(Curve& projectedCurve,Image2D& imageIn, Curve& boundingBox)
    {
        ImageFlowData imageFlowData(imageIn);
        imageFlowData.init(ImageFlowData::DilationOnly,gluedCurveLength);

        InnerProjection IP(imageFlowData.getKSpace());
        IP(projectedCurve,imageFlowData.getMostOuterCurve(),boundingBox);
    }

private:
    int gluedCurveLength;
    double _cutValue;
};

#endif //SEGBYCUT_BOUNDEDFLOW_H
