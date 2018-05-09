#ifndef SEGBYCUT_ADAPTATIVESCALEFLOW_H
#define SEGBYCUT_ADAPTATIVESCALEFLOW_H

#include "../../FlowGraph/ImageFlowData.h"
#include "../../FlowGraph/FlowGraphBuilder.h"
#include "../../FlowGraph/weightSettings.h"
#include "../../FlowGraph/FlowGraphQuery.h"
#include "../../FlowGraph/ImageUpdater.h"

class AdaptativeScaleFlow
{
public:
    enum IterationType{EvenIteration=0,OddIteration=1};

    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
    typedef FlowGraphBuilder::LinelWeightMap LinelWeightMap;

public:
    AdaptativeScaleFlow(int gluedCurveLength):gluedCurveLength(gluedCurveLength),
                                              iteration(EvenIteration),
                                              _cutValue(0){}


    void operator()(Image2D& imageIn, Image2D& imageOut)
    {

        switch(iteration)
        {
            case EvenIteration:
                evenIteration(imageIn,imageOut);
                iteration = OddIteration;
                break;
            case OddIteration:
                oddIteration(imageIn,imageOut);
                iteration = EvenIteration;
                break;
        }

    }

    inline double cutValue(){return _cutValue;}

private:
    void evenIteration(Image2D& imageIn, Image2D& imageOut)
    {
        //DilationFlow -> Cut
        ImageFlowData imageFlowData(imageIn);
        imageFlowData.init(ImageFlowData::DilationOnly,gluedCurveLength);

        LinelWeightMap weightMap;
        setArcsWeight(imageFlowData,weightMap);

        FlowGraph fg;
        FlowGraphBuilder fgb(fg,imageFlowData,weightMap);

        _cutValue = FlowGraphQuery::cutValue(fg);

        ImageUpdater imageUpdater(imageFlowData.getKSpace());
        imageUpdater(imageOut,fg,imageIn);

    }

    void oddIteration(Image2D& imageIn, Image2D& imageOut)
    {
        //2x Scaling -> ErosionFlow -> Cut -> 0.5 Scaling
        Domain domainIn = imageIn.domain();
        Domain domainResized( Point(domainIn.lowerBound()),Point(domainIn.upperBound()*2));

        Image2D resizedImageIn(domainResized);

        ImageProc::resize(imageIn,resizedImageIn,2);

        ImageFlowData imageFlowData(resizedImageIn);
        imageFlowData.init(ImageFlowData::ErosionOnly,gluedCurveLength*2);

        LinelWeightMap weightMap;
        setArcsWeight(imageFlowData,weightMap);

        FlowGraph fg;
        FlowGraphBuilder fgb(fg, imageFlowData,weightMap);

        _cutValue = FlowGraphQuery::cutValue(fg);

        Image2D resizedImageOut(domainResized);
        ImageUpdater imageUpdater(imageFlowData.getKSpace());
        imageUpdater(resizedImageOut,fg,resizedImageIn);

        ImageProc::resize(resizedImageOut,imageOut,0.5);

    };

private:
    int iteration;
    int gluedCurveLength;

    double _cutValue;
};

#endif //SEGBYCUT_ADAPTATIVESCALEFLOW_H
