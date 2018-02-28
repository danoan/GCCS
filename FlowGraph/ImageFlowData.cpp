#include "ImageFlowData.h"


ImageFlowData::ImageFlowData(Image2D image):originalImage(image),
                                            dilatedImage(image.domain()),
                                            erodedImage(image.domain()),
                                            itIsInitialized(false),
                                            stgcF()
{
    KImage.init(image.domain().lowerBound(),image.domain().upperBound(),true);
}

ImageFlowData::ImageFlowData(const ImageFlowData& other):originalImage(other.originalImage),
                                                         dilatedImage(other.dilatedImage),
                                                         erodedImage(other.erodedImage),
                                                         KImage(other.KImage),
                                                         itIsInitialized(other.itIsInitialized),
                                                         stgcF()
{
    if(other.itIsInitialized){
        init(other.flowMode,other.gluedCurveLength);
    }

}


ImageFlowData& ImageFlowData::operator=(const ImageFlowData& other)
{
    originalImage = other.originalImage;
    erodedImage = other.erodedImage;
    dilatedImage = other.dilatedImage;
    KImage = other.KImage;

    init(other.flowMode,other.gluedCurveLength);

    return *this;
}

ImageFlowData::CurveData& ImageFlowData::addNewCurve(CurveType ct)
{
    CurveData cd;
    cd.curve = Curve();
    cd.curveType = ct;

    curvesVector.push_back( cd );
    return curvesVector.back();
}

ImageFlowData::CurveData& ImageFlowData::getCurveData(CurveType ct)
{
    int i=0;
    for(;i<curvesVector.size() && curvesVector[i].curveType!=ct;++i);

    if(i==curvesVector.size()) throw "No curve of requested type was found";

    return curvesVector.at(i);
}

void ImageFlowData::init(FlowMode fm, int gcLength)
{
    itIsInitialized = true;
    flowMode = fm;
    gluedCurveLength = gcLength;
    consecutiveGluedPairDistance = 4;
    diffDistance = 8;

    stgcF.init(gcLength);


    /* Since I am using references to vector elements, I must make sure
     * the vector holds enough places for all the elements I am going to
     * need, otherwise, the vector class will reallocate the elements
     * during a push_back, and the reference I saved will be lost.
     * */
    curvesVector.reserve(4);
    curvesPairVector.reserve(4);
    {
        CurveData &cd = addNewCurve(CurveType::OriginalCurve);
        ImageProc::computeBoundaryCurve(originalImage, cd.curve, 100);
        registerCirculator(cd);
    }

    if(fm==FlowMode::DilationOnly || fm==FlowMode::DilationErosion)
    {
        CurveData &cd = addNewCurve(CurveType::DilatedCurve);
        ImageProc::dilate(dilatedImage, originalImage, 1);

        ImageProc::computeBoundaryCurve(dilatedImage, cd.curve, 100);
        registerCirculator(cd);
    }

    if(fm==FlowMode::ErosionOnly || fm==FlowMode::DilationErosion)
    {
        CurveData &cd = addNewCurve(CurveType::ErodedCurve);
        ImageProc::erode(erodedImage, originalImage, 1);

        ImageProc::computeBoundaryCurve(erodedImage, cd.curve, 100);
        registerCirculator(cd);
    }

    if(fm==DilationOnly){
        CurveData& a = getCurveData(CurveType::OriginalCurve);
        CurveData& b = getCurveData(CurveType::DilatedCurve);

        initRange(getCurveData(CurveType::OriginalCurve),
                  getCurveData(CurveType::DilatedCurve));

    }else if(fm==ErosionOnly){
        initRange(getCurveData(CurveType::ErodedCurve),
                  getCurveData(CurveType::OriginalCurve));

    }else{
        initRange(getCurveData(CurveType::ErodedCurve),
                  getCurveData(CurveType::OriginalCurve));

        initRange(getCurveData(CurveType::OriginalCurve),
                  getCurveData(CurveType::DilatedCurve));
    }

}

void ImageFlowData::registerCirculator(CurveData& curveData)
{
    SCellCirculator circulator( curveData.curve.begin(),
                                curveData.curve.begin(),
                                curveData.curve.end() );

    curveData.curveCirculator = circulator;

}

void ImageFlowData::initRange(CurveData& intCurveData, CurveData& extCurveData)
{
    CurvePair cp(intCurveData,
                 extCurveData,
                 KImage,
                 stgcF);

    curvesPairVector.push_back(cp);
}

Curve& ImageFlowData::getMostInnerCurve()
{
    switch(flowMode){
        case FlowMode::DilationOnly:
            return getCurveData(CurveType::OriginalCurve).curve;
        case FlowMode::ErosionOnly:
            return getCurveData(CurveType::ErodedCurve).curve;
        case FlowMode::DilationErosion:
            return getCurveData(CurveType::ErodedCurve).curve;
        default:
            throw "No flow mode defined. Impossible to define the most inner curve";
    }
}

Curve& ImageFlowData::getMostOuterCurve()
{
    switch(flowMode){
        case FlowMode::DilationOnly:
            return getCurveData(CurveType::DilatedCurve).curve;
        case FlowMode::ErosionOnly:
            return getCurveData(CurveType::OriginalCurve).curve;
        case FlowMode::DilationErosion:
            return getCurveData(CurveType::DilatedCurve).curve;
        default:
            throw "No flow mode defined. Impossible to define the most inner curve";
    }
}