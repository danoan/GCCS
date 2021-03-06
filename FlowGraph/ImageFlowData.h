#ifndef SEGBYCUT_IMAGEFLOWDATA_H
#define SEGBYCUT_IMAGEFLOWDATA_H

#include <DGtal/images/ImageContainerBySTLVector.h>
#include "DGtal/helpers/StdDefs.h"

#include "ConnectorSeedRange.h"
#include "ConnectorFunctor.h"

#include "../ImageProc/imageProc.h"

class ImageFlowData{

public:
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

    typedef Curve::ConstIterator SCellIteratorType;
    typedef DGtal::Circulator<SCellIteratorType> SCellCirculator;
    typedef ConnectorSeedRange<KSpace,SCellCirculator> ConnectorSeedRangeType;
    typedef ConnectorSeedToGluedCurveRange< ConnectorSeedRangeType::ConnectorSeedType > SeedToGluedCurveRangeFunctor;
    typedef SeedToGluedCurveRangeFunctor::GCIterator SCellGluedCurveIterator;

    typedef ConstRangeAdapter<  ConnectorSeedRangeType::ConnectorSeedIteratorType,
                                SeedToGluedCurveRangeFunctor,
                                SeedToGluedCurveRangeFunctor::Output > GluedCurveSetRange;

    typedef GluedCurveSetRange::ConstIterator GluedCurveIteratorPair;


    enum CurveType{
        OriginalCurve,DilatedCurve,ErodedCurve
    };

    enum FlowMode{
        DilationOnly,
        ErosionOnly,
        DilationErosion
    };

    typedef struct _CurveData{
        Curve curve;
        CurveType curveType;

        SCellCirculator curveCirculator;
    } CurveData;


    class CurvePair{
    public:
        CurvePair(CurveData& intCurveData,
                  CurveData& extCurveData,
                  KSpace& KImage,
                  SeedToGluedCurveRangeFunctor& stgcF):intCurveData(intCurveData),
                                                       extCurveData(extCurveData),
                                                       KImage(KImage),
                                                       stgcF(stgcF),
                                                       seedRange(KImage,
                                                                 intCurveData.curveCirculator,
                                                                 extCurveData.curveCirculator),
                                                       gcsRange(seedRange.begin(),
                                                                seedRange.end(),
                                                                stgcF){};

        CurvePair(const CurvePair& other):intCurveData(other.intCurveData),
                                          extCurveData(other.extCurveData),
                                          KImage(other.KImage),
                                          stgcF(other.stgcF),
                                          seedRange(KImage,
                                                    intCurveData.curveCirculator,
                                                    extCurveData.curveCirculator),
                                          gcsRange(seedRange.begin(),
                                                   seedRange.end(),
                                                   stgcF){};


        /*The right thing to do is to implement an iterator at the level of
         * ImageFlowData. Tha iterator would traverse all the gcsRanges.
         * */
        CurvePair& operator=(const CurvePair& other)
        {
            throw "CurvePair assignment is forbidden";
        }


        GluedCurveSetRange::ConstIterator gcsRangeBegin() const{ return gcsRange.begin(); }
        GluedCurveSetRange::ConstIterator gcsRangeEnd() const{ return gcsRange.end(); }

        GluedCurveSetRange getGCSRange() const{ return gcsRange; }


        CurveData& intCurveData;
        CurveData& extCurveData;

        KSpace& KImage;
        SeedToGluedCurveRangeFunctor& stgcF;

    private:
        ConnectorSeedRangeType seedRange;
        GluedCurveSetRange gcsRange;
    };

    typedef std::vector<CurvePair>::const_iterator CurvePairIterator;


public:

    ImageFlowData(Image2D image,Image2D refImage);
    ImageFlowData(const ImageFlowData& other);
    ImageFlowData& operator=(const ImageFlowData& other);

    void setInitParameters(FlowMode fm, int gcLength);
    void setCurveData(FlowMode fm);

    void init(FlowMode fm,int gluedCurveLength);
    void init(FlowMode fm,int gluedCurveLength,Curve outerCurve);

    std::vector<CurveData>::const_iterator curveDataBegin(){ return curvesVector.begin(); }
    std::vector<CurveData>::const_iterator curveDataEnd(){ return curvesVector.end(); }

    CurvePairIterator curvePairBegin(){ return curvesPairVector.begin(); }
    CurvePairIterator curvePairEnd(){ return curvesPairVector.end(); }

    KSpace& getKSpace(){return KImage;}

    FlowMode getFlowMode(){return flowMode;}

    Curve& getMostInnerCurve();
    Curve& getMostOuterCurve();
    Curve& getBaseCurve();

    int getGluedCurveLength(){ return gluedCurveLength; }
    int getConsecutiveGluedPairDistance(){ return consecutiveGluedPairDistance; }
    int getDiffDistance(){ return diffDistance; }

    Image2D& getRefImage(){ return refImage; };
    Image2D& getOriginalImage(){ return originalImage; }
    Image2D& getDilatedImage(){ return dilatedImage; }
    Image2D& getErodedImage(){ return erodedImage; }

    Image2D& getMostInnerImage(){
        switch(flowMode) {
            case DilationOnly:
                return originalImage;
            default:
                return erodedImage;
        }
    }

private:
    friend class ImageFlowDataDebug;

    bool itIsInitialized;
    FlowMode flowMode;
    int gluedCurveLength;
    int consecutiveGluedPairDistance;
    int diffDistance;

    KSpace KImage;
    Image2D refImage,originalImage,dilatedImage,erodedImage;

    std::vector<CurveData> curvesVector;
    std::vector<CurvePair> curvesPairVector;

    SeedToGluedCurveRangeFunctor stgcF; //Must exist during the lifetime of GluedCurveSetRange


    CurveData& addNewCurve(CurveType ct);
    CurveData& getCurveData(CurveType ct);

    void initRange(CurveData& intCurveData, CurveData& extCurveData);
    void registerCirculator(CurveData& curveData);
};


#endif //SEGBYCUT_IMAGEFLOWDATA_H
