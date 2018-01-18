#include "utils.h"

void drawCurvatureMap(std::vector<Z2i::SCell>::const_iterator begin,
                      std::vector<Z2i::SCell>::const_iterator end,
                      double& cmin,
                      double& cmax,
                      Board2D& board,
                      std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> values;

    auto it=begin;
    do{
        values.push_back( weightMap[*it] );
        cmin = weightMap[*it]<cmin?weightMap[*it]:cmin;
        cmax = weightMap[*it]>cmax?weightMap[*it]:cmax;
        ++it;
    }while(it!=end);

    draw(values,begin,end,board,cmin,cmax);
}

void curvatureEstimatorsGridCurve(UtilsTypes::Curve::ConstIterator begin,
                                  UtilsTypes::Curve::ConstIterator end,
                                  UtilsTypes::KSpace& KImage,
                                  std::vector<double>& estimations)
{
    typedef UtilsTypes::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

    typedef ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
                               AdapterFunctor,
                               AdapterFunctor::Output > RangeAdapter;

    AdapterFunctor myFunctor(KImage);
    RangeAdapter range(begin,end,myFunctor);


    if(Patch::useDGtal){
        Patch::estimationsDGtalMCMSECurvature(range.c(),range.c(),estimations);
    }else{
        Patch::estimationsPatchMCMSECurvature(range.c(),range.c(),estimations);
    }

}

void curvatureEstimatorsGluedCurve(UtilsTypes::SCellGluedCurveIterator begin,
                                   UtilsTypes::SCellGluedCurveIterator end,
                                   UtilsTypes::KSpace& KImage,
                                   std::vector<double>& estimations)
{
    typedef UtilsTypes::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

    typedef UtilsTypes::ConstRangeAdapter< UtilsTypes::SCellGluedCurveIterator,
                                           AdapterFunctor,
                                           AdapterFunctor::Output > RangeAdapter;

    AdapterFunctor myFunctor(KImage);
    RangeAdapter range(begin,end,myFunctor);

    if(Patch::useDGtal){
        Patch::estimationsDGtalMCMSECurvature(range.begin(),range.end(),estimations);
    }else{
        Patch::estimationsPatchMCMSECurvature(range.begin(),range.end(),estimations);
    }

}


void curvatureEstimatorsConnections(UtilsTypes::GluedCurveIteratorPair begin,
                                    UtilsTypes::GluedCurveIteratorPair end,
                                    UtilsTypes::KSpace& KImage,
                                    unsigned int gluedCurveLength,
                                    std::vector<double>& estimations)
{

    typedef UtilsTypes::functors::SCellToIncidentPoints<UtilsTypes::KSpace> SCellToIncPointsFunctor;
    SCellToIncPointsFunctor myFunc(KImage);

    std::vector<double> partialEstimations;
    UtilsTypes::GluedCurveIteratorPair it= begin;
    do{
        UtilsTypes::SCellGluedCurveIterator gcBegin = it->first;
        UtilsTypes::SCellGluedCurveIterator gcEnd = it->second;

        UtilsTypes::GluedCurveIncidentPointsRange gcipRange(gcBegin,
                                                                   gcEnd,
                                                                   myFunc );

        if(Patch::useDGtal){
            Patch::estimationsDGtalMCMSECurvature(gcipRange.begin(),
                                                  gcipRange.end(),
                                                  partialEstimations);
        }else{
            Patch::estimationsPatchMCMSECurvature(gcipRange.begin(),
                                                  gcipRange.end(),
                                                  partialEstimations);
        }

        estimations.push_back( partialEstimations[gluedCurveLength] );
        partialEstimations.clear();

        ++it;
    }while(it!=end);

}



void tangentEstimatorsGridCurve(UtilsTypes::Curve::ConstIterator begin,
                                UtilsTypes::Curve::ConstIterator end,
                                UtilsTypes::KSpace& KImage,
                                std::vector< UtilsTypes::TangentVector >& estimationsTangent )
{
    typedef ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
                               functors::SCellToPoint<UtilsTypes::KSpace>,
                               UtilsTypes::KSpace::Point > RangeAdapter;

    functors::SCellToPoint<UtilsTypes::KSpace> myFunctor = functors::SCellToPoint<UtilsTypes::KSpace>(KImage);
    RangeAdapter range(begin,end,myFunctor);


    if(Patch::useDGtal){
        Patch::estimationsDGtalMCMSETangent(range.c(),range.c(),estimationsTangent);
    }else{
        Patch::estimationsPatchMCMSETangent(range.c(),range.c(),estimationsTangent);
    }


}

void tangentEstimatorsGluedCurve(UtilsTypes::Curve::ConstIterator begin,
                                 UtilsTypes::Curve::ConstIterator end,
                                 UtilsTypes::KSpace& KImage,
                                std::vector< UtilsTypes::TangentVector >& estimationsTangent )
{
    typedef UtilsTypes::ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
                                           UtilsTypes::functors::SCellToPoint<KSpace>,
                                           UtilsTypes::KSpace::Point> RangeAdapter;

    UtilsTypes::functors::SCellToPoint<KSpace> myFunctor = UtilsTypes::functors::SCellToPoint<KSpace>(KImage);
    RangeAdapter range(begin,end,myFunctor);

    if(Patch::useDGtal){
        Patch::estimationsDGtalMCMSETangent(range.begin(),range.end(),estimationsTangent);
    }else{
        Patch::estimationsPatchMCMSETangent(range.begin(),range.end(),estimationsTangent);
    }

}


void tangentEstimatorsConnections(UtilsTypes::GluedCurveIteratorPair begin,
                                  UtilsTypes::GluedCurveIteratorPair end,
                                  KSpace& KImage,
                                  unsigned int gluedCurveLength,
                                  std::vector< UtilsTypes::TangentVector >& estimations)
{
    typedef UtilsTypes::functors::SCellToPoint<UtilsTypes::KSpace> SCellToPointFunctor;

    typedef UtilsTypes::ConstRangeAdapter<  UtilsTypes::SCellGluedCurveIterator,
                                            SCellToPointFunctor,
                                            UtilsTypes::KSpace::Point > GluedCurvePointsRange;


    SCellToPointFunctor myFunc(KImage);
    std::vector< UtilsTypes::TangentVector > partialEstimations;

    UtilsTypes::GluedCurveIteratorPair it= begin;
    do{
        UtilsTypes::SCellGluedCurveIterator begin = it->first;
        UtilsTypes::SCellGluedCurveIterator end = it->second;

        GluedCurvePointsRange gcipRange(begin,
                                        end,
                                        myFunc );

        if(Patch::useDGtal){
            Patch::estimationsDGtalMCMSETangent(gcipRange.begin(),
                                                gcipRange.end(),
                                                partialEstimations);
        }else{
            Patch::estimationsPatchMCMSETangent(gcipRange.begin(),
                                                gcipRange.end(),
                                                partialEstimations);
        }

        estimations.push_back( partialEstimations[gluedCurveLength] );
        partialEstimations.clear();

        ++it;
    }while(it!=end);

}


void setKImage(std::string imgFilePath,KSpace& KImage)
{
    UtilsTypes::Image2D image = UtilsTypes::GenericReader<UtilsTypes::Image2D>::import(imgFilePath);
    KImage.init(image.domain().lowerBound(),image.domain().upperBound(),true);
}

void setCurves(std::string imgFilePath,
               Curve& intCurve,
               Curve& extCurve)
{
    UtilsTypes::Image2D image = GenericReader<UtilsTypes::Image2D>::import(imgFilePath);
    KSpace KImage;

    computeBoundaryCurve(intCurve,KImage,image);

    UtilsTypes::Image2D dilatedImage(image.domain());
    dilate(dilatedImage,imgFilePath,1);
    computeBoundaryCurve(extCurve,KImage,dilatedImage);
}

UtilsTypes::ConnectorSeedRangeType getSeedRange(KSpace& KImage,
                                                         Curve& intCurve,
                                                         Curve& extCurve)
{
    UtilsTypes::SCellCirculator intCirculator( intCurve.begin(),
                                                        intCurve.begin(),
                                                        intCurve.end() );

    UtilsTypes::SCellCirculator extCirculator( extCurve.begin(),
                                                        extCurve.begin(),
                                                        extCurve.end() );


    UtilsTypes::ConnectorSeedRangeType seedRange(  KImage,
                                                            intCirculator,
                                                            extCirculator );

    return seedRange;
}


