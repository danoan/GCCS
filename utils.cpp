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
                                  std::vector<double>& estimations,
                                  bool closedCurve)
{
    typedef UtilsTypes::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

    typedef ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
                               AdapterFunctor,
                               AdapterFunctor::Output > RangeAdapter;

    Curve negativeCurve;
    invertCurve(KImage,
                begin,
                end,
                negativeCurve);

    AdapterFunctor myFunctor(KImage);
    RangeAdapter rangePositiveCurve(begin,
                                    end,
                                    myFunctor);

    RangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                    negativeCurve.end(),
                                    myFunctor);


    std::vector<double> positiveEstimations;
    std::vector<double> negativeEstimations;

    if(closedCurve) {
        Patch::estimationsDGtalMCMSECurvature(rangePositiveCurve.c(),
                                              rangePositiveCurve.c(),
                                              positiveEstimations);

        Patch::estimationsDGtalMCMSECurvature(rangeNegativeCurve.c(),
                                              rangeNegativeCurve.c(),
                                              negativeEstimations);
    }else{
        Patch::estimationsDGtalMCMSECurvature(rangePositiveCurve.begin(),
                                              rangePositiveCurve.end(),
                                              positiveEstimations);

        Patch::estimationsDGtalMCMSECurvature(rangeNegativeCurve.begin(),
                                              rangeNegativeCurve.end(),
                                              negativeEstimations);
    }

    if(Patch::solveShift){
        //Solve Shift
        positiveEstimations.push_back(positiveEstimations[0]);
        positiveEstimations.erase(positiveEstimations.begin());

        negativeEstimations.push_back(negativeEstimations[0]);
        negativeEstimations.erase(negativeEstimations.begin());
    }


    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( ( fabs(positiveEstimations[ip]) + fabs(negativeEstimations[nL-ip]) )/2.0 );
        ++ip;
    }while(ip<=nL);

}

void curvatureEstimatorsGluedCurve(UtilsTypes::SCellGluedCurveIterator begin,
                                   UtilsTypes::SCellGluedCurveIterator end,
                                   UtilsTypes::KSpace& KImage,
                                   std::vector<double>& estimations)
{
    typedef UtilsTypes::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

    typedef UtilsTypes::ConstRangeAdapter< UtilsTypes::SCellGluedCurveIterator,
                                           AdapterFunctor,
                                           AdapterFunctor::Output > GluedRangeAdapter;


    typedef ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
            AdapterFunctor,
            AdapterFunctor::Output > RangeAdapter;


    AdapterFunctor myFunctor(KImage);
    GluedRangeAdapter rangePositiveCurve(begin,
                                    end,
                                    myFunctor);

    Curve negativeCurve;
    invertCurve(KImage,
                begin,
                end,
                negativeCurve);

    RangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                    negativeCurve.end(),
                                    myFunctor);

    std::vector<double> positiveEstimations;
    std::vector<double> negativeEstimations;

    Patch::estimationsDGtalMCMSECurvature(rangePositiveCurve.begin(),
                                          rangePositiveCurve.end(),
                                          positiveEstimations);

    Patch::estimationsDGtalMCMSECurvature(rangeNegativeCurve.begin(),
                                          rangeNegativeCurve.end(),
                                          negativeEstimations);

    if(Patch::solveShift){
        //Solve Shift
        positiveEstimations.push_back(positiveEstimations[0]);
        positiveEstimations.erase(positiveEstimations.begin());

        negativeEstimations.push_back(negativeEstimations[0]);
        negativeEstimations.erase(negativeEstimations.begin());
    }


    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( ( fabs(positiveEstimations[ip]) + fabs(negativeEstimations[nL-ip]) )/2.0 );
        ++ip;
    }while(ip<=nL);
}

void normalizeAroundNeighbors(std::vector<double>& v, int radius)
{
    std::vector<double> temp = v;
    double s;
    for(int i=radius;i<v.size()-radius;++i){
        s=0;
        for(int j=-radius;j<=radius;++j){
            s+=temp[i+j];
        }
        v[i] = s/(2*radius+1);
    }
}

void normalizeAroundNeighbors(std::vector<UtilsTypes::TangentVector>& v, int radius)
{
    std::vector<UtilsTypes::TangentVector> temp = v;
    UtilsTypes::TangentVector s;
    for(int i=radius;i<v.size()-radius;++i){
        s = UtilsTypes::TangentVector(0,0);
        for(int j=-radius;j<=radius;++j){
            s+=temp[i+j];
            s = s.getNormalized();
        }
        v[i] = s;
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

        curvatureEstimatorsGluedCurve(gcBegin,
                                      gcEnd,
                                      KImage,
                                      partialEstimations);

        for(int i=0;i<gcBegin.numberOfConnectors();i++){
            estimations.push_back( partialEstimations[gluedCurveLength+i] );
        }

        partialEstimations.clear();

        ++it;
    }while(it!=end);

}



void tangentEstimatorsGridCurve(UtilsTypes::Curve::ConstIterator begin,
                                UtilsTypes::Curve::ConstIterator end,
                                UtilsTypes::KSpace& KImage,
                                std::vector< UtilsTypes::TangentVector >& estimationsTangent,
                                bool closedCurve)
{
    typedef ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
                               functors::SCellToPoint<UtilsTypes::KSpace>,
                               UtilsTypes::KSpace::Point > RangeAdapter;

    functors::SCellToPoint<UtilsTypes::KSpace> myFunctor = functors::SCellToPoint<UtilsTypes::KSpace>(KImage);
    RangeAdapter rangePositiveCurve(begin,
                                    end,
                                    myFunctor);

    Curve negativeCurve;
    invertCurve(KImage,begin,end,negativeCurve);

    RangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                    negativeCurve.end(),
                                    myFunctor);

    std::vector<UtilsTypes::TangentVector> positiveEstimations;
    std::vector<UtilsTypes::TangentVector> negativeEstimations;
    if(closedCurve) {

        Patch::estimationsDGtalMCMSETangent(rangePositiveCurve.c(),
                                            rangePositiveCurve.c(),
                                            positiveEstimations);

        Patch::estimationsDGtalMCMSETangent(rangeNegativeCurve.c(),
                                            rangeNegativeCurve.c(),
                                            negativeEstimations);
    }else{
        Patch::estimationsDGtalMCMSETangent(rangePositiveCurve.begin(),
                                            rangePositiveCurve.end(),
                                            positiveEstimations);

        Patch::estimationsDGtalMCMSETangent(rangeNegativeCurve.begin(),
                                            rangeNegativeCurve.end(),
                                            negativeEstimations);
    }

    if(Patch::solveShift){
        //Solve Shift
        positiveEstimations.push_back(positiveEstimations[0]);
        positiveEstimations.erase(positiveEstimations.begin());

        negativeEstimations.push_back(negativeEstimations[0]);
        negativeEstimations.erase(negativeEstimations.begin());
    }

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimationsTangent.push_back( (positiveEstimations[ip]-negativeEstimations[nL-ip]).getNormalized() );
        ++ip;
    }while(ip<=nL);


}

void tangentEstimatorsGluedCurve(UtilsTypes::SCellGluedCurveIterator begin,
                                 UtilsTypes::SCellGluedCurveIterator end,
                                 UtilsTypes::KSpace& KImage,
                                std::vector< UtilsTypes::TangentVector >& estimationsTangent )
{
    typedef UtilsTypes::ConstRangeAdapter< UtilsTypes::SCellGluedCurveIterator,
            UtilsTypes::functors::SCellToPoint<KSpace>,
            UtilsTypes::KSpace::Point> GluedRangeAdapter;

    UtilsTypes::functors::SCellToPoint<KSpace> myFunctor = UtilsTypes::functors::SCellToPoint<KSpace>(KImage);
    GluedRangeAdapter rangePositiveCurve(begin,
                                         end,
                                         myFunctor);

    typedef UtilsTypes::ConstRangeAdapter< UtilsTypes::Curve::ConstIterator,
            UtilsTypes::functors::SCellToPoint<KSpace>,
            UtilsTypes::KSpace::Point> RangeAdapter;

    Curve negativeCurve;
    invertCurve(KImage,
                begin,
                end,
                negativeCurve);


    RangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                    negativeCurve.end(),
                                    myFunctor);

    std::vector<UtilsTypes::TangentVector> positiveEstimations;
    Patch::estimationsDGtalMCMSETangent(rangePositiveCurve.begin(),
                                        rangePositiveCurve.end(),
                                        positiveEstimations);

    std::vector<UtilsTypes::TangentVector> negativeEstimations;
    Patch::estimationsDGtalMCMSETangent(rangeNegativeCurve.begin(),
                                        rangeNegativeCurve.end(),
                                        negativeEstimations);

    if(Patch::solveShift){
        //Solve Shift
        positiveEstimations.push_back(positiveEstimations[0]);
        positiveEstimations.erase(positiveEstimations.begin());

        negativeEstimations.push_back(negativeEstimations[0]);
        negativeEstimations.erase(negativeEstimations.begin());
    }

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimationsTangent.push_back( (positiveEstimations[ip]-negativeEstimations[nL-ip]).getNormalized() );
        ++ip;
    }while(ip<=nL);

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

        tangentEstimatorsGluedCurve(begin,
                                    end,
                                    KImage,
                                    partialEstimations);

        for(int i=0;i<begin.numberOfConnectors();i++){
            estimations.push_back( partialEstimations[gluedCurveLength+i] );
        }

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




