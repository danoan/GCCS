#include "utils.h"

namespace UtilsTypes
{
    std::function< double(double) > toDouble = [](double x){return x;};
};


void drawCurvatureMap(std::vector<Z2i::SCell>::const_iterator begin,
                      std::vector<Z2i::SCell>::const_iterator end,
                      double& cmin,
                      double& cmax,
                      Board2D& board,
                      std::map<Z2i::SCell,double>& weightMap)
{
    if(begin==end) return;
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
        Patch::curvatureEstimation(rangePositiveCurve.c(),
                                   rangePositiveCurve.c(),
                                   positiveEstimations);

        Patch::curvatureEstimation(rangeNegativeCurve.c(),
                                   rangeNegativeCurve.c(),
                                   negativeEstimations);
    }else{
        Patch::curvatureEstimation(rangePositiveCurve.begin(),
                                   rangePositiveCurve.end(),
                                   positiveEstimations);

        Patch::curvatureEstimation(rangeNegativeCurve.begin(),
                                   rangeNegativeCurve.end(),
                                   negativeEstimations);
    }

    if(Development::solveShift){
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

    Patch::curvatureEstimation(rangePositiveCurve.begin(),
                               rangePositiveCurve.end(),
                               positiveEstimations);

    Patch::curvatureEstimation(rangeNegativeCurve.begin(),
                               rangeNegativeCurve.end(),
                               negativeEstimations);

    if(Development::solveShift){
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

        Patch::tangentEstimation(rangePositiveCurve.c(),
                                 rangePositiveCurve.c(),
                                 positiveEstimations);

        Patch::tangentEstimation(rangeNegativeCurve.c(),
                                 rangeNegativeCurve.c(),
                                 negativeEstimations);
    }else{
        Patch::tangentEstimation(rangePositiveCurve.begin(),
                                 rangePositiveCurve.end(),
                                 positiveEstimations);

        Patch::tangentEstimation(rangeNegativeCurve.begin(),
                                 rangeNegativeCurve.end(),
                                 negativeEstimations);
    }

    if(Development::solveShift){
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
    Patch::tangentEstimation(rangePositiveCurve.begin(),
                             rangePositiveCurve.end(),
                             positiveEstimations);

    std::vector<UtilsTypes::TangentVector> negativeEstimations;
    Patch::tangentEstimation(rangeNegativeCurve.begin(),
                             rangeNegativeCurve.end(),
                             negativeEstimations);

    if(Development::solveShift){
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

    ImageProc::computeBoundaryCurve(image,intCurve,0);

    UtilsTypes::Image2D dilatedImage(image.domain());
    ImageProc::dilate(dilatedImage,imgFilePath,1);
    ImageProc::computeBoundaryCurve(dilatedImage,extCurve,0);
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

void randomizeCurve(std::vector<Z2i::SCell>::const_iterator cBegin,
                    std::vector<Z2i::SCell>::const_iterator cEnd,
                    int size,
                    std::vector<Z2i::SCell>& scellsRand)
{
    DGtal::Circulator<std::vector<Z2i::SCell>::const_iterator> C(cBegin,
                                                                 cBegin,
                                                                 cEnd);

    int inc = rand()%size;
    while(inc>0){
        ++C;
        --inc;
    }

    DGtal::Circulator<std::vector<Z2i::SCell>::const_iterator> it = C;
    do{
        scellsRand.push_back(*it);
        ++it;
    }while(it!=C);
}


void draw(std::map<Z2i::SCell,double>& weightMap,
          std::vector<Z2i::SCell>& highlightedArcs,
          UtilsTypes::Board2D& board,
          double cmin,
          double cmax)
{
    if(cmin==cmax) cmax = cmin + 0.0000001;

    UtilsTypes::GradientColorMap<double,
            UtilsTypes::ColorGradientPreset::CMAP_JET> cmap_jet(cmin,cmax);

    board << UtilsTypes::SetMode(highlightedArcs[0].className(), "Paving" );
    std::string specificStyle = highlightedArcs[0].className() + "/Paving";


    for(auto it = highlightedArcs.begin();it!=highlightedArcs.end();++it)
    {
        board << UtilsTypes::CustomStyle( specificStyle,
                                          new UtilsTypes::CustomColors(UtilsTypes::Color::Black,
                                                                       cmap_jet( weightMap[*it] )));
        board << *it;
    }
}

void eliminateLoops(Curve& curveOut,
                    KSpace& KImage,
                    Curve& curveIn)
{
    std::vector<KSpace::SCell> vectorOfSCell;

    std::map<KSpace::SCell,KSpace::SCell> appearanceTable;
    KSpace::SCell toReconnectSCell;

    for(auto it=curveIn.begin();it!=curveIn.end();++it)
    {
        KSpace::SCell pointel = KImage.sDirectIncident(*it, *KImage.sDirs(*it));
        if (appearanceTable.find(pointel) == appearanceTable.end()) {
            appearanceTable[pointel] = *it;
            vectorOfSCell.push_back(*it);
        } else {
            toReconnectSCell = appearanceTable[pointel];

            KSpace::SCell backOfVector = vectorOfSCell.back();
            while (backOfVector != toReconnectSCell) {
                vectorOfSCell.pop_back();
                backOfVector = vectorOfSCell.back();
            }
        }
    }

    curveOut.initFromSCellsVector(vectorOfSCell);
}

