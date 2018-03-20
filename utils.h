#ifndef SEGBYCUT_CURVATUREWEIGHTMETHODS_H
#define SEGBYCUT_CURVATUREWEIGHTMETHODS_H


#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"

#include "typeDefs.h"
#include "imageProc.h"
#include "Patch/patch.h"

namespace Development{
    extern bool solveShift;
}

namespace UtilsTypes{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    using namespace SegCut;

    extern std::function< double(double) > toDouble;
}

void drawCurvatureMap(std::vector<Z2i::SCell>::const_iterator begin,
                      std::vector<Z2i::SCell>::const_iterator end,
                      double& cmin,
                      double& cmax,
                      Board2D& board,
                      std::map<Z2i::SCell,double>& weightMap);

void curvatureEstimatorsGridCurve(UtilsTypes::Curve::ConstIterator begin,
                                  UtilsTypes::Curve::ConstIterator end,
                                  UtilsTypes::KSpace& KImage,
                                  std::vector<double>& estimations,
                                  bool closedCurve=true);

void curvatureEstimatorsGluedCurve(UtilsTypes::SCellGluedCurveIterator begin,
                                   UtilsTypes::SCellGluedCurveIterator end,
                                   UtilsTypes::KSpace& KImage,
                                   std::vector<double>& estimations);

void curvatureEstimatorsConnections(UtilsTypes::GluedCurveIteratorPair begin,
                                    UtilsTypes::GluedCurveIteratorPair end,
                                    UtilsTypes::KSpace& KImage,
                                    unsigned int gluedCurveLength,
                                    std::vector<double>& estimations);

void tangentEstimatorsGridCurve(UtilsTypes::Curve::ConstIterator begin,
                                UtilsTypes::Curve::ConstIterator end,
                                UtilsTypes::KSpace& KImage,
                                std::vector< UtilsTypes::TangentVector >& estimationsTangent,
                                bool closedCurve=true);

void tangentEstimatorsGluedCurve(UtilsTypes::SCellGluedCurveIterator begin,
                                 UtilsTypes::SCellGluedCurveIterator end,
                                 UtilsTypes::KSpace& KImage,
                                 std::vector< UtilsTypes::TangentVector >& estimationsTangent );

void tangentEstimatorsConnections(UtilsTypes::GluedCurveIteratorPair begin,
                                  UtilsTypes::GluedCurveIteratorPair end,
                                  UtilsTypes::KSpace& KImage,
                                  unsigned int gluedCurveLength,
                                  std::vector< UtilsTypes::TangentVector >& estimations);


void setKImage(std::string imgFilePath,KSpace& KImage);

void setCurves(std::string imgFilePath,
               UtilsTypes::Curve& intCurve,
               UtilsTypes::Curve& extCurve);

UtilsTypes::ConnectorSeedRangeType getSeedRange(UtilsTypes::KSpace& KImage,
                                                UtilsTypes::Curve& intCurve,
                                                UtilsTypes::Curve& extCurve);


void randomizeCurve(std::vector<Z2i::SCell>::const_iterator cBegin,
                    std::vector<Z2i::SCell>::const_iterator cEnd,
                    int size,
                    std::vector<Z2i::SCell>& scellsRand);


template<typename SCellIteratorType>
void invertCurve(KSpace& KImage,
                 SCellIteratorType begin,
                 SCellIteratorType end,
                 UtilsTypes::Curve& c2)
{
    std::vector<UtilsTypes::Z2i::SCell> SCells;
    auto it=begin;
    do{
        SCells.push_back(*it);
        ++it;
    }while(it!=end);

    std::vector<UtilsTypes::Z2i::SCell> newSCells;
    {
        auto it = SCells.rbegin();
        do{
            UtilsTypes::Z2i::SCell newLinel = KImage.sCell( *it);
            KImage.sSetSign(newLinel,!KImage.sSign(*it));

            newSCells.push_back(newLinel);
            ++it;;
        }while(it!=SCells.rend());
    }

    c2.initFromSCellsVector(newSCells);
}

template<typename Value>
void setWeight(UtilsTypes::Curve::ConstIterator begin,
               UtilsTypes::Curve::ConstIterator end,
               UtilsTypes::KSpace& KImage,
               std::map<UtilsTypes::SCell,Value> &weightMap,
               std::vector<Value>& addValues)
{
    auto it = begin;
    int i = 0;
    do {
        weightMap[*it] += addValues[i];
        ++it;
        ++i;
    } while (it != end);
}

template<typename Value>
void max_and_min(std::vector<Value>& V,
                 double& cmin,
                 double& cmax,
                 typename std::function< double(Value) >& toDouble = UtilsTypes::toDouble
)
{
    for(auto itV=V.begin();itV!=V.end();itV++){
        if( toDouble(*itV) < cmin ) cmin = toDouble(*itV);
        if( toDouble(*itV) > cmax ) cmax = toDouble(*itV);
    }

}

template<typename IteratorType>
void updateToSquared(IteratorType begin, IteratorType end)
{
    IteratorType it = begin;
    do{
        *it = pow(*it,2);
        ++it;
    }while(it!=end);
}

template<typename Value, typename IteratorType>
void draw(std::vector<Value>& V,
          IteratorType begin,
          IteratorType end,
          UtilsTypes::Board2D& board,
          double cmin,
          double cmax,
          typename std::function< double(Value) >& toDouble = UtilsTypes::toDouble)
{
    if(cmin==cmax) cmax = cmin + 0.0000001;

    UtilsTypes::GradientColorMap<double,
                                 UtilsTypes::ColorGradientPreset::CMAP_JET> cmap_jet(cmin,cmax);

    board << UtilsTypes::SetMode( begin->className(), "Paving" );
    std::string specificStyle = begin->className() + "/Paving";

    {
        auto itSCell = begin;
        auto itV = V.begin();
        do{
            board << UtilsTypes::CustomStyle( specificStyle,
                                              new UtilsTypes::CustomColors(UtilsTypes::Color::Black,
                                                                           cmap_jet( toDouble(*itV) )));
            board << *itSCell;

            ++itSCell;
            ++itV;
        }while(itSCell!=end);
    }
}


void draw(std::map<Z2i::SCell,double>& weightMap,
          std::vector<Z2i::SCell>& highlightedArcs,
          UtilsTypes::Board2D& board,
          double cmin,
          double cmax);

void eliminateLoops(Curve& curveOut,
                    KSpace& KImage,
                    Curve& curveIn);

#endif //SEGBYCUT_CURVATUREWEIGHTMETHODS_H
