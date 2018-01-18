#ifndef SEGBYCUT_CURVATUREWEIGHTMETHODS_H
#define SEGBYCUT_CURVATUREWEIGHTMETHODS_H


#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"

#include "typeDefs.h"
#include "imageProc.h"
#include "Patch/patch.h"

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
                                  std::vector<double>& estimations);

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
                                std::vector< UtilsTypes::TangentVector >& estimationsTangent );

void tangentEstimatorsGluedCurve(UtilsTypes::Curve::ConstIterator begin,
                                 UtilsTypes::Curve::ConstIterator end,
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
    cmax = cmax>cmin?cmax:cmin+1;
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


#endif //SEGBYCUT_CURVATUREWEIGHTMETHODS_H
