#ifndef SEGBYCUT_ENERGYEVALUATION_UTILS_H
#define SEGBYCUT_ENERGYEVALUATION_UTILS_H

#include <DGtal/helpers/StdDefs.h>
#include "Patch/estimatorMethods.h"

namespace Development
{
    extern bool solveShift;
}

namespace EnergyEvaluation
{
    typedef DGtal::Z2i::DigitalSet DigitalSet;
    typedef DGtal::Z2i::Domain Domain;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

    typedef DGtal::PointVector<2,double> TangentVector;
    typedef std::map<SCell,double> WeightMap;


    void tangentWeight(Curve::ConstIterator begin,
                       Curve::ConstIterator end,
                       KSpace& KImage,
                       std::vector< TangentVector >& estimationsTangent,
                       std::vector< double >& tangentWeightVector);

    void setGridCurveWeight(Curve curvePriorGS,
                            KSpace& KImage,
                            WeightMap& weightMap,
                            double flength=1);

    void curvatureEstimatorsGridCurve(Curve::ConstIterator begin,
                                      Curve::ConstIterator end,
                                      KSpace& KImage,
                                      std::vector<double>& estimations,
                                      bool closedCurve=true);

    void tangentEstimatorsGridCurve(Curve::ConstIterator begin,
                                    Curve::ConstIterator end,
                                    KSpace& KImage,
                                    std::vector< TangentVector >& estimationsTangent,
                                    bool closedCurve=true);

    template<typename IteratorType>
    void updateToSquared(IteratorType begin, IteratorType end)
    {
        IteratorType it = begin;
        do{
            *it = pow(*it,2);
            ++it;
        }while(it!=end);
    }

    template<typename SCellIteratorType>
    void invertCurve(KSpace& KImage,
                     SCellIteratorType begin,
                     SCellIteratorType end,
                     Curve& c2)
    {
        std::vector<SCell> SCells;
        auto it=begin;
        do{
            SCells.push_back(*it);
            ++it;
        }while(it!=end);

        std::vector<SCell> newSCells;
        {
            auto it = SCells.rbegin();
            do{
                SCell newLinel = KImage.sCell( *it);
                KImage.sSetSign(newLinel,!KImage.sSign(*it));

                newSCells.push_back(newLinel);
                ++it;
            }while(it!=SCells.rend());
        }

        c2.initFromSCellsVector(newSCells);
    }

}
#endif