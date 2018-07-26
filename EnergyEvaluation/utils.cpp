#include "utils.h"

namespace EnergyEvaluation
{
    void tangentWeight(Curve::ConstIterator begin,
                       Curve::ConstIterator end,
                       KSpace &KImage,
                       std::vector<TangentVector> &estimationsTangent,
                       std::vector<double> &tangentWeightVector)
    {
        auto it = begin;
        int i = 0;
        do {

            //Dot product estimator (projection) Suboptimal property

            Point pTarget = KImage.sCoords(
                    KImage.sDirectIncident(*it, *KImage.sDirs(*it)));
           Point pSource = KImage.sCoords(
                    KImage.sIndirectIncident(*it, *KImage.sDirs(*it)));

            Point scellVector = pTarget - pSource;

            tangentWeightVector.push_back(fabs(estimationsTangent[i].dot(scellVector)));


            ++it;
            ++i;
        } while (it != end);
    }

    void setGridCurveWeight(Curve curvePriorGS,
                            KSpace &KImage,
                            WeightMap &weightMap,
                            double flength)
    {
        std::vector<double> curvatureEstimations;
        curvatureEstimatorsGridCurve(curvePriorGS.begin(),
                                     curvePriorGS.end(),
                                     KImage,
                                     curvatureEstimations);


        updateToSquared(curvatureEstimations.begin(), curvatureEstimations.end());


        std::vector<TangentVector> tangentEstimations;
        tangentEstimatorsGridCurve(curvePriorGS.begin(),
                                   curvePriorGS.end(),
                                   KImage,
                                   tangentEstimations);


        std::vector<double> tangentWeightVector;
        tangentWeight(curvePriorGS.begin(),
                      curvePriorGS.end(),
                      KImage,
                      tangentEstimations,
                      tangentWeightVector);


        {
            int i = 0;
            for (auto it = curvePriorGS.begin(); it != curvePriorGS.end(); ++it) {
                weightMap[*it] = curvatureEstimations[i];
                ++i;
            }
        }

        {
            int i = 0;
            for (auto it = curvePriorGS.begin(); it != curvePriorGS.end(); ++it) {
                weightMap[*it] *= tangentWeightVector[i] * flength;
                //weightMap[*it] += 0.001*tangentWeightVector[i];
                ++i;
            }
        }
    }

    void curvatureEstimatorsGridCurve(Curve::ConstIterator begin,
                                      Curve::ConstIterator end,
                                      KSpace& KImage,
                                      std::vector<double>& estimations,
                                      bool closedCurve)
    {
        typedef DGtal::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;
        typedef DGtal::ConstRangeAdapter< Curve::ConstIterator,
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

    void tangentEstimatorsGridCurve(Curve::ConstIterator begin,
                                    Curve::ConstIterator end,
                                    KSpace& KImage,
                                    std::vector< TangentVector >& estimationsTangent,
                                    bool closedCurve)
    {
        typedef DGtal::ConstRangeAdapter< Curve::ConstIterator,
        DGtal::functors::SCellToPoint<KSpace>,
        KSpace::Point > RangeAdapter;

        DGtal::functors::SCellToPoint<KSpace> myFunctor = DGtal::functors::SCellToPoint<KSpace>(KImage);
        RangeAdapter rangePositiveCurve(begin,
                                        end,
                                        myFunctor);

        Curve negativeCurve;
        invertCurve(KImage,begin,end,negativeCurve);

        RangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                        negativeCurve.end(),
                                        myFunctor);

        std::vector<TangentVector> positiveEstimations;
        std::vector<TangentVector> negativeEstimations;
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
}