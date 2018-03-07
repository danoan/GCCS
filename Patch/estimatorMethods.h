#ifndef SEGBYCUT_ESTIMATORMETHODS_H
#define SEGBYCUT_ESTIMATORMETHODS_H

#include "DGtal/helpers/StdDefs.h"

#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/estimation/LambdaMST2D.h"

#include "DGtal/geometry/curves/StabbingCircleComputer.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"

#include "MostCenteredMaximalSegmentEstimatorWrapper.h"
#include "LambdaMaximalSegmentEstimator.h"
#include "PessimistMaximalSegmentEstimator.h"

namespace Development
{
    extern bool lambdaEstimator;
    extern bool pessimistEstimator;
}

namespace Patch{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef DGtal::PointVector<2,double> TangentVector;

    template<typename IteratorType>
    void curvatureEstimation(IteratorType itb,
                             IteratorType ite,
                             std::vector<double>& estimations);


    template<typename IteratorType>
    void tangentEstimation(IteratorType itb,
                           IteratorType ite,
                           std::vector< TangentVector >& estimations);

    template<typename IteratorType>
    void estimationsPessimistMDCACurvature(IteratorType itb,
                                           IteratorType ite,
                                           std::vector<double>& estimations);

    template<typename IteratorType>
    void estimationsLambdaMDCACurvature(IteratorType itb,
                                        IteratorType ite,
                                        std::vector<double>& estimations);

    template<typename IteratorType>
    void estimationsMCMDCACurvature(IteratorType itb,
                                    IteratorType ite,
                                    std::vector<double>& estimations);


    template<typename IteratorType>
    void estimationsLambdaMDSSTangent(IteratorType itb,
                                      IteratorType ite,
                                      std::vector<TangentVector>& estimations);

    template<typename IteratorType>
    void estimationsMCMDSSTangent(IteratorType itb,
                                  IteratorType ite,
                                  std::vector<TangentVector>& estimations);
}

#include "estimatorMethods.ih"

#endif //SEGBYCUT_ESTIMATORMETHODS_H
