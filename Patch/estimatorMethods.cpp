#include "estimatorMethods.h"

namespace Patch{

    void initializeCurveCurvatureEstimator(KSpace &KImage,
                                           Curve &initialCurve,
                                           Curve &newCurve)
    {

        Backstage::Functors::IsGoodStarterStabbingCircles<Curve::ConstIterator> gsStabbingFunctor;

        Backstage::initializeCurveWithAGoodStarter(KImage,
                                                   gsStabbingFunctor,
                                                   initialCurve,
                                                   newCurve);
    };

    void initializeCurveTangentEstimator(KSpace &KImage,
                                         Curve &initialCurve,
                                         Curve &newCurve)
    {
        Backstage::Functors::IsGoodStarterArithmeticalDSSComputer<Curve::ConstIterator> gsArithmeticalFunctor;

        Backstage::initializeCurveWithAGoodStarter(KImage,
                                                   gsArithmeticalFunctor,
                                                   initialCurve,
                                                   newCurve);
    };

}