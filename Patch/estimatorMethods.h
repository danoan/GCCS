#ifndef SEGBYCUT_ESTIMATORMETHODS_H
#define SEGBYCUT_ESTIMATORMETHODS_H

#include "DGtal/helpers/StdDefs.h"

#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/StabbingCircleComputer.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"

#include "MostCenteredMaximalSegmentEstimatorWrapper.h"


namespace Patch{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef DGtal::PointVector<2,double> TangentVector;

    namespace Backstage {

        namespace Functors{
            template<typename CurveSCellIteratorType>
            class IsGoodStarterStabbingCircles{
            public:
                bool operator()(const CurveSCellIteratorType& itb,
                                const CurveSCellIteratorType& ite,
                                KSpace& KImage/*,
                                std::pair< KSpace::Point, KSpace::Point>& incidentPixels*/);
            };

            template<typename CurveSCellIteratorType>
            class IsGoodStarterArithmeticalDSSComputer{
            public:
                bool operator()(const CurveSCellIteratorType& itb,
                                const CurveSCellIteratorType& ite,
                                KSpace& KImage/*,
                                std::pair< KSpace::Point, KSpace::Point>& incidentPixels*/);
            };
        }

        template<typename FunctorGoodStarter>
        void initializeCurveWithAGoodStarter(KSpace &KImage,
                                             FunctorGoodStarter fgs,
                                             Curve &initialCurve,
                                             Curve &newCurve);

        template<typename IteratorType, typename FunctorGoodStarter>
        void findGoodStarter(KSpace &KImage,
                             FunctorGoodStarter isGoodStarter,
                             IteratorType &curveBegin,
                             IteratorType &curveEnd,
                             std::vector<SCell>& scells);


        template<typename IteratorType, typename SegmentComputer, typename Segmentation, typename SegmentIterator>
        bool IsGoodStarter(const IteratorType &itb,
                           const IteratorType &ite/*,
                           std::pair <KSpace::Point, KSpace::Point> &incidentPixels*/);

    }

    /*Closed Curves Only*/
    void initializeCurveCurvatureEstimator(KSpace &KImage,
                                           Curve &initialCurve,
                                           Curve &newCurve);

    /*Closed Curves Only*/
    void initializeCurveTangentEstimator(KSpace &KImage,
                                         Curve &initialCurve,
                                         Curve &newCurve);


    template<typename IteratorType>
    void estimationsPatchMCMSECurvature(IteratorType itb,
                                        IteratorType ite,
                                        std::vector<double>& estimations);

    template<typename IteratorType>
    void estimationsDGtalMCMSECurvature(IteratorType itb,
                                        IteratorType ite,
                                        std::vector<double>& estimations);


    template<typename IteratorType>
    void estimationsPatchMCMSETangent(IteratorType itb,
                                      IteratorType ite,
                                      std::vector< TangentVector >& estimations);

    template<typename IteratorType>
    void estimationsDGtalMCMSETangent(IteratorType itb,
                                      IteratorType ite,
                                      std::vector< TangentVector >& estimations);
}

#include "estimatorMethods.ih"

#endif //SEGBYCUT_ESTIMATORMETHODS_H
