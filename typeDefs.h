#ifndef SEGBYCUT_TYPEDEFS_H
#define SEGBYCUT_TYPEDEFS_H

#include "DGtal/helpers/StdDefs.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"

#include "GluedCurve/ConnectorSeedRange.h"
#include "GluedCurve/ConnectorFunctor.h"


namespace SegCut{
    /*Morphology Operations and Boundary Computation*/
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
    typedef DGtal::functors::SimpleThresholdForegroundPredicate<Image2D> ThreshPredicate;


    /*Curves Iterator*/
    typedef Curve::ConstIterator SCellIteratorType;
    typedef DGtal::Circulator<Curve::ConstIterator> SCellCirculatorType;
    typedef DGtal::Circulator<SCellIteratorType> SCellCirculator;


    typedef DGtal::PointVector<2,double> TangentVector;


    /*Glued Curves Iterators*/
    typedef ConnectorSeedRange<KSpace,SCellCirculator> ConnectorSeedRangeType;
    typedef ConnectorSeedToGluedCurveRange< ConnectorSeedRangeType::ConnectorSeedType > SeedToGluedCurveRangeFunctor;
    typedef SeedToGluedCurveRangeFunctor::GCIterator SCellGluedCurveIterator;

    typedef ConstRangeAdapter<  ConnectorSeedRangeType::ConnectorSeedIteratorType,
                                SeedToGluedCurveRangeFunctor,
                                SeedToGluedCurveRangeFunctor::Output > GluedCurveSetRange;

    typedef GluedCurveSetRange::ConstIterator GluedCurveIteratorPair;

    typedef ConstRangeAdapter<  SCellGluedCurveIterator,
                                DGtal::functors::SCellToIncidentPoints<KSpace>,
                                std::pair< DGtal::GridCurve<DGtal::Z2i::KSpace>::Point,
                                           DGtal::GridCurve<DGtal::Z2i::KSpace>::Point >
                             > GluedCurveIncidentPointsRange;
};

#endif //SEGBYCUT_TYPEDEFS_H
