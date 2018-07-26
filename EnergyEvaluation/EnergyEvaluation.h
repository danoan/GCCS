#ifndef SEGBYCUT_ENERGYEVALUATION_H
#define SEGBYCUT_ENERGYEVALUATION_H

#include <DGtal/helpers/StdDefs.h>
#include "DGtal/topology/helpers/Surfaces.h"
#include "ImageProc/ImageProc.h"

#include "utils.h"

namespace EnergyEvaluation
{
    typedef DGtal::Z2i::DigitalSet DigitalSet;
    typedef DGtal::Z2i::Domain Domain;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

    typedef DGtal::Surfaces<KSpace> Surfaces;
    typedef DGtal::SurfelAdjacency<KSpace::dimension> SurfelAdjacency;

    typedef std::map<SCell,double> WeightMap;

    typedef ImageProc::DigitalBoundary<ImageProc::EightNeighborhoodPredicate<DigitalSet>> MyBoundary;

    double integratedSquaredCurvature(DigitalSet& ds);
}

#endif //SEGBYCUT_ENERGYEVALUATION_H
