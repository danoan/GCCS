#ifndef SEGBYCUT_WEIGHTSETTINGS_H
#define SEGBYCUT_WEIGHTSETTINGS_H

#include <DGtal/io/boards/Board2D.h>
#include <boost/filesystem/path.hpp>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "../utils/utils.h"

namespace WeightSettingsTypes
{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    using namespace UtilsTypes;
}

void tangentWeight(WeightSettingsTypes::Curve::ConstIterator begin,
                   WeightSettingsTypes::Curve::ConstIterator end,
                   WeightSettingsTypes::KSpace& KImage,
                   std::vector< WeightSettingsTypes::TangentVector >& estimationsTangent,
                   std::vector< double >& tangentWeightVector);

void setGridCurveWeight(WeightSettingsTypes::Curve curvePriorGS,
                        WeightSettingsTypes::KSpace& KImage,
                        std::map<WeightSettingsTypes::Z2i::SCell,double>& weightMap,double flength=1);

void setGluedCurveWeight(WeightSettingsTypes::GluedCurveSetRange::ConstIterator gcsRangeBegin,
                         WeightSettingsTypes::GluedCurveSetRange::ConstIterator gcsRangeEnd,
                         KSpace& KImage,
                         unsigned int gluedCurveLength,
                         std::map<Z2i::SCell,double>& weightMap,double flength=1);

double computeLength(Curve& curve, KSpace& KImage);


#endif //SEGBYCUT_WEIGHTSETTINGS_H
