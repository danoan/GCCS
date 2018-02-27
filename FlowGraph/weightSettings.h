#ifndef SEGBYCUT_WEIGHTSETTINGS_H
#define SEGBYCUT_WEIGHTSETTINGS_H

#include <DGtal/io/boards/Board2D.h>
#include <boost/filesystem/path.hpp>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "../utils.h"
#include "ImageFlowData.h"

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
                        std::map<WeightSettingsTypes::Z2i::SCell,double>& weightMap);

void setGluedCurveWeight(WeightSettingsTypes::GluedCurveSetRange::ConstIterator gcsRangeBegin,
                         WeightSettingsTypes::GluedCurveSetRange::ConstIterator gcsRangeEnd,
                         KSpace& KImage,
                         unsigned int gluedCurveLength,
                         std::map<Z2i::SCell,double>& weightMap);

void setArcsWeight(ImageFlowData& imageFlowData,std::map<Z2i::SCell,double>& weightMap);

#endif //SEGBYCUT_WEIGHTSETTINGS_H
