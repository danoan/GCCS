#include "weightSettings.h"

void tangentWeight(WeightSettingsTypes::Curve::ConstIterator begin,
                   WeightSettingsTypes::Curve::ConstIterator end,
                   WeightSettingsTypes::KSpace& KImage,
                   std::vector< WeightSettingsTypes::TangentVector >& estimationsTangent,
                   std::vector< double >& tangentWeightVector)
{
    auto it = begin;
    int i =0;
    do{
        WeightSettingsTypes::KSpace::Point pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
        WeightSettingsTypes::KSpace::Point pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );

        WeightSettingsTypes::KSpace::Point scellVector = pTarget-pSource;

        tangentWeightVector.push_back( fabs( estimationsTangent[i].dot(scellVector) ) );
        ++it;
        ++i;
    }while(it!=end);
}

void setGridCurveWeight(Curve curvePriorGS,
                        KSpace& KImage,
                        std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> curvatureEstimations;
    curvatureEstimatorsGridCurve(curvePriorGS.begin(),
                                 curvePriorGS.end(),
                                 KImage,
                                 curvatureEstimations);


    updateToSquared(curvatureEstimations.begin(),curvatureEstimations.end());


    std::vector<WeightSettingsTypes::TangentVector> tangentEstimations;
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
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] = curvatureEstimations[i];
            ++i;
        }
    }

    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] *=tangentWeightVector[i];
            weightMap[*it] += 0.01*tangentWeightVector[i];
            ++i;
        }
    }

}

void setGluedCurveWeight(WeightSettingsTypes::GluedCurveSetRange::ConstIterator gcsRangeBegin,
                         WeightSettingsTypes::GluedCurveSetRange::ConstIterator gcsRangeEnd,
                         KSpace& KImage,
                         unsigned int gluedCurveLength,
                         std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> estimationsCurvature;
    curvatureEstimatorsConnections(gcsRangeBegin,gcsRangeEnd,KImage,gluedCurveLength,estimationsCurvature);

    updateToSquared(estimationsCurvature.begin(),estimationsCurvature.end());

    {
        int i = 0;
        for (WeightSettingsTypes::GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {

            auto itC = it->first.connectorsBegin();
            do{
                weightMap[*itC] = estimationsCurvature[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);

        }
    }

    std::vector<WeightSettingsTypes::TangentVector> estimationsTangent;
    tangentEstimatorsConnections(gcsRangeBegin,gcsRangeEnd,KImage,gluedCurveLength,estimationsTangent);


    std::vector<double> tangentWeightVector;
    {
        WeightSettingsTypes::GluedCurveIteratorPair it = gcsRangeBegin;
        int i = 0;
        for (WeightSettingsTypes::GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {

            auto itC = it->first.connectorsBegin();
            do {

                KSpace::SCell linel = *itC;

                UtilsTypes::KSpace::Point pTarget = KImage.sCoords(KImage.sDirectIncident(linel, *KImage.sDirs(linel)));
                UtilsTypes::KSpace::Point pSource = KImage.sCoords(KImage.sIndirectIncident(linel, *KImage.sDirs(linel)));

                UtilsTypes::KSpace::Point scellVector = pTarget - pSource;

                tangentWeightVector.push_back(fabs(estimationsTangent[i].dot(scellVector)));

                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);

        }

    }

    {
        int i = 0;
        for (WeightSettingsTypes::GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {
            auto itC = it->first.connectorsBegin();
            do {
                weightMap[*itC]*= tangentWeightVector[i];
                weightMap[*itC]+= 0.01*tangentWeightVector[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);
        }
    }

}