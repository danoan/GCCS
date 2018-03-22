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

        //Dot product estimator (projection)

//        WeightSettingsTypes::KSpace::Point pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
//        WeightSettingsTypes::KSpace::Point pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );
//
//        WeightSettingsTypes::KSpace::Point scellVector = pTarget-pSource;
//
//        tangentWeightVector.push_back( fabs( estimationsTangent[i].dot(scellVector) ) );



        //1.0/(cos+sin) Length estimation
        tangentWeightVector.push_back( 1.0/( fabs(estimationsTangent[i][0]) + fabs(estimationsTangent[i][1]) ) );


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


    int totalLinels=0;
    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] = curvatureEstimations[i];
            ++i;
            ++totalLinels;
        }
    }

    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] *= tangentWeightVector[i];
//            weightMap[*it] *= 100;//+100 only to avoid small numbers
//            weightMap[*it] /= totalLinels;
            weightMap[*it] += 0.001*tangentWeightVector[i];
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

    int totalLinels = estimationsCurvature.size();
    {
        int i = 0;
        for (WeightSettingsTypes::GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {
            auto itC = it->first.connectorsBegin();
            do {
                weightMap[*itC]*= tangentWeightVector[i];
//                weightMap[*itC]*= 100; //+100 only to avoid small numbers
//                weightMap[*itC]/= totalLinels;
                weightMap[*itC]+= 0.1*tangentWeightVector[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);
        }
    }

}

void setArcsWeight(ImageFlowData& imageFlowData,std::map<Z2i::SCell,double>& weightMap)
{

    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        setGridCurveWeight(it->curve,
                           imageFlowData.getKSpace(),
                           weightMap
        );
    }

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        setGluedCurveWeight( it->gcsRangeBegin(),
                             it->gcsRangeEnd(),
                             imageFlowData.getKSpace(),
                             imageFlowData.getGluedCurveLength(),
                             weightMap);
    }
}
