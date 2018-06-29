#include "API.h"

namespace BestLinel
{

    void setCurveWeight(SegCut::SCellGluedCurveIterator begin,
                        SegCut::SCellGluedCurveIterator end,
                        KSpace& KImage,
                        WeightMap& wm)
    {
        std::vector<double> estimationsCurvature;
        curvatureEstimatorsGluedCurve(begin,
                                      end,
                                      KImage,
                                      estimationsCurvature);

        updateToSquared(estimationsCurvature.begin(),estimationsCurvature.end());


        std::vector<WeightSettingsTypes::TangentVector> estimationsTangent;
        tangentEstimatorsGluedCurve(begin,
                                    end,
                                    KImage,
                                    estimationsTangent);

        std::vector<double> lengthWeightVector;
        auto it = begin;
        int i =0;
        do{
            lengthWeightVector.push_back( 1.0/( fabs(estimationsTangent[i][0]) + fabs(estimationsTangent[i][1]) ) );
            ++it;
            ++i;
        }while(it!=end);


        {
            int i =0;
            for(auto it=begin;it!=end;++it){
              wm[*it] = estimationsCurvature[i]*lengthWeightVector[i];
//            wm[*it] += 0.01*estimationsTangent[i];
                ++i;
            }
        }

    }

    double gdet(GluedCurveSetRange::ConstIterator gcRange, WeightMap& wm, KSpace& KImage,int g)
    {
        WeightMap gluedCurveWeights;

        setCurveWeight(gcRange->first,
                       gcRange->second,
                       KImage,
                       gluedCurveWeights);

        SegCut::SCellGluedCurveIterator it = gcRange->first;
        double gdetV=0;
        for(;it!=gcRange->second;++it)
        {
            gdetV+=wm[*it] - gluedCurveWeights[*it];
        }

        return gdetV;
    }

    void gneighborhood(GluedCurveSetRange::ConstIterator gcRange)
    {

    }
}