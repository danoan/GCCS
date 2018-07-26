#include "FlowGraphUtils.h"

void setArcsWeight(ImageFlowData& imageFlowData,std::map<Z2i::SCell,double>& weightMap)
{

    double flength = 1;
    double innerLength = computeLength(imageFlowData.getMostInnerCurve(),imageFlowData.getKSpace());
    double outerLength = computeLength(imageFlowData.getMostOuterCurve(),imageFlowData.getKSpace());
    
    flength = innerLength/outerLength;
    
    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        setGridCurveWeight(it->curve,
                           imageFlowData.getKSpace(),
                           weightMap,
                           flength
        );
    }

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        setGluedCurveWeight( it->gcsRangeBegin(),
                             it->gcsRangeEnd(),
                             imageFlowData.getKSpace(),
                             imageFlowData.getGluedCurveLength(),
                             weightMap,
                             flength);
    }
}