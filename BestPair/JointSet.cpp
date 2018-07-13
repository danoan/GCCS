#include "JointSet.h"

namespace BestLinel
{
    JointSet::JointSet(ImageFlowData& ifd)
    {

        std::map<SCell,int> indexOrder;
        int index=0;
        for(auto it=ifd.getMostOuterCurve().begin();it!=ifd.getMostOuterCurve().end();++it,++index)
        {
            indexOrder[*it]=index;
        }

        for(auto it=ifd.curvePairBegin();it!=ifd.curvePairEnd();++it)
        {
            ImageFlowData::GluedCurveSetRange gcsRange = it->getGCSRange();
            for(auto iu=gcsRange.begin();iu!=gcsRange.end();++iu)
            {


                if(iu->first.connectorType()==ConnectorType::internToExtern)
                {
                    int jointIndex = indexOrder[ *iu->first.curveSegment2Begin() ];
                    inJointsVector.push_back(Joint(*iu->first.connectorsBegin(),
                                                   iu->first.curveSegment1End(),
                                                   iu->first.curveSegment2Begin(),
                                                   jointIndex));
                }
                else
                {
                    int jointIndex = indexOrder[ *iu->first.curveSegment1End() ];
                    outJointsVector.push_back( Joint(*iu->first.connectorsBegin(),
                                                     iu->first.curveSegment1End(),
                                                     iu->first.curveSegment2Begin(),
                                                     jointIndex ) );
                }


            }
        }

        std::sort(inJointsVector.begin(),inJointsVector.end(),[](const Joint& j1, const Joint& j2){ return j1.index < j2.index;});
        std::sort(outJointsVector.begin(),outJointsVector.end(),[](const Joint& j1, const Joint& j2){ return j1.index < j2.index;});
    }

}