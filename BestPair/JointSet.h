#ifndef SEGBYCUT_JOINTSET_H
#define SEGBYCUT_JOINTSET_H

#include "../FlowGraph/ImageFlowData.h"

namespace BestLinel
{
    class JointSet
    {
    public:
        typedef DGtal::Z2i::SCell SCell;
        typedef DGtal::Z2i::Curve::const_iterator SCellIterator;
        typedef DGtal::Circulator<DGtal::Z2i::Curve::const_iterator> CurveCirculator;
        typedef ImageFlowData::SCellGluedCurveIterator GluedCurveIterator;

        struct Joint
        {
            Joint(SCell s,
                  CurveCirculator antCurve,
                  CurveCirculator sucCurve,
                  int index):scell(s),
                             antCurve(antCurve),
                             sucCurve(sucCurve),
                             index(index){}

            SCell scell;
            CurveCirculator antCurve;
            CurveCirculator sucCurve;
            int index;
        };

        typedef std::vector<Joint>::const_iterator JointIterator;


    public:
        JointSet(ImageFlowData &ifd);

        std::vector<Joint> inJointsVector;
        std::vector<Joint> outJointsVector;

    };
}
#endif //SEGBYCUT_JOINTSET_H
