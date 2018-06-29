#ifndef SEGBYCUT_API_H
#define SEGBYCUT_API_H

#include <DGtal/helpers/StdDefs.h>

#include "../FlowGraph/weightSettings.h"
#include <library.h>
#include "../utils/typeDefs.h"

namespace BestLinel
{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef std::map<Z2i::SCell,double> WeightMap;

    typedef SegCut::GluedCurveSetRange GluedCurveSetRange;

    //Internal use only
    void setCurveWeight(SegCut::SCellGluedCurveIterator begin,
                        SegCut::SCellGluedCurveIterator end,
                        KSpace& KImage,
                        WeightMap& wm);

    /*
     * Return the weight difference of linels in glued curves and
     * linels in original curve. A positive value means an improvement
     * in the energy.
     */
    double gdet(GluedCurveSetRange::ConstIterator gcRange, WeightMap& wm, KSpace& KImage,int g);


    void gneighborhood(GluedCurveSetRange::ConstIterator gcRange);
}


#endif //SEGBYCUT_API_H
