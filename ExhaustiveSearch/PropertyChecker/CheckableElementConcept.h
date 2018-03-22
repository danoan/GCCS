#ifndef SEGBYCUT_CHECKABLEELEMENT_H
#define SEGBYCUT_CHECKABLEELEMENT_H


#include <boost/concept/assert.hpp>
#include <ConnectorSeedRange.h>
#include <DGtal/helpers/StdDefs.h>

template<typename T>
class CheckableElementConcept
{
    typedef typename T::CheckedType CheckedType;
    typedef typename T::MarkedType MarkedType;
    typedef typename T::ComparisonClass ComparisonClass;
};


#endif //SEGBYCUT_CHECKABLEELEMENT_H
