#ifndef FLOWGRAPH_GLUEDCURVEITERATOR_H
#define FLOWGRAPH_GLUEDCURVEITERATOR_H

#include "DGtal/helpers/StdDefs.h"


using namespace DGtal;
using namespace DGtal::Z2i;


template <typename CurveCirculator>
class GluedCurveIterator
        : public boost::iterator_facade<
                GluedCurveIterator<CurveCirculator>,
                typename CurveCirculator::value_type,
                boost::bidirectional_traversal_tag
        >
{
private:
    CurveCirculator myIt1b,myIt1e,myIt2b,myIt2e;
    CurveCirculator currentIterator;

    ConnectorType cType;

    Z2i::SCell myLinkSurfel;
    Z2i::SCell *element;

    char iteratorStage;

    bool myFlagIsValid;

public:

    GluedCurveIterator();

    GluedCurveIterator(Z2i::SCell linkSurfel,
                         CurveCirculator it1b,
                         CurveCirculator it1e,
                         CurveCirculator it2b,
                         CurveCirculator it2e,
                         ConnectorType cType,
                         bool theEnd=false);

    GluedCurveIterator(const GluedCurveIterator& other);
    GluedCurveIterator<CurveCirculator>& operator =(const GluedCurveIterator& other);

    ~GluedCurveIterator(){delete element;};

    inline Z2i::SCell linkSurfel(){return myLinkSurfel;};
    inline ConnectorType connectorType(){return cType;};

private:
    friend class boost::iterator_core_access;

    void increment();

    bool equal(const GluedCurveIterator& other) const;

    typename CurveCirculator::value_type& dereference() const;

    void decrement();
};

#include "GluedCurveIterator.cpp"


#endif //FLOWGRAPH_GLUEDCURVEITERATOR_H
