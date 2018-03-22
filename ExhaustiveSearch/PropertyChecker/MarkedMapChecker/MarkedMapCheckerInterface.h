#ifndef SEGBYCUT_CONSISTENCYCHECKER_H
#define SEGBYCUT_CONSISTENCYCHECKER_H

#include <DGtal/helpers/StdDefs.h>
#include <ConnectorSeedRange.h>
#include "../CheckableElementConcept.h"

template<typename CheckableElement>
class MarkedMapCheckerInterface
{
public:
    BOOST_CONCEPT_ASSERT( (CheckableElementConcept<CheckableElement>) );

    typedef typename CheckableElement::CheckedType CheckedType;
    typedef typename CheckableElement::MarkedType MarkedType;

public:
    virtual bool operator()(const CheckableElement& e)const{return true;};
    virtual void mark(const CheckableElement& e){};
    virtual void unmark(const CheckableElement& e){};
};




#endif //SEGBYCUT_CONSISTENCYCHECKER_H
