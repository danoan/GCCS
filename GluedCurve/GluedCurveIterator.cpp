#include "GluedCurveIterator.h"

template <typename CurveCirculator,typename LinkIteratorType>
GluedCurveIterator<CurveCirculator,LinkIteratorType>::GluedCurveIterator()
    :myItLb(),
     myItLe(),
     myIt1b(),
     myIt1e(),
     myIt2b(),
     myIt2e(),
     currentIterator(),
     cType(),
     myFlagIsValid(false),
     iteratorStage(0)
{
    element = new Z2i::SCell();
}

template <typename CurveCirculator,typename LinkIteratorType>
GluedCurveIterator<CurveCirculator,LinkIteratorType>::GluedCurveIterator(LinkIteratorType itLb,
                                                                         LinkIteratorType itLe,
                                                                         CurveCirculator it1b,
                                                                         CurveCirculator it1e,
                                                                         CurveCirculator it2b,
                                                                         CurveCirculator it2e,
                                                                         ConnectorType cType,
                                                                         bool theEnd):
            myItLb(itLb),
            myItLe(itLe),
            myIt1b(it1b),
            myIt1e(it1e),
            myIt2b(it2b),
            myIt2e(it2e),
            cType(cType),
            myFlagIsValid(true)
{
    element = new Z2i::SCell();

    if(theEnd){
        iteratorStage = 2;
        currentIterator.set(it2e);
    }else{
        iteratorStage = 0;
        currentIterator.set(it1b);
    }

    myFlagIsValid = true;

}

template <typename CurveCirculator,typename LinkIteratorType>
GluedCurveIterator<CurveCirculator,LinkIteratorType>::GluedCurveIterator(const GluedCurveIterator& other):
            myItLb(other.myItLb),
            myItLe(other.myItLe),
            myIt1b(other.myIt1b),
            myIt1e(other.myIt1e),
            myIt2b(other.myIt2b),
            myIt2e(other.myIt2e),
            cType(other.cType),
            myFlagIsValid(other.myFlagIsValid),
            iteratorStage(other.iteratorStage),
            currentIterator(other.currentIterator)
{
    element = new Z2i::SCell();
}

template <typename CurveCirculator,typename LinkIteratorType>
GluedCurveIterator<CurveCirculator,LinkIteratorType>& GluedCurveIterator<CurveCirculator,LinkIteratorType>::operator =(const GluedCurveIterator& other){
    myItLb = other.myItLb;
    myItLe = other.myItLe;

    myIt1b = other.myIt1b;
    myIt1e = other.myIt1e;
    myIt2b = other.myIt2b;
    myIt2e = other.myIt2e;

    cType = other.cType;

    myFlagIsValid = other.myFlagIsValid;
    iteratorStage = other.iteratorStage;
    currentIterator = other.currentIterator;

    element = new Z2i::SCell();

    return *this;
}

template <typename CurveCirculator,typename LinkIteratorType>
void GluedCurveIterator<CurveCirculator,LinkIteratorType>::increment()
{
    if( iteratorStage==0 || iteratorStage==2 ){
        if( iteratorStage==0 ) {
            if (currentIterator == myIt1e) {
                iteratorStage = 1;
                currentIterator.set(myItLb);
                return;
            }
        }else{
            if (currentIterator == myIt2e) {
//                iteratorStage=0;
//                currentIterator=myIt1b;
                return;
            }
        }
    }else {
        if(currentIterator == myItLe){
            currentIterator.set(myIt2b);
            iteratorStage = 2;
            return;
        }
    }
    ++(currentIterator);

}

template <typename CurveCirculator,typename LinkIteratorType>
bool GluedCurveIterator<CurveCirculator,LinkIteratorType>::equal(const GluedCurveIterator& other) const
{
    if( iteratorStage == other.iteratorStage ){
        return currentIterator==other.currentIterator;
    }else{
        return false;
    }
}


template <typename CurveCirculator,typename LinkIteratorType>
typename CurveCirculator::value_type& GluedCurveIterator<CurveCirculator,LinkIteratorType>::dereference() const
{
    *element = *currentIterator;
    return *element;
}

template <typename CurveCirculator,typename LinkIteratorType>
void GluedCurveIterator<CurveCirculator,LinkIteratorType>::decrement()
{
    if( iteratorStage==0 || iteratorStage==2 ){
        if( iteratorStage==2 ) {
            if ( currentIterator == myIt2b) {
                iteratorStage = 1;
                currentIterator.set(myItLe);
                return;
            }
        }else{
            if ( currentIterator == myIt1b) {
//                iteratorStage = 2;
//                currentIterator = myIt2e;
                return;
            }
        }
    }else {
        if(currentIterator == myItLb){
            currentIterator.set(myIt1e);
            iteratorStage = 0;
            return;
        }
    }
    --currentIterator;

}

