#include "GluedCurveIterator.h"

template <typename CurveCirculator>
GluedCurveIterator<CurveCirculator>::GluedCurveIterator()
    :myIt1b(),
     myIt1e(),
     myIt2b(),
     myIt2e(),
     currentIterator(),
     cType(),
     myLinkSurfel(),
     myFlagIsValid(false),
     iteratorStage(0)
{
    element = new Z2i::SCell();
}

template <typename CurveCirculator>
GluedCurveIterator<CurveCirculator>::GluedCurveIterator(Z2i::SCell linkSurfel,
                                           CurveCirculator it1b,
                                           CurveCirculator it1e,
                                           CurveCirculator it2b,
                                           CurveCirculator it2e,
                                           ConnectorType cType,
                                           bool theEnd):

            myIt1b(it1b),
            myIt1e(it1e),
            myIt2b(it2b),
            myIt2e(it2e),
            cType(cType),
            myLinkSurfel(linkSurfel),
            myFlagIsValid(true)
{
    element = new Z2i::SCell();

    if(theEnd){
        iteratorStage = 2;
        currentIterator = it2e;
    }else{
        iteratorStage = 0;
        currentIterator = it1b;
    }

    myFlagIsValid = true;

}

template <typename CurveCirculator>
GluedCurveIterator<CurveCirculator>::GluedCurveIterator(const GluedCurveIterator& other):
            myIt1b(other.myIt1b),
            myIt1e(other.myIt1e),
            myIt2b(other.myIt2b),
            myIt2e(other.myIt2e),
            myLinkSurfel(other.myLinkSurfel),
            cType(other.cType),
            myFlagIsValid(other.myFlagIsValid),
            iteratorStage(other.iteratorStage),
            currentIterator(other.currentIterator)
{
    element = new Z2i::SCell();
}

template <typename CurveCirculator>
GluedCurveIterator<CurveCirculator>& GluedCurveIterator<CurveCirculator>::operator =(const GluedCurveIterator& other){
    myIt1b = other.myIt1b;
    myIt1e = other.myIt1e;
    myIt2b = other.myIt2b;
    myIt2e = other.myIt2e;

    cType = other.cType;

    myLinkSurfel = other.myLinkSurfel;
    myFlagIsValid = other.myFlagIsValid;
    iteratorStage = other.iteratorStage;
    currentIterator = other.currentIterator;

    element = new Z2i::SCell();

    return *this;
}

template <typename CurveCirculator>
void GluedCurveIterator<CurveCirculator>::increment()
{
    if( iteratorStage==0 || iteratorStage==2 ){
        if( iteratorStage==0 ) {
            if (currentIterator == myIt1e) {
                iteratorStage = 1;
                return;
            }
        }else{
            if (currentIterator == myIt2e) {
//                iteratorStage=0;
//                currentIterator=myIt1b;
                return;
            }
        }
        (currentIterator)++;
    }else {
        currentIterator = myIt2b;
        iteratorStage = 2;
    }

}

template <typename CurveCirculator>
bool GluedCurveIterator<CurveCirculator>::equal(const GluedCurveIterator& other) const
{
    if( iteratorStage == other.iteratorStage ){
        if(iteratorStage==1){
            return myLinkSurfel==other.myLinkSurfel;
        }else{
            return currentIterator==other.currentIterator;
        }

    }else{
        return false;
    }
}


template <typename CurveCirculator>
typename CurveCirculator::value_type& GluedCurveIterator<CurveCirculator>::dereference() const
{
    *element = iteratorStage==1?myLinkSurfel:*currentIterator;
    return *element;
}

template <typename CurveCirculator>
void GluedCurveIterator<CurveCirculator>::decrement()
{
    if( iteratorStage==0 || iteratorStage==2 ){
        if( iteratorStage==2 ) {
            if ( currentIterator == myIt2b) {
                iteratorStage = 1;
                return;
            }
        }else{
            if ( currentIterator == myIt1b) {
//                iteratorStage = 2;
//                currentIterator = myIt2e;
                return;
            }
        }
        currentIterator--;
    }else {
        currentIterator = myIt1e;
        iteratorStage = 0;
    }

}

