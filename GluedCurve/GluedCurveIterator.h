#ifndef FLOWGRAPH_GLUEDCURVEITERATOR_H
#define FLOWGRAPH_GLUEDCURVEITERATOR_H

#include "DGtal/helpers/StdDefs.h"


using namespace DGtal;
using namespace DGtal::Z2i;

template<typename CurveCirculator,typename LinkIteratorType>
class CurrentIteratorType{

public:
    CurrentIteratorType():cit(),
                          lit(),
                          gridIterator(false){};

    CurrentIteratorType(const CurrentIteratorType& other):cit(other.cit),
                                                          lit(other.lit),
                                                          gridIterator(other.gridIterator){};

    CurrentIteratorType& operator=(const CurrentIteratorType& other){
        cit = other.cit;
        lit = other.lit;
        gridIterator = other.gridIterator;
        return *this;
    };

    bool operator==(const CurveCirculator& c) const{
        return gridIterator && c==cit;
    }

    bool operator==(const LinkIteratorType& l) const{
        return !gridIterator && l==lit;
    }

    bool operator==(const CurrentIteratorType& other) const{
        return gridIterator==other.gridIterator && ( (gridIterator && cit==other.cit) || (!gridIterator && lit==other.lit) ) ;
    }

    void set(CurveCirculator& c){
        gridIterator = true;
        cit = c;
    }

    void set(LinkIteratorType& l){
        gridIterator = false;
        lit = l;
    }

    void operator++(){
        if(gridIterator) ++cit;
        else ++lit;
    }

    void operator--(){
        if(gridIterator) --cit;
        else --lit;
    }

    SCell operator*() const{
        if(gridIterator) return *cit;
        else return *lit;
    }

private:
    CurveCirculator cit;
    LinkIteratorType lit;
    bool gridIterator;

};

template <typename CurveCirculator, typename LinkIteratorType>
class GluedCurveIterator
        : public boost::iterator_facade<
                GluedCurveIterator<CurveCirculator,LinkIteratorType>,
                typename CurveCirculator::value_type,
                boost::bidirectional_traversal_tag
        >
{
private:
    CurveCirculator myIt1b,myIt1e,myIt2b,myIt2e;
    LinkIteratorType myItLb,myItLe;

    CurrentIteratorType<CurveCirculator,LinkIteratorType> currentIterator;

    ConnectorType cType;

    Z2i::SCell *element;

    char iteratorStage;

    bool myFlagIsValid;

public:

    GluedCurveIterator();

    GluedCurveIterator(  LinkIteratorType itLb,
                         LinkIteratorType itLe,
                         CurveCirculator it1b,
                         CurveCirculator it1e,
                         CurveCirculator it2b,
                         CurveCirculator it2e,
                         ConnectorType cType,
                         bool theEnd=false);

    GluedCurveIterator(const GluedCurveIterator& other);
    GluedCurveIterator<CurveCirculator,LinkIteratorType>& operator =(const GluedCurveIterator& other);

    ~GluedCurveIterator(){delete element;};

    inline Z2i::SCell linkSurfel(){return *myItLb;};
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
