#ifndef FLOWGRAPH_GLUEDCURVEITERATOR_H
#define FLOWGRAPH_GLUEDCURVEITERATOR_H

#include "DGtal/helpers/StdDefs.h"
#include "ConnectorSeed.h"

using namespace DGtal;
using namespace DGtal::Z2i;

template<typename CurveCirculator>
class CurrentIteratorType{
public:

    typedef std::vector<SCell>::const_iterator LinkIteratorType;

    CurrentIteratorType(){
        element = new SCell();
    }

    CurrentIteratorType(const CurrentIteratorType& other){
        grid = other.grid;
        gridIt = other.gridIt;
        connectorIt = other.connectorIt;

        element = new SCell();
    }

    CurrentIteratorType& operator=(const CurrentIteratorType& other){
        grid = other.grid;
        gridIt = other.gridIt;
        connectorIt = other.connectorIt;

        element = new SCell();

        return *this;
    }

    ~CurrentIteratorType(){delete element;}

    void set(CurveCirculator& it){
        grid=true;
        gridIt = it;
    }

    void set(LinkIteratorType& it){
        grid=false;
        connectorIt = it;
    }


    bool operator==(const CurrentIteratorType& other) const{
        return (grid && grid==other.grid && gridIt==other.gridIt) ||
               (!grid && grid==other.grid && connectorIt==other.connectorIt) ;
    }

    bool operator==(const CurveCirculator& otherGridIt) const{
        return grid && gridIt==otherGridIt;
    }

    bool operator==(const LinkIteratorType& otherLinkIt) const{
        return !grid && connectorIt==otherLinkIt;
    }

    SCell& operator*() const{

        if(grid){
            *element = *gridIt;
        } else{
            *element = *connectorIt;
        }
        return *element;
    }

    void operator++(){
        if(grid) ++gridIt;
        else ++connectorIt;
    }

    void operator--(){
        if(grid) --gridIt;
        else --connectorIt;
    }

private:
    SCell* element;
    CurveCirculator gridIt;
    LinkIteratorType connectorIt;
    bool grid;
};

template <typename CurveCirculator>
class GluedCurveIterator
        : public boost::iterator_facade<
                GluedCurveIterator<CurveCirculator>,
                typename CurveCirculator::value_type,
                boost::bidirectional_traversal_tag
        >
{
private:
    typedef std::vector<SCell>::const_iterator LinkIteratorType;

    CurveCirculator myIt1b,myIt1e,myIt2b,myIt2e;
    CurrentIteratorType<CurveCirculator> currentIterator;

    ConnectorType cType;


    LinkIteratorType myItLb,myItLe;
    Z2i::SCell *element;

    char iteratorStage;

    bool myFlagIsValid;

public:

    GluedCurveIterator();

    GluedCurveIterator(LinkIteratorType itLb,
                       LinkIteratorType itLe,
                         CurveCirculator it1b,
                         CurveCirculator it1e,
                         CurveCirculator it2b,
                         CurveCirculator it2e,
                         ConnectorType cType,
                         bool theEnd=false);

    GluedCurveIterator(const GluedCurveIterator& other);
    GluedCurveIterator<CurveCirculator>& operator =(const GluedCurveIterator& other);

    ~GluedCurveIterator(){delete element;};

    inline SCell linkSurfel(){ return *myItLb; }

    inline LinkIteratorType linkSurfelsBegin(){return myItLb;};
    inline LinkIteratorType linkSurfelsEnd(){return myItLe;};

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
