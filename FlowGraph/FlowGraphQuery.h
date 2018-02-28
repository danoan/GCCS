#ifndef FLOWGRAPHQUERY_H
#define FLOWGRAPHQUERY_H

#include <set>
#include <queue>
#include <stack>

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

#include "DGtal/base/ConstIteratorAdapter.h"

#include "FlowGraph.h"


class SCellToArc
{
public:
    typedef FlowGraph::SCell SCell;

public:
    SCellToArc(FlowGraph& fg):fg(fg){};

    lemon::ListDigraph::Arc operator()(const SCell& scell) const
    {
        return fg.arc(scell);
    }

private:
    FlowGraph& fg;
};




template<typename MapType,typename IteratorType>
class FilterComposer
{
public:
    typedef MapType FilterType;
    typedef lemon::ListDigraph ListDigraph;
public:
    FilterComposer(ListDigraph& graph,
                   FilterType& initialMap):baseGraph(graph),
                                           initialMap(initialMap){};

    FilterComposer& operator*(FilterType& otherMap)
    {
        for(IteratorType it(baseGraph);it!=lemon::INVALID;++it)
        {
            if(otherMap[it] && initialMap[it])
            {
                initialMap[it]=true;
            }else
            {
                initialMap[it]=false;
            }

        }

        return *this;
    }

    FilterComposer& operator+(FilterType& otherMap)
    {
        for(IteratorType it(baseGraph);it!=lemon::INVALID;++it)
        {
            if(otherMap[it])
            {
                initialMap[it]=true;
            }
        }

        return *this;
    }

    FilterComposer& operator-(FilterType& otherMap)
    {
        for(IteratorType it(baseGraph);it!=lemon::INVALID;++it)
        {
            if(otherMap[it])
            {
                initialMap[it]=false;
            }
        }

        return *this;
    }

    FilterType& operator()(){return initialMap;}

public:
    ListDigraph& baseGraph;
    FilterType&  initialMap;
};





class FlowGraphQuery
{
public:
    typedef lemon::ListDigraph ListDigraph;

    typedef std::pair< ListDigraph::Arc, ListDigraph::Arc > ArcPair;
    typedef std::set<ArcPair> ArcPairSet;
    typedef std::set<ArcPair>::const_iterator ArcPairIterator;

    typedef std::map< ArcPair, std::set< ListDigraph::Arc > > DetourArcMap;
    typedef DetourArcMap::const_iterator DetourArcMapIterator;
    typedef std::set< ListDigraph::Arc >::const_iterator DetourArcIterator;


    typedef lemon::SubDigraph< ListDigraph,ListDigraph::NodeMap<bool> > SubGraph;

    typedef FlowGraph::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef std::pair< Point, Point > ArcKey;

public:
    
    static void sourceComponent(FlowGraph& fg,
                                ListDigraph::ArcMap<bool>& cutFilter);

    
    static void sourceComponent(FlowGraph& fg,
                                ListDigraph::NodeMap<bool>& nodeFilter,
                                ListDigraph::ArcMap<bool>& arcFilter,
                                ListDigraph::ArcMap<bool>& cutFilter);

    static void selectArcsInCut(FlowGraph& fg,
                                ListDigraph::ArcMap<bool>& selectedArcs,
                                FlowGraph::ArcType at);


    static void coreArcs(FlowGraph& fg,
                         ListDigraph::ArcMap<bool>& coreArcs);

    static void pixelsFilter(FlowGraph& fg,
                             ListDigraph::NodeMap<bool>& pixelsFilter);


    static void detourArcFilter(FlowGraph& fg,
                                ListDigraph::ArcMap<bool>& arcFilter,
                                DetourArcMap& detourArcMap);


    template<typename MappedType>
    static void filterArcs(FlowGraph& fg,
                           ListDigraph::ArcMap<MappedType>& retArcFilter,
                           FlowGraph::ArcType at,
                           MappedType valueForTrue);

    static int arcDistance(FlowGraph& fg,
                           ListDigraph::Arc& a1,
                           ListDigraph::Arc& a2);

    static void gluedArcPairSet(FlowGraph& fg,
                                std::set<ArcPair>& gluedArcPairSet);

    static void detourArcMap(FlowGraph& fg,
                             DetourArcMap& detourArcMap);

    static double computeEnergyValue(FlowGraph& fg,
                                     ListDigraph::ArcMap<bool>& arcFilter);

    static FlowGraphQuery::ArcKey arcKey(FlowGraph& fg,
                                         ListDigraph::Arc& arc);

    static void arcFilterConversion(ListDigraph::ArcMap<bool>& arcMap,
                                    FlowGraph& fg,
                                    ListDigraph::ArcMap<bool>& newArcMap,
                                    FlowGraph& newFg);


    static double cutValue(FlowGraph& fg);


private:
    template<typename IteratorType>
    static void insertSCellFromArc(FlowGraph& fg,
                                   std::set<ListDigraph::Arc>& v,
                                   SCellToArc& staFunctor,
                                   IteratorType begin,
                                   int items);

    template<typename IteratorType,typename LimitIteratorType>
    static void insertSCellFromArc(FlowGraph& fg,
                                   std::set<ListDigraph::Arc>& v,
                                   SCellToArc& staFunctor,
                                   IteratorType begin,
                                   LimitIteratorType limitIterator,
                                   int items);

    template<typename MyCirculator>
    static MyCirculator moveIterator(const MyCirculator c, int w);


};


template<typename MappedType>
void FlowGraphQuery::filterArcs(FlowGraph& fg,
                                ListDigraph::ArcMap<MappedType>& retArcFilter,
                                FlowGraph::ArcType at,
                                MappedType valueForTrue)
{
    for(ListDigraph::ArcIt a(fg.graph());a!=lemon::INVALID;++a)
    {
        if(fg.arcType(a)==at)
        {
            retArcFilter[a] =valueForTrue;
        }

    }
}

template<typename IteratorType>
void FlowGraphQuery::insertSCellFromArc(FlowGraph& fg,
                                        std::set<ListDigraph::Arc>& v,
                                        SCellToArc& staFunctor,
                                        IteratorType begin,
                                        int items)
{
    DGtal::ConstIteratorAdapter<IteratorType,SCellToArc,ListDigraph::Arc> itAdapter(begin,
                                                                                    staFunctor);

    unsigned long int initialSize = v.size();
    v.insert(itAdapter,FlowGraphQuery::moveIterator(itAdapter,items) );
    if(initialSize+items!=v.size())
    {
//        std::cout << "Collision of " << initialSize+items-v.size() << " items" << std::endl;
    }
}

template<typename IteratorType,typename LimitIteratorType>
void FlowGraphQuery::insertSCellFromArc(FlowGraph& fg,
                                        std::set<ListDigraph::Arc>& v,
                                        SCellToArc& staFunctor,
                                        IteratorType begin,
                                        LimitIteratorType limitIterator,
                                        int items)
{
    DGtal::ConstIteratorAdapter<IteratorType,SCellToArc,ListDigraph::Arc> itAdapter(begin,
                                                                                    staFunctor);

    unsigned long int initialSize = v.size();
    ListDigraph::Arc a = staFunctor(*limitIterator);
    int n=0;
    while(n<items)
    {
        if(*itAdapter==a) break;

        std::cout << fg.pixel( fg.source(*itAdapter) ) << " -> " 
                  <<  fg.pixel( fg.target(*itAdapter) ) << "::ID " 
                  << fg.id(*itAdapter) << std::endl;

        v.insert(*itAdapter);
        ++itAdapter;
        n++;
    }

    if(initialSize+items!=v.size()){
        std::cout << "(*)Collision of " << initialSize+items-v.size() << " items" << std::endl;
    }
}

template<typename MyCirculator>
MyCirculator FlowGraphQuery::moveIterator(const MyCirculator c, int w){
    MyCirculator cc = c;
    while(w!=0){
        w>0?cc++:cc--;
        w>0?w--:w++;
    }

    return cc;
}


#endif