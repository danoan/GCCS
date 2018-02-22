#ifndef FLOWGRAPHQUERY_H
#define FLOWGRAPHQUERY_H

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;

#include "FlowGraphBuilder.h"

class SCellToArc
{
public:
    typedef DGtal::Z2i::SCell SCell;

    SCellToArc(FlowGraphBuilder::SCellArcMap& map):scellArc(map){};

    ListDigraph::Arc operator()(const SCell& a) const
    {
        return scellArc[a];
    }

private:
    FlowGraphBuilder::SCellArcMap& scellArc;
};

class FlowGraphQuery
{
public:
    typedef std::pair< ListDigraph::Arc, ListDigraph::Arc > ArcPair;
    typedef std::vector<ArcPair>::const_iterator ArcPairIterator;

    typedef std::map< ArcPair, std::set< ListDigraph::Arc > >::const_iterator DetourArcMapIterator;
    typedef std::set< ListDigraph::Arc >::const_iterator DetourArcIterator;
public:
    FlowGraphQuery(FlowGraphBuilder& fgb):fgb(fgb){};

    void sourceComponent(ListDigraph::ArcMap<bool>& cutFilter);

    void sourceComponent(ListDigraph::NodeMap<bool>& nodeFilter,
                         ListDigraph::ArcMap<bool>& arcFilter,
                         ListDigraph::ArcMap<bool>& cutFilter);

    template<typename MappedType>
    void filterArcs(ListDigraph::ArcMap<MappedType>& retArcFilter,
                    FlowGraphBuilder::ArcType at,
                    MappedType valueForTrue);

    template<typename MappedType>
    static void filterArcs(FlowGraphBuilder& fgb,
                           ListDigraph::ArcMap<MappedType>& retArcFilter,
                           FlowGraphBuilder::ArcType at,
                           MappedType valueForTrue);

    ArcPairIterator gluedEdgePairsBegin();
    ArcPairIterator gluedEdgePairsEnd();

    DetourArcMapIterator detourArcsBegin();
    DetourArcMapIterator detourArcsEnd();

private:
    void setGluedEdgePairsOnCut();
    void setDetourArcs();

    template<typename IteratorType>
    void insertSCellFromArc(std::set<ListDigraph::Arc>& v,
                            SCellToArc& staFunctor,
                            IteratorType begin,
                            int items);

    template<typename IteratorType,typename LimitIteratorType>
    void insertSCellFromArc(std::set<ListDigraph::Arc>& v,
                            SCellToArc& staFunctor,
                            IteratorType begin,
                            LimitIteratorType limitIterator,
                            int items);

    //Deprecated
    void setDetourArcsAsExternalArcs();

private:
    FlowGraphBuilder& fgb;

    std::vector< ArcPair > gluedEdgePairsOnCut;
    std::map< ArcPair, std::set< ListDigraph::Arc > > detourArcs;
};


template<typename MappedType>
void FlowGraphQuery::filterArcs(ListDigraph::ArcMap<MappedType>& retArcFilter,
                                FlowGraphBuilder::ArcType at,
                                MappedType valueForTrue) {
    FlowGraphQuery::filterArcs(fgb,
                               retArcFilter,
                               at,
                               valueForTrue);
}

template<typename MappedType>
void FlowGraphQuery::filterArcs(FlowGraphBuilder& fgb,
                                ListDigraph::ArcMap<MappedType>& retArcFilter,
                                FlowGraphBuilder::ArcType at,
                                MappedType valueForTrue)
{
    for(ListDigraph::ArcIt a(fgb.graph());a!=INVALID;++a)
    {
        if(fgb.getArcType(a)==at)
        {
            retArcFilter[a] =valueForTrue;
        }

    }
}

template<typename IteratorType>
void FlowGraphQuery::insertSCellFromArc(std::set<ListDigraph::Arc>& v,
                                        SCellToArc& staFunctor,
                                        IteratorType begin,
                                        int items)
{
    ConstIteratorAdapter<IteratorType,SCellToArc,ListDigraph::Arc> itAdapter(begin,
                                                                             staFunctor);

    unsigned long int initialSize = v.size();
    v.insert(itAdapter,walkCirculator(itAdapter,items) );
    if(initialSize+items!=v.size()){
        std::cout << "Collision of " << initialSize+items-v.size() << " items" << std::endl;
    }
}

template<typename IteratorType,typename LimitIteratorType>
void FlowGraphQuery::insertSCellFromArc(std::set<ListDigraph::Arc>& v,
                                        SCellToArc& staFunctor,
                                        IteratorType begin,
                                        LimitIteratorType limitIterator,
                                        int items)
{
    ConstIteratorAdapter<IteratorType,SCellToArc,ListDigraph::Arc> itAdapter(begin,
                                                                             staFunctor);

    int n =0;
    ListDigraph::Arc a = staFunctor(*limitIterator);
    while(n<items)
    {
        if(*itAdapter==a) break;
//        std::cout << fgb.graph().id(*itAdapter) << "?" << fgb.arcType[*itAdapter] << "::" << fgb.graph().id(a) << "?" << fgb.arcType[a] << std::endl;
        v.insert(*itAdapter);
        ++itAdapter;
        ++n;
    }

    if(n!=items){
        std::cout << "Limit condition reached: " << n << " of " << items << " items included" << std::endl;
    }
}


#endif