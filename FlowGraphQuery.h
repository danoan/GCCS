#ifndef FLOWGRAPHQUERY_H
#define FLOWGRAPHQUERY_H

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;

#include "FlowGraphBuilder.h"

class FlowGraphQuery
{
public:
    typedef std::pair< ListDigraph::Arc, ListDigraph::Arc > ArcPair;
    typedef std::vector<ArcPair>::const_iterator ArcPairIterator;

    typedef std::map< ArcPair, std::vector< ListDigraph::Arc > >::const_iterator DetourArcMapIterator;
    typedef std::vector< ListDigraph::Arc >::const_iterator DetourArcIterator;
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

private:
    FlowGraphBuilder& fgb;

    std::vector< ArcPair > gluedEdgePairsOnCut;
    std::map< ArcPair, std::vector< ListDigraph::Arc > > detourArcs;
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


#endif