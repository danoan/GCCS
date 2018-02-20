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

    ArcPairIterator gluedEdgePairsBegin();
    ArcPairIterator gluedEdgePairsEnd();

    DetourArcMapIterator detourArcsBegin();
    DetourArcMapIterator detourArcsEnd();

private:
    void setGluedEdgePairsOnCut();
    void setDetourArcs();

    FlowGraphBuilder::SubGraph sourceComponnent(ListDigraph::NodeMap<bool>& nodeFilter,
                                                ListDigraph::ArcMap<bool>& arcFilter,
                                                ListDigraph::ArcMap<bool>& cutFilter);

private:
    FlowGraphBuilder& fgb;

    std::vector< ArcPair > gluedEdgePairsOnCut;
    std::map< ArcPair, std::vector< ListDigraph::Arc > > detourArcs;
};

#endif