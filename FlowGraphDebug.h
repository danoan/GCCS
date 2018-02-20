#ifndef SEGBYCUT_FLOWGRAPHDEBUG_H
#define SEGBYCUT_FLOWGRAPHDEBUG_H

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

#include "FlowGraphBuilder.h"

using namespace lemon;

class FlowGraphDebug{
public:
    FlowGraphDebug(FlowGraphBuilder& fgb):fgb(fgb){};

    void drawCutGraph(std::string outputFolder,
                      std::string suffix);

    void highlightArcs(ListDigraph::ArcMap<int>& arcColors,
                       std::string imageOutputPath);


private:

    void getCutFilter(Preflow<ListDigraph,ListDigraph::ArcMap<double> >& flow,
                      ListDigraph::ArcMap<bool>& arcsInTheCut);


private:
    FlowGraphBuilder& fgb;
};

#endif //SEGBYCUT_FLOWGRAPHDEBUG_H
