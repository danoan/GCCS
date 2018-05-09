#ifndef SEGBYCUT_FLOWGRAPHDEBUG_H
#define SEGBYCUT_FLOWGRAPHDEBUG_H

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

#include "FlowGraphBuilder.h"
#include "../Flows/RefundFlow/RefundFlow.h"

using namespace lemon;

class FlowGraphDebug{
public:
    FlowGraphDebug(FlowGraph& fg):fg(fg){};

    void drawFlowGraph(std::string outputFolder,
                       std::string suffix);

    void colorizeArcs(ListDigraph::ArcMap<int>& arcColors,
                      ListDigraph::ArcMap<bool>& arcFilter,
                      int color);

    void drawCutGraph(std::string outputFolder,
                      std::string suffix);

    void drawRefundGraph(std::string outputFolder,
                         std::string suffix);
    
    
    void highlightArcs(ListDigraph::ArcMap<bool>& arcFilter,
                       std::string imageOutputFolder,
                       std::string suffix);

    void highlightArcs(ListDigraph::ArcMap<int>& arcColors,
                       std::string imageOutputFolder,
                       std::string suffix);

    void highlightArcs(ListDigraph::NodeMap<bool>& nodeFilter,
                       ListDigraph::ArcMap<bool>& arcFilter,
                       ListDigraph::ArcMap<int>& arcColors,
                       std::string imageOutputFolder,
                       std::string suffix);

    double energyValue(ListDigraph::ArcMap<bool>& arcFilter);

private:
    FlowGraph& fg;
};


#endif //SEGBYCUT_FLOWGRAPHDEBUG_H
