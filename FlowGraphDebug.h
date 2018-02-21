#ifndef SEGBYCUT_FLOWGRAPHDEBUG_H
#define SEGBYCUT_FLOWGRAPHDEBUG_H

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

#include "FlowGraphBuilder.h"
#include "Flow.h"

using namespace lemon;

class FlowGraphDebug{
public:
    FlowGraphDebug(Flow& flow):fgb(flow.fgb),
                               fgq(fgb){};

    FlowGraphDebug(FlowGraphBuilder& fgb):fgb(fgb),
                                          fgq(fgb){};

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

private:
    FlowGraphBuilder& fgb;
    FlowGraphQuery fgq;
};

template<typename MapType,typename IteratorType>
class FilterComposer
{
public:
    typedef MapType FilterType;

public:
    FilterComposer(ListDigraph& graph):baseGraph(graph),
                                       initialMap(graph,false){};

    FilterComposer& operator*(FilterType& otherMap)
    {
        for(IteratorType it(baseGraph);it!=INVALID;++it)
        {
            if(otherMap[it] && initialMap[it])
            {
                initialMap[it]=true;
            }
        }

        return *this;
    }

    FilterComposer& operator+(FilterType& otherMap)
    {
        for(IteratorType it(baseGraph);it!=INVALID;++it)
        {
            if(otherMap[it])
            {
                initialMap[it]=true;
            }
        }

        return *this;
    }

    FilterType& operator()(){return initialMap;}

public:
    ListDigraph& baseGraph;
    FilterType initialMap;
};


#endif //SEGBYCUT_FLOWGRAPHDEBUG_H
