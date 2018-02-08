
#include <boost/filesystem/operations.hpp>
#include "FlowGraphDebug.h"

void FlowGraphDebug::drawCutGraph(std::string outputFolder,
                                  std::string suffix)
{
    FlowGraphBuilder::Flow flow = fgb.preparePreFlow();
    flow.run();

    ListDigraph::ArcMap<bool> arcsInTheCut(fgb.graph(),false);
    getCutFilter(flow,
                 arcsInTheCut);

    ListDigraph::NodeMap<bool> nodeFilter(fgb.graph(),true);
    ListDigraph::ArcMap<bool> arcFilter(fgb.graph(),true);

    nodeFilter[fgb.source()] = false;
    nodeFilter[fgb.target()] = false;

    FlowGraphBuilder::SubGraph subgraph(fgb.graph(),
                                        nodeFilter,
                                        arcFilter);

    Palette palette;
    FlowGraphBuilder::SubGraph::ArcMap<int> arcColors(subgraph,0);

    for(FlowGraphBuilder::SubGraph::ArcIt a(subgraph);a!=INVALID;++a)
    {
        if(arcsInTheCut[a]){
            arcColors[a] = 1;
        }
    }


    std::string imageOutputPath = outputFolder + "/cutGraph" + suffix + ".eps";
    boost::filesystem::path p1(imageOutputPath.c_str());
    p1.remove_filename();
    boost::filesystem::create_directories(p1);

    graphToEps(subgraph,imageOutputPath.c_str())
            .coords(fgb.coordsMap())
            .nodeScale(.0005)
            .arcWidthScale(0.0005)
            .arcColors((composeMap(palette,arcColors)))
            .drawArrows()
            .arrowWidth(0.25)
            .run();
}



void FlowGraphDebug::getCutFilter(Preflow<ListDigraph,ListDigraph::ArcMap<double> >& flow,
                                  ListDigraph::ArcMap<bool>& arcsInTheCut)
{
    ListDigraph::NodeMap<bool> nodesInTheCut(fgb.graph());

    for(ListDigraph::NodeIt n(fgb.graph());n!=INVALID;++n){
        if( flow.minCut(n) ){
            nodesInTheCut[n] = true;
        }else{
            nodesInTheCut[n] = false;
        }
    }
    nodesInTheCut[fgb.source()]=false;

    for(ListDigraph::ArcIt a(fgb.graph());a!=INVALID;++a){
        ListDigraph::Node s =fgb.graph().source(a);
        ListDigraph::Node t =fgb.graph().target(a);

        if( (nodesInTheCut[s] && !nodesInTheCut[t])  )
        {
            arcsInTheCut[a] = true;
        }
    }

}



