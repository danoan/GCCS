
#include <boost/filesystem/operations.hpp>
#include "FlowGraphDebug.h"

void FlowGraphDebug::colorizeArcs(ListDigraph::ArcMap<int>& arcColors,
                                  ListDigraph::ArcMap<bool>& arcFilter,
                                  int color)
{
    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(arcFilter[a])
        {
            arcColors[a] = color;
        }
    }
}

void FlowGraphDebug::drawFlowGraph(std::string outputFolder,
                                    std::string suffix)
{
    ListDigraph::ArcMap<bool> internalArcs(fg.graph(),false);
    FlowGraphQuery::filterArcs(fg,internalArcs,FlowGraph::InternalCurveArc,true);

    ListDigraph::ArcMap<bool> externalArcs(fg.graph(),false);
    FlowGraphQuery::filterArcs(fg,externalArcs,FlowGraph::ExternalCurveArc,true);

    ListDigraph::ArcMap<bool> escapeArcs(fg.graph(),false);
    FlowGraphQuery::filterArcs(fg,escapeArcs,FlowGraph::EscapeArc,true);

    ListDigraph::NodeMap<bool> noTerminal(fg.graph(),true);
    noTerminal[fg.source()]=false;
    noTerminal[fg.target()]=false;

    ListDigraph::ArcMap<bool> allArcs(fg.graph(),true);

//    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> FC(fg.graph(),internalArcs);
//    FC+externalArcs+escapeArcs;

    ListDigraph::ArcMap<int> arcColors(fg.graph(),0);
    colorizeArcs(arcColors,allArcs,1);
    colorizeArcs(arcColors,internalArcs,2);
    colorizeArcs(arcColors,externalArcs,3);
    colorizeArcs(arcColors,escapeArcs,4);

    suffix = "FlowGraph-" + suffix;
    highlightArcs(noTerminal,allArcs,arcColors,outputFolder,suffix);
}

void FlowGraphDebug::drawCutGraph(std::string outputFolder,
                                  std::string suffix)
{
    FlowGraph::FlowComputer flowComputer = fg.prepareFlow();
    flowComputer.run();

    ListDigraph::ArcMap<bool> arcsInTheCut(fg.graph(),false);
    FlowGraphQuery::sourceComponent(fg,arcsInTheCut);

    suffix = "CutGraph-" + suffix;
    highlightArcs(arcsInTheCut,outputFolder,suffix);
}

void FlowGraphDebug::drawRefundGraph(std::string outputFolder,
                                     std::string suffix)
{
    ListDigraph::ArcMap<bool> noArcs(fg.graph(),false);
    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> filterComposer(fg.graph(),noArcs);

    ListDigraph::ArcMap<bool> refundArcsFilter(fg.graph(),false);
    FlowGraphQuery::filterArcs(fg,
                               refundArcsFilter,
                               FlowGraph::ArcType::RefundArc,true);

    ListDigraph::ArcMap<bool> externalArcsFilter(fg.graph(),
                                                 false);

    FlowGraphQuery::filterArcs(fg,
                               externalArcsFilter,
                               FlowGraph::ArcType::ExternalCurveArc,
                               true);

    filterComposer+refundArcsFilter+externalArcsFilter;

    suffix = suffix + "-RefundArcs";
    highlightArcs(filterComposer.initialMap,outputFolder,suffix);
}


void FlowGraphDebug::highlightArcs(ListDigraph::ArcMap<bool>& arcFilter,
                                   std::string imageOutputFolder,
                                   std::string suffix)
{

    ListDigraph::NodeMap<bool> allNodes(fg.graph(),true);
    ListDigraph::ArcMap<bool> allArcs(fg.graph(),true);

    ListDigraph::ArcMap<int> arcColors(fg.graph(),0);

    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(arcFilter[a])
        {
            arcColors[a] = 1;
        }
    }

    highlightArcs(allNodes,
                  arcFilter,
                  arcColors,
                  imageOutputFolder,
                  suffix);
}

void FlowGraphDebug::highlightArcs(ListDigraph::ArcMap<int>& arcColors,
                                   std::string imageOutputFolder,
                                   std::string suffix)
{

    ListDigraph::NodeMap<bool> allNodes(fg.graph(),true);
    ListDigraph::ArcMap<bool> allArcs(fg.graph(),true);

    highlightArcs(allNodes,
                  allArcs,
                  arcColors,
                  imageOutputFolder,
                  suffix);
}


void FlowGraphDebug::highlightArcs(ListDigraph::NodeMap<bool>& nodeFilter,
                                   ListDigraph::ArcMap<bool>& arcFilter,
                                   ListDigraph::ArcMap<int>& arcColors,
                                   std::string imageOutputFolder,
                                   std::string suffix)
{

    Palette palette;

    std::string imageOutputPath = imageOutputFolder + "/" + suffix + ".eps";

    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    FlowGraphQuery::SubGraph sg(fg.graph(),
                                nodeFilter,
                                arcFilter);

    graphToEps(sg,imageOutputPath.c_str())
            .coords(fg.coordsMap())
            .nodeScale(.0005)
            .arcWidthScale(0.0005)
            .arcColors((composeMap(palette,arcColors)))
            .drawArrows()
            .arrowWidth(0.25)
            .run();
}

double FlowGraphDebug::energyValue(ListDigraph::ArcMap<bool>& arcFilter)
{
    double s = 0;
    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(arcFilter[a])
        {
            s+=fg.weight(a);
        }

    }
    return s;
}




