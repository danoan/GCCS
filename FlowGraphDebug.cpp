
#include <boost/filesystem/operations.hpp>
#include "FlowGraphDebug.h"

void FlowGraphDebug::drawCutGraph(std::string outputFolder,
                                  std::string suffix)
{
    FlowGraphBuilder::FlowComputer flowComputer = fgb.preparePreFlow();
    flowComputer.run();

    ListDigraph::ArcMap<bool> arcsInTheCut(fgb.graph(),false);
    fgq.sourceComponent(arcsInTheCut);

    suffix = "CutGraph-" + suffix;
    highlightArcs(arcsInTheCut,outputFolder,suffix);
}

void FlowGraphDebug::drawRefundGraph(std::string outputFolder,
                                     std::string suffix)
{
    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> filterComposer(fgb.graph());

    ListDigraph::ArcMap<bool> refundArcsFilter(fgb.graph(),false);
    fgq.filterArcs(refundArcsFilter,FlowGraphBuilder::ArcType::RefundArc,true);

    ListDigraph::ArcMap<bool> externalArcsFilter(fgb.graph(),false);
    fgq.filterArcs(externalArcsFilter,FlowGraphBuilder::ArcType::ExternalCurveArc,true);

    filterComposer+refundArcsFilter+externalArcsFilter;

    suffix = suffix + "-RefundArcs";
    highlightArcs(filterComposer(),outputFolder,suffix);
}


void FlowGraphDebug::highlightArcs(ListDigraph::ArcMap<bool>& arcFilter,
                                   std::string imageOutputFolder,
                                   std::string suffix)
{

    ListDigraph::NodeMap<bool> allNodes(fgb.graph(),true);
    ListDigraph::ArcMap<bool> allArcs(fgb.graph(),true);

    ListDigraph::ArcMap<int> arcColors(fgb.graph(),0);

    for(ListDigraph::ArcIt a(fgb.graph());a!=INVALID;++a)
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

    ListDigraph::NodeMap<bool> allNodes(fgb.graph(),true);
    ListDigraph::ArcMap<bool> allArcs(fgb.graph(),true);

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

    FlowGraphBuilder::SubGraph sg(fgb.graph(),
                                  nodeFilter,
                                  arcFilter);

    graphToEps(sg,imageOutputPath.c_str())
            .coords(fgb.coordsMap())
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
    for(ListDigraph::ArcIt a(fgb.graph());a!=INVALID;++a)
    {
        if(arcFilter[a])
        {
            s+=fgb.arcWeight[a];
        }

    }
    return s;
}




