#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "utils/utils.h"
#include "FlowGraph/FlowGraphUtils.h"

#include "FlowGraph/FlowGraphBuilder.h"
#include "FlowGraph/ImageFlowDataDebug.h"
#include "FlowGraph/FlowGraphDebug.h"
#include "utils/Artist.h"

using namespace UtilsTypes;

namespace Development{
    bool solveShift=false;
    bool crossElement=false;

    bool lambdaEstimator=false;
    bool pessimistEstimator=false;

    bool makeConvexArcs=false;
    bool invertGluedArcs=false;
};

int main()
{
    Image2D imageOriginal = GenericReader<Image2D>::import("../images/comedic-mi-parcours/line/line-artifacts.pgm");
    Image2D imageIteration = GenericReader<Image2D>::import("../images/comedic-mi-parcours/line/8.pgm");

    int gluedCurveLenght = 5;

    ImageFlowData imageFlowData(imageOriginal);
    imageFlowData.init(ImageFlowData::DilationOnly,gluedCurveLenght);

    std::string outputFolder = "../output/comedic-mi-parcours";


    Board2D board;
    board << imageFlowData.getMostInnerCurve();
    board << imageFlowData.getMostOuterCurve();

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        for(auto itR=it->gcsRangeBegin();itR!=it->gcsRangeEnd();++itR)
        {
            board << *itR->first.connectorsBegin();
        }
    }
    board.save( (outputFolder+"/base-grid.eps").c_str() );





    std::map<DGtal::Z2i::SCell,double> weightMap;
    setArcsWeight(imageFlowData,weightMap);
    FlowGraph fg;
    FlowGraphBuilder fgb(fg,imageFlowData,weightMap);

    FlowGraphDebug fgd(fg);
    ListDigraph::ArcMap<bool> sourceArcs(fg.graph(),false);
    ListDigraph::ArcMap<bool> targetArcs(fg.graph(),false);
    ListDigraph::ArcMap<bool> notTerminalArcs(fg.graph(),true);
    FlowGraphQuery::filterArcs(fg,sourceArcs,FlowGraph::ArcType::SourceArc,true);
    FlowGraphQuery::filterArcs(fg,sourceArcs,FlowGraph::ArcType::TargetArc,true);

    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> fc(fg.graph(),notTerminalArcs);
    fc-sourceArcs-targetArcs;

    fgd.highlightArcs(fc(),outputFolder,"flow-graph");


    ListDigraph::ArcMap<bool> cutFilter(fg.graph(),false);
    FlowGraphQuery::sourceComponent(fg,cutFilter);
    ListDigraph::ArcMap<int> arcColors(fg.graph(),0);
    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(cutFilter[a])
        {
            arcColors[a]=1;
        }
    }

    ListDigraph::NodeMap<bool> allNodes(fg.graph(),true);
    fgd.highlightArcs(allNodes,fc(),arcColors,outputFolder,"highlight-cut-Arcs");


    ImageFlowData imf2(imageIteration);
    imf2.init(ImageFlowData::DilationOnly,gluedCurveLenght);

    board.clear();
    board << imf2.getMostInnerCurve();
    board.save("../output/comedic-mi-parcours/iteration.eps");


    {
        ImageFlowData imf(imageOriginal);
        imf.init(ImageFlowData::DilationOnly,5);

        Artist EA(imf.getKSpace(), board);
        /*Maximal Stabbing Circles*/
        EA.board.clear(DGtal::Color::White);
        EA.drawMaximalStabbingCircles(imf.getMostInnerCurve());
        std::string outputFilePath = outputFolder + "/Original-MaximalStabbingCircles.eps";
        EA.board.save(outputFilePath.c_str());
    }

    {
        ImageFlowData imf(imageIteration);
        imf.init(ImageFlowData::DilationOnly,5);

        Artist EA(imf.getKSpace(), board);
        /*Maximal Stabbing Circles*/
        EA.board.clear(DGtal::Color::White);
        EA.drawMaximalStabbingCircles(imf.getMostInnerCurve());
        std::string outputFilePath = outputFolder + "/Iteration-MaximalStabbingCircles.eps";
        EA.board.save(outputFilePath.c_str());
    }

    return 0;
}