#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;


#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include <boost/filesystem/path.hpp>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "../utils.h"
#include "../FlowGraph/weightSettings.h"

#include "../FlowGraph/FlowGraphBuilder.h"
#include "../FlowGraph/FlowGraphQuery.h"
#include "../FlowGraph/FlowGraphDebug.h"

using namespace UtilsTypes;

namespace Development{
    bool solveShift;
    bool crossElement;

    bool makeConvexArcs;
    bool invertGluedArcs;
};

double computeEnergyValue(Image2D& image,
                          ImageFlowData& model)
{
    ImageFlowData imf(image);
    imf.init(model.getFlowMode(),model.getGluedCurveLength());

    std::map<Z2i::SCell,double> weightMap;


    FlowGraph f;
    FlowGraphBuilder fgb(f,imf,weightMap);


    ListDigraph::ArcMap<bool> internalArcs(f.graph(),false);
    FlowGraphQuery::filterArcs(f,
                               internalArcs,
                               FlowGraph::ArcType::InternalCurveArc,
                               true);


    double s = 0;
    for(ListDigraph::ArcIt a(f.graph());a!=INVALID;++a)
    {
        if(internalArcs[a]){
            s+= f.weight(a);
        }
    }

    return s;

}


void computeFlow(SegCut::Image2D& image,
                 unsigned int gluedCurveLength,
                 std::map<Z2i::SCell,double>& weightMap,
                 std::string outputFolder,
                 std::string suffix)
{
    ImageFlowData imageFlowData(image);

    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,
                       gluedCurveLength);


    setArcsWeight(imageFlowData,weightMap);

    FlowGraph fg;
    FlowGraphBuilder fgb(fg,imageFlowData,weightMap);


    FlowGraphDebug fgd(fg);
    fgd.drawCutGraph("../cmake-build-debug","cutDebug");

    ListDigraph& myGraph = fg.graph();

    {
        std::cout << "GluedEdgePairs" << std::endl;

        FlowGraphQuery::ArcPairSet gluedArcPairSet;
        FlowGraphQuery::gluedArcPairSet(fg,gluedArcPairSet);

        ListDigraph::ArcMap<int> highligthGluedEdgePair(myGraph, 0);
        int i = 1;
        for (FlowGraphQuery::ArcPairIterator api = gluedArcPairSet.begin(); api != gluedArcPairSet.end(); ++api) {
            highligthGluedEdgePair[api->first] = i;
            highligthGluedEdgePair[api->second] = i;

            std::cout << myGraph.id( api->first ) << "::" << myGraph.id(api->second) << std::endl;
            i++;
        }
        fgd.highlightArcs(highligthGluedEdgePair, outputFolder,"highlightGluedEdgePair");
    }


    {
        std::cout << "DetourArcs" << std::endl;

        FlowGraphQuery::DetourArcMap detourArcMap;
        FlowGraphQuery::detourArcMap(fg,detourArcMap);

        ListDigraph::ArcMap<int> highlightDetourArcs(myGraph, 0);
        ListDigraph::ArcMap<int> highligthGluedEdgePairSecond(myGraph, 0);
        int i=1;
        for (FlowGraphQuery::DetourArcMapIterator dami = detourArcMap.begin(); dami != detourArcMap.end(); ++dami) {
            FlowGraphQuery::ArcPair key = dami->first;

            for (FlowGraphQuery::DetourArcIterator dai = dami->second.begin(); dai != dami->second.end(); ++dai) {
                highlightDetourArcs[*dai] = i;
                std::cout << myGraph.id(*dai) << std::endl;
            }

            highligthGluedEdgePairSecond[key.first] = i;
            highligthGluedEdgePairSecond[key.second] = i;

            i++;
        }
        fgd.highlightArcs(highlightDetourArcs, outputFolder,"highlightDetourArcs");

        fgd.highlightArcs(highligthGluedEdgePairSecond, outputFolder,"highlightGluedEdgePairSecond");


    }


}

int main() {
    Development::solveShift = false;
    Development::crossElement = false;

    Development::makeConvexArcs = false;
    Development::invertGluedArcs = false;

    unsigned int gluedCurveLength = 5;

    std::string imgPath = "../images/flow-evolution/single_square.pgm";
    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(imgPath);

    std::string outputFolder = "../output/testModules/FlowGraphQuery";


    KSpace KImage;
    setKImage(imgPath, KImage);

    std::map<Z2i::SCell,double> weightMap;
    computeFlow(image,
                gluedCurveLength,
                weightMap,
                outputFolder,
                "");


    return 0;
}