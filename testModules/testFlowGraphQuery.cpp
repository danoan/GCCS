
#include "../weightSettings.h"

#include "../FlowGraphBuilder.h"
#include "../FlowGraphQuery.h"
#include "../FlowGraphDebug.h"

namespace Development{
    bool solveShift;
    bool crossElement;

    bool makeConvexArcs;
    bool invertGluedArcs;
};

void computeFlow(SegCut::Image2D& image,
                 unsigned int gluedCurveLength,
                 std::map<Z2i::SCell,double>& weightMap,
                 std::string outputFolder,
                 std::string suffix)
{
    ImageFlowData imageFlowData(image);

    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,
                       gluedCurveLength);


    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        setGridCurveWeight(it->curve,
                           imageFlowData.getKSpace(),
                           weightMap
        );
    }

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        setGluedCurveWeight( it->gcsRangeBegin(),
                             it->gcsRangeEnd(),
                             imageFlowData.getKSpace(),
                             gluedCurveLength,
                             weightMap);
    }

    FlowGraphBuilder fgb(imageFlowData);
    fgb(weightMap);

    FlowGraphDebug fgd(fgb);
    fgd.drawCutGraph("../cmake-build-debug","cutDebug");

    FlowGraphQuery fgq(fgb);
    ListDigraph& myGraph = fgb.graph();

    {
        std::cout << "GluedEdgePairs" << std::endl;
        ListDigraph::ArcMap<int> highligthGluedEdgePair(myGraph, 0);
        int i = 1;
        for (FlowGraphQuery::ArcPairIterator api = fgq.gluedEdgePairsBegin(); api != fgq.gluedEdgePairsEnd(); ++api) {
            highligthGluedEdgePair[api->first] = i;
            highligthGluedEdgePair[api->second] = i;

            std::cout << myGraph.id( api->first ) << "::" << myGraph.id(api->second) << std::endl;
            i++;
        }
        fgd.highlightArcs(highligthGluedEdgePair, "highlightGluedEdgePair.eps");
    }


    {
        std::cout << "DetourArcs" << std::endl;
        ListDigraph::ArcMap<int> highlightDetourArcs(myGraph, 0);
        int i=1;
        for (FlowGraphQuery::DetourArcMapIterator dami = fgq.detourArcsBegin(); dami != fgq.detourArcsEnd(); ++dami) {
            for (FlowGraphQuery::DetourArcIterator dai = dami->second.begin(); dai != dami->second.end(); ++dai) {
                highlightDetourArcs[*dai] = i;
                std::cout << myGraph.id(*dai) << std::endl;
            }
            i++;
        }
        fgd.highlightArcs(highlightDetourArcs, "highlightDetourArcs.eps");
    }


}

int main() {
    Development::solveShift = false;
    Development::crossElement = false;

    Development::makeConvexArcs = false;
    Development::invertGluedArcs = false;

    unsigned int gluedCurveLength = 7;

    std::string imgPath = "../images/flow-evolution/single_square.pgm";
    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(imgPath);

    std::string outputFolder = "../output/testFlowGraphQuery/single-square";


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