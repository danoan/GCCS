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

#include "utils.h"
#include "weightSettings.h"

#include "FlowGraphBuilder.h"
#include "ImageFlowData.h"
#include "ImageFlowDataDebug.h"
#include "FlowGraphDebug.h"
#include "FlowGraphQuery.h"

using namespace UtilsTypes;


void updateImage(FlowGraphBuilder& fgb,
                 ImageFlowData& imageFlowData,
                 Image2D& out)
{
    FlowGraphBuilder::Flow flow = fgb.preparePreFlow();
    flow.run();

    ListDigraph::NodeMap<bool> node_filter(fgb.graph(),false);
    ListDigraph::ArcMap<bool> arc_filter(fgb.graph(),true);

    for(ListDigraph::NodeIt n(fgb.graph());n!=INVALID;++n){
        if( flow.minCut(n) ){
            node_filter[n] = true;
        }else{
            node_filter[n] = false;
        }
    }

    node_filter[fgb.source()]=false;



    FlowGraphBuilder::SubGraph subgraph(fgb.graph(),
                                        node_filter,
                                        arc_filter);


    KSpace& KImage = imageFlowData.getKSpace();
    std::vector<Z2i::Point> coordPixelsSourceSide;
    std::vector<Z2i::SCell> pixelsInTheGraph;

    for(FlowGraphBuilder::SubGraph::NodeIt n(subgraph);n!=INVALID;++n)
    {
        Z2i::SCell pixel = fgb.pixelsMap()[n];
        Z2i::Point p = KImage.sCoords(pixel);

        coordPixelsSourceSide.push_back(p);
    }

    for(ListDigraph::NodeIt n(fgb.graph());n!=INVALID;++n)
    {
        if(fgb.source()==n || fgb.target()==n) continue;
        Z2i::SCell pixel = fgb.pixelsMap()[n];

        pixelsInTheGraph.push_back(pixel);
    }


    for(auto it=coordPixelsSourceSide.begin();it!=coordPixelsSourceSide.end();++it)
    {
        out.setValue(*it,255);
    }


    /*Fill the holes*/
    for(auto it=pixelsInTheGraph.begin();it!=pixelsInTheGraph.end();++it)
    {
        Z2i::SCell pixel = *it;
        Z2i::Point p = KImage.sCoords(pixel);
        Z2i::SCells N = KImage.sNeighborhood(pixel);

        unsigned char v = out(p);
        unsigned char nv;
        bool hole=true;
        for(int j=1;j<5;++j){
            nv = out(KImage.sCoords(N[j]));
            if(v==nv){
                hole=false;
                break;
            }
        }
        if(hole){
            out.setValue(p,nv);
        }
    }
}

void singleFlow(FlowGraphBuilder& fgb,
                ImageFlowData& imageFlowData,
                int gluedCurveLength,
                std::map<Z2i::SCell,double>& weightMap)
{

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
}

void updateAndCompare(Image2D& image,
                      int gluedCurveLength,
                      FlowGraphBuilder& fgb,
                      ImageFlowData& imageFlowData,
                      ListDigraph::ArcMap<bool>& arcFilter,
                      ListDigraph::ArcMap<double>& arcDiff)
{
    Image2D out = image;

    updateImage(fgb,
                imageFlowData,
                out);

    typedef std::pair< Z2i::RealVector,
                       Z2i::RealVector > MyKey;

    {
        std::set< MyKey > betweenGraphSet;

        for(ListDigraph::ArcIt a(fgb.graph());a!=INVALID;++a)
        {
            if(arcFilter[a]){
                Z2i::SCell source = fgb.pixelsMap()[ fgb.graph().source(a) ];
                Z2i::SCell target = fgb.pixelsMap()[ fgb.graph().target(a) ];

                betweenGraphSet.insert( MyKey( source.preCell().coordinates,
                                               source.preCell().coordinates );
            }
        }
        

        std::map<Z2i::SCell, double> weightMap;
        ImageFlowData tempImageFlowData(image);

        imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,
                           gluedCurveLength);

        FlowGraphBuilder tempFGB(imageFlowData);
        singleFlow(tempFGB,
                   imageFlowData,
                   gluedCurveLength,
                   weightMap);

        tempFGB(weightMap);

        for(ListDigraph::ArcIt a(tempFGB.graph());a!=INVALID;++a)
        {
            tempFGB.e
            arcsSet.insert(*a);
        }
    }



}

void computeFlow(SegCut::Image2D& image,
                 unsigned int gluedCurveLength,
                 std::map<Z2i::SCell,double>& weightMap,
                 Image2D& out)
{
    ImageFlowData imageFlowData(image);

    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,
                       gluedCurveLength);

    FlowGraphBuilder fgb(imageFlowData);
    singleFlow(fgb,
               imageFlowData,
               gluedCurveLength,
               weightMap);

    fgb(weightMap);

    do {

        FlowGraphQuery fgq(fgb);
        ListDigraph::ArcMap<bool> arcFilter(fgb.graph(),false);
        for (FlowGraphQuery::DetourArcMapIterator dami = fgq.detourArcsBegin(); dami != fgq.detourArcsEnd(); ++dami) {
            FlowGraphQuery::ArcPair key = dami->first;
            const std::vector<ListDigraph::Arc> &values = dami->second;

            for (FlowGraphQuery::DetourArcIterator dai = values.begin(); dai != values.end(); ++dai) {
                arcFilter[*dai] =true;
            }
        }

        ListDigraph::ArcMap<double> arcDiff(fgb.graph());
        updateAndCompare(image,
                         fgb,
                         imageFlowData,
                         arcDiff);


    }while();

}


void saveImage(Image2D& out,
               std::string outputFolder,
               std::string suffix)
{

    std::string imageOutputPath = outputFolder + "/" + suffix + ".pgm";

    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    GenericWriter<Image2D>::exportFile(imageOutputPath.c_str(),out);

}

void drawCurvatureMaps(Image2D& image,
                       std::map<Z2i::SCell,double>& weightMap,
                       int gluedCurveLength,
                       std::string outputFolder,
                       std::string suffix)
{
    ImageFlowData imageFlowData(image);
    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,gluedCurveLength);

    ImageFlowDataDebug imageFlowDataDebug(imageFlowData);

    imageFlowDataDebug.drawNoConnectionsCurvatureMap(weightMap,
                                                     outputFolder,
                                                     suffix);

    imageFlowDataDebug.drawInteriorConnectionsCurvatureMap(weightMap,
                                                           outputFolder,
                                                           suffix);

    imageFlowDataDebug.drawExteriorConnectionsCurvatureMap(weightMap,
                                                           outputFolder,
                                                           suffix);

    imageFlowDataDebug.drawMakeConvexConnectionsCurvatureMap(weightMap,
                                                             outputFolder,
                                                             suffix);

}


namespace Development{
    bool solveShift;
    bool crossElement;

    bool makeConvexArcs;
    bool invertGluedArcs;
};

int main(){
    Development::solveShift = false;
    Development::crossElement = false;

    Development::makeConvexArcs = false;
    Development::invertGluedArcs = true;

    unsigned int gluedCurveLength = 5;

    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/single_square.pgm");
    SegCut::Image2D imageOut = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/single_square.pgm");

    std::string outputFolder = "../output/flow/square/square-5-Invert";

    std::map<Z2i::SCell,double> weightMap;
    for(int i=0;i<200;++i)
    {
        weightMap.clear();


        computeFlow(image,
                    gluedCurveLength,
                    weightMap,
                    imageOut);


        drawCurvatureMaps(image,
                          weightMap,
                          gluedCurveLength,
                          outputFolder,
                          std::to_string(i));


        saveImage(imageOut,
                  outputFolder,
                  std::to_string(i));

        image=imageOut;

        std::cout << "OK " << i << std::endl;
    }

    return 0;
}