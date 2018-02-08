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
#include "FlowGraphBuilder.h"
#include "ImageFlowData.h"
#include "ImageFlowDataDebug.h"
#include "FlowGraphDebug.h"

using namespace UtilsTypes;

void tangentWeight(Curve::ConstIterator begin,
                   Curve::ConstIterator end,
                   KSpace& KImage,
                   std::vector< TangentVector >& estimationsTangent,
                   std::vector< double >& tangentWeightVector)
{
    auto it = begin;
    int i =0;
    do{
        UtilsTypes::KSpace::Point pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
        UtilsTypes::KSpace::Point pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );

        UtilsTypes::KSpace::Point scellVector = pTarget-pSource;

        tangentWeightVector.push_back( fabs( estimationsTangent[i].dot(scellVector) ) );
        ++it;
        ++i;
    }while(it!=end);
}

void setGridCurveWeight(Curve curvePriorGS,
                        KSpace& KImage,
                        std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> curvatureEstimations;
    curvatureEstimatorsGridCurve(curvePriorGS.begin(),
                                 curvePriorGS.end(),
                                 KImage,
                                 curvatureEstimations);


    updateToSquared(curvatureEstimations.begin(),curvatureEstimations.end());


    std::vector<TangentVector> tangentEstimations;
    tangentEstimatorsGridCurve(curvePriorGS.begin(),
                               curvePriorGS.end(),
                               KImage,
                               tangentEstimations);


    std::vector<double> tangentWeightVector;
    tangentWeight(curvePriorGS.begin(),
                  curvePriorGS.end(),
                  KImage,
                  tangentEstimations,
                  tangentWeightVector);


    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] = curvatureEstimations[i];
            ++i;
        }
    }

    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] *=tangentWeightVector[i];
            ++i;
        }
    }

}

void setGluedCurveWeight(GluedCurveSetRange::ConstIterator gcsRangeBegin,
                         GluedCurveSetRange::ConstIterator gcsRangeEnd,
                         KSpace& KImage,
                         unsigned int gluedCurveLength,
                         std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> estimationsCurvature;
    curvatureEstimatorsConnections(gcsRangeBegin,gcsRangeEnd,KImage,gluedCurveLength,estimationsCurvature);

    updateToSquared(estimationsCurvature.begin(),estimationsCurvature.end());

    {
        int i = 0;
        for (GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {

            auto itC = it->first.connectorsBegin();
            do{
                weightMap[*itC] = estimationsCurvature[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);

        }
    }

    std::vector<TangentVector> estimationsTangent;
    tangentEstimatorsConnections(gcsRangeBegin,gcsRangeEnd,KImage,gluedCurveLength,estimationsTangent);


    std::vector<double> tangentWeightVector;
    {
        GluedCurveIteratorPair it = gcsRangeBegin;
        int i = 0;
        for (GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {

            auto itC = it->first.connectorsBegin();
            do {

                KSpace::SCell linel = *itC;

                UtilsTypes::KSpace::Point pTarget = KImage.sCoords(KImage.sDirectIncident(linel, *KImage.sDirs(linel)));
                UtilsTypes::KSpace::Point pSource = KImage.sCoords(KImage.sIndirectIncident(linel, *KImage.sDirs(linel)));

                UtilsTypes::KSpace::Point scellVector = pTarget - pSource;

                tangentWeightVector.push_back(fabs(estimationsTangent[i].dot(scellVector)));

                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);

        }

    }

    {
        int i = 0;
        for (GluedCurveIteratorPair it = gcsRangeBegin; it != gcsRangeEnd; ++it) {
            auto itC = it->first.connectorsBegin();
            do {
                weightMap[*itC]*= tangentWeightVector[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);
        }
    }

}


void computeFlow(SegCut::Image2D& image,
                 unsigned int gluedCurveLength,
                 std::map<Z2i::SCell,double>& weightMap,
                 std::vector<Z2i::Point>& coordPixelsSourceSide,
                 std::vector<Z2i::SCell>& pixelsInTheGraph,
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


    FlowGraphDebug flowGraphDebug(fgb);
    flowGraphDebug.drawCutGraph(outputFolder,suffix );

}


void updateImage(std::vector<Z2i::Point>& coordPixelsSourceSide,
                 std::vector<Z2i::SCell>& pixelsInTheGraph,
                 Image2D& out,
                 std::string outputFolder,
                 std::string suffix)
{
    KSpace KImage;
    KImage.init(out.domain().lowerBound(),out.domain().upperBound(),true);
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

    std::string outputFolder = "../output/flow/square/square-5-Invert";

    std::vector<Z2i::Point> coordPixelsSourceSide;
    std::vector<Z2i::SCell> pixelsInTheGraph;
    std::map<Z2i::SCell,double> weightMap;

    for(int i=0;i<200;++i)
    {
        coordPixelsSourceSide.clear();
        pixelsInTheGraph.clear();
        weightMap.clear();


        computeFlow(image,
                    gluedCurveLength,
                    weightMap,
                    coordPixelsSourceSide,
                    pixelsInTheGraph,
                    outputFolder,
                    std::to_string(i));


        drawCurvatureMaps(image,
                          weightMap,
                          gluedCurveLength,
                          outputFolder,
                          std::to_string(i));


        updateImage(coordPixelsSourceSide,
                    pixelsInTheGraph,
                    image,
                    outputFolder,
                    std::to_string(i));



        std::cout << "OK " << i << std::endl;
    }

    return 0;
}