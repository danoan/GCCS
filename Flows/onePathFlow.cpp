#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;


#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include <boost/filesystem/path.hpp>
#include <SeedCapture/CapturedData.h>
#include <SeedCapture/seedCapture.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"


#include "../utils/utils.h"
#include "../FlowGraph/weightSettings.h"

#include "../FlowGraph/Patch/OnePathFlowGraphBuilder.h"
#include "../FlowGraph/ImageFlowData.h"
#include "../FlowGraph/ImageFlowDataDebug.h"
#include "../FlowGraph/FlowGraphDebug.h"
#include "../utils/io.h"
#include "PreprocessImage.h"
#include "../utils/Artist.h"


namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;

    bool iteractive = false;
    std::string windowName = "simpleFlow";
};

void computeFlow(SegCut::Image2D& image,
                 SegCut::Image2D& refImage,
                 unsigned int gluedCurveLength,
                 std::map<Z2i::SCell,double>& weightMap,
                 std::vector<Z2i::Point>& coordPixelsSourceSide,
                 std::vector<Z2i::SCell>& pixelsInTheGraph,
                 std::string outputFolder,
                 std::string suffix,
                 ImageFlowData::FlowMode flowMode,
                 double* distrFrg,
                 double* distrBkg)
{
    ImageFlowData imageFlowData(image,refImage);

    imageFlowData.init(flowMode,
                       gluedCurveLength);

    setArcsWeight(imageFlowData,weightMap);

    FlowGraph fg;
    Development::OnePathFlowGraphBuilder fgb(fg,imageFlowData,weightMap,distrFrg,distrBkg);


    FlowGraph::FlowComputer flow = fg.prepareFlow();
    flow.run();

    ListDigraph::NodeMap<bool> node_filter(fg.graph(),false);
    ListDigraph::ArcMap<bool> arc_filter(fg.graph(),true);

    for(ListDigraph::NodeIt n(fg.graph());n!=INVALID;++n){
        if( flow.minCut(n) ){
            node_filter[n] = true;
        }else{
            node_filter[n] = false;
        }
    }

    node_filter[fg.source()]=false;



    FlowGraphQuery::SubGraph subgraph(fg.graph(),
                                      node_filter,
                                      arc_filter);


    KSpace& KImage = imageFlowData.getKSpace();

    for(FlowGraphQuery::SubGraph::NodeIt n(subgraph);n!=INVALID;++n)
    {
        Z2i::SCell pixel = fg.pixel(n);
        Z2i::Point p = KImage.sCoords(pixel);

        coordPixelsSourceSide.push_back(p);
    }

    for(ListDigraph::NodeIt n(fg.graph());n!=INVALID;++n)
    {
        if(fg.source()==n || fg.target()==n) continue;
        Z2i::SCell pixel = fg.pixel(n);

        pixelsInTheGraph.push_back(pixel);
    }


    FlowGraphDebug flowGraphDebug(fg);
    flowGraphDebug.drawCutGraph(outputFolder,suffix);
    flowGraphDebug.drawFlowGraph(outputFolder,suffix);

}


void updateImage(std::vector<Z2i::Point>& coordPixelsSourceSide,
                 std::vector<Z2i::SCell>& pixelsInTheGraph,
                 SegCut::Image2D& out,
                 std::string outputFolder,
                 std::string suffix)
{
    KSpace KImage;
    KImage.init(out.domain().lowerBound(),out.domain().upperBound(),true);


    SegCut::Image2D temp = out;
    ImageProc::Curve boundary;
    ImageProc::computeBoundaryCurve(temp,boundary,100);

    for(auto it=boundary.getInnerPointsRange().begin();it!=boundary.getInnerPointsRange().end();++it)
    {
        out.setValue(*it,0);
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

    std::string imageOutputPath = outputFolder + "/" + suffix + ".pgm";

    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    GenericWriter<SegCut::Image2D>::exportFile(imageOutputPath.c_str(),out);

//    Board2D board;
//    Artist EA(KImage,board);
//    EA.drawCurvesAndConnectionsCurvatureMap(imageOutputPath,outputFolder + "/curvmap" + suffix + ".eps");

}

void segmentImage(std::string originalImagePath,
                  std::string outputFolder,
                  int gluedCurveLength,
                  int maxIterations)
{
    typedef Z2i::Point Point;

    ImageData ID(originalImagePath);
    SegCut::Image2D refImage = ID.preprocessedImage;
    SegCut::Image2D imageOut = refImage;

    std::string preprocessedFilepath = IO::saveImage(refImage,outputFolder,"preprocessing");

    if(Development::iteractive) {
        IO::displayImage(Development::windowName, originalImagePath);
        IO::displayImage(Development::windowName, preprocessedFilepath);
    }


    std::vector<Z2i::Point> coordPixelsSourceSide;
    std::vector<Z2i::SCell> pixelsInTheGraph;
    std::map<Z2i::SCell,double> weightMap;

    SeedCapture::CapturedData cd = SeedCapture::capture(preprocessedFilepath);
    std::cout << cd.foregroundRegion.size() << ";" << std::endl;

    SegCut::Image2D seedImage(refImage.domain());
    for(auto it=refImage.domain().begin();it!=refImage.domain().end();++it)
    {
        seedImage.setValue(*it,0);
    }

    for(auto it=cd.foregroundRegion.begin();it!=cd.foregroundRegion.end();++it)
    {
        seedImage.setValue(Point( it->x,refImage.domain().upperBound()[1] - it->y),255);
    }


    ImageFlowData::FlowMode fm;
    for(int i=0;i<maxIterations;++i)
    {
        coordPixelsSourceSide.clear();
        pixelsInTheGraph.clear();
        weightMap.clear();

        //if(i%2!=0) fm = ImageFlowData::FlowMode::DilationOnly;
        //else fm = ImageFlowData::FlowMode::ErosionOnly;

        fm = ImageFlowData::FlowMode::DilationOnly;

        computeFlow(seedImage,
                    refImage,
                    gluedCurveLength,
                    weightMap,
                    coordPixelsSourceSide,
                    pixelsInTheGraph,
                    outputFolder,
                    std::to_string(i),
                    fm,
                    cd.foregroundDistrib,
                    cd.backgroundDistrib);


        updateImage(coordPixelsSourceSide,
                    pixelsInTheGraph,
                    seedImage,
                    outputFolder,
                    std::to_string(i));

        std::cout << "OK " << i << std::endl;
    }

}

int main()
{
    //cvNamedWindow(Development::windowName.c_str(), CV_WINDOW_AUTOSIZE);

    int gluedCurveLength = 5;
    std::string outputFolder = "../output/onepath_flow/binary_images";
    std::string datasetFolder = "../images/binary_images";


    typedef boost::filesystem::path path;
    typedef boost::filesystem::directory_iterator directory_iterator;

    path p( datasetFolder );
    for(directory_iterator it(p);it!=directory_iterator{};++it)
    {
        if( boost::filesystem::is_regular_file(*it) )
        {
            std::string filename = it->path().stem().generic_string();
            std::cout << "Segmentation of image:" << filename << std::endl;

            if(filename!="low_contrast") continue;

            try {
                segmentImage(it->path().generic_string(), outputFolder + "/binary_images_NOC/" + filename, gluedCurveLength, 200);
            }catch (Exception ex)
            {
                std::cout << "Segmentation could not be finished." << std::endl;
                std::cout << ex.what() << std::endl;
            }
        }
    }


}