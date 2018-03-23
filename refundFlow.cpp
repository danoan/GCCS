#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;


#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include <boost/filesystem/path.hpp>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "utils.h"
#include "FlowGraph/weightSettings.h"

#include "FlowGraph/ImageFlowData.h"
#include "FlowGraph/FlowGraphBuilder.h"
#include "FlowGraph/ImageFlowDataDebug.h"
#include "FlowGraph/FlowGraphDebug.h"


using namespace UtilsTypes;

std::string saveImage(Image2D& out,
                      std::string outputFolder,
                      std::string suffix)
{

    std::string imageOutputPath = outputFolder + "/" + suffix + ".pgm";

    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    GenericWriter<Image2D>::exportFile(imageOutputPath.c_str(),out);
    return imageOutputPath;
}

void drawCurvatureMaps(Image2D& image,
                       int gluedCurveLength,
                       std::string outputFolder,
                       std::string suffix)
{
    ImageFlowData imageFlowData(image);
    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,gluedCurveLength);

    ImageFlowDataDebug imageFlowDataDebug(imageFlowData);
    FlowGraphBuilder::LinelWeightMap weightMap;
    setArcsWeight(imageFlowData,weightMap);

    FlowGraph flow;
    FlowGraphBuilder fgb(flow,imageFlowData,weightMap);

    imageFlowDataDebug.drawNoConnectionsCurvatureMap(weightMap,
                                                     outputFolder,
                                                     suffix);

//    imageFlowDataDebug.drawInteriorConnectionsCurvatureMap(weightMap,
//                                                           outputFolder,
//                                                           suffix);
//
//    imageFlowDataDebug.drawExteriorConnectionsCurvatureMap(weightMap,
//                                                           outputFolder,
//                                                           suffix);

//    imageFlowDataDebug.drawMakeConvexConnectionsCurvatureMap(weightMap,
//                                                             outputFolder,
//                                                             suffix);

}


namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};

void displayImage(std::string windowName,
                  std::string imagePath)
{
    cv::imshow( windowName,
                cv::imread(imagePath,CV_8U)
    );
    cv::waitKey(0);
}

int main()
{

    unsigned int gluedCurveLength = 5;
    std::string originalImagePath = "../images/flow-evolution/straight-line-bruit-small.pgm";
    std::string outputFolder = "../output/refundFlow/linel/line-DE-MC-5";


    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(originalImagePath);
    SegCut::Image2D imageOut = image;

    std::string windowName = "flowEvolution";
    cvNamedWindow(windowName.c_str(),CV_WINDOW_AUTOSIZE);


    displayImage(windowName,originalImagePath);
    ImageProc::closing(image,image,1);
    displayImage(windowName,saveImage(image,
                                      outputFolder,
                                      "closing")
    );


    double currentEnergyValue;
    for(int i=0;i<200;++i)
    {
        ImageFlowData imageFlowData(image);

        imageFlowData.init(ImageFlowData::FlowMode::DilationErosion,gluedCurveLength);


        RefundFlow refundFlow(imageFlowData);

        currentEnergyValue = refundFlow.run(i);

        drawCurvatureMaps(image,
                          gluedCurveLength,
                          outputFolder,
                          std::to_string(i));

        imageOut = refundFlow.outputImage();

        std::cout << "Press key to continue" << std::endl;
        displayImage(windowName,saveImage(imageOut,
                                          outputFolder,
                                          std::to_string(i))
        );


        image=imageOut;

        std::cout << "OK " << i
                  << "    Energy Value: " << currentEnergyValue
                  << std::endl;

//        exit(1);
    }

    return 0;
}