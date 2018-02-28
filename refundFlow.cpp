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
#include "FlowGraph/weightSettings.h"

#include "FlowGraph/ImageFlowData.h"
#include "FlowGraph/FlowGraphBuilder.h"
#include "FlowGraph/ImageFlowDataDebug.h"
#include "FlowGraph/FlowGraphDebug.h"


using namespace UtilsTypes;

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
    bool solveShift;
    bool crossElement;

    bool makeConvexArcs;
    bool invertGluedArcs;
};

int main(){
    Development::solveShift = false;
    Development::crossElement = false;

    Development::makeConvexArcs = false;
    Development::invertGluedArcs = false;

    unsigned int gluedCurveLength = 8;

    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/single_triangle.pgm");
    SegCut::Image2D imageOut = image;

    std::string outputFolder = "../output/refundFlow/triangle/triangle-5";
    double currentEnergyValue;
    for(int i=0;i<200;++i)
    {
        ImageFlowData imageFlowData(image);
        imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,gluedCurveLength);

        RefundFlow refundFlow(imageFlowData);

        currentEnergyValue = refundFlow.run();

        drawCurvatureMaps(image,
                          gluedCurveLength,
                          outputFolder,
                          std::to_string(i));

        imageOut = refundFlow.outputImage();

        saveImage(imageOut,
                  outputFolder,
                  std::to_string(i));

        image=imageOut;

        std::cout << "OK " << i
                  << "    Energy Value: " << currentEnergyValue
                  << std::endl;

        //if(i==2) exit(1);
    }

    return 0;
}