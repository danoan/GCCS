

#include <string>
#include <DGtal/io/readers/GenericReader.h>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "typeDefs.h"
#include "imageProc.h"
#include "ProjectionFlow.h"

using namespace DGtal;

std::string saveImage(SegCut::Image2D& out,
                      std::string outputFolder,
                      std::string suffix)
{

    std::string imageOutputPath = outputFolder + "/" + suffix + ".pgm";

    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    GenericWriter<SegCut::Image2D>::exportFile(imageOutputPath.c_str(),out);
    return imageOutputPath;
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
    std::string originalImagePath = "../images/flow-evolution/single_square.pgm";
//    std::string originalImagePath = "../output/projectionFlow/linel/line-DE-MC-5/solution.pgm";
    std::string outputFolder = "../output/projectionFlow/linel/line-DE-MC-5";


    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(originalImagePath);
    SegCut::Image2D imageOut = image;

    std::string windowName = "projectionFlow";
    cvNamedWindow(windowName.c_str(),CV_WINDOW_AUTOSIZE);


    displayImage(windowName,originalImagePath);
    ImageProc::closing(image,image,1);
    displayImage(windowName,saveImage(image,
                                      outputFolder,
                                      "closing")
    );

    int borderWidth = 5;

    Domain newDomain( image.domain().lowerBound(),
                      image.domain().upperBound() + DGtal::Z2i::Point(2*borderWidth,2*borderWidth)
    );

    SegCut::Image2D paddedImage(newDomain);
    ImageProc::createBorder(paddedImage,image,borderWidth);
    image = paddedImage;
    for(int i=0;i<10;++i)
    {
        std::cout << "Iteration: " << i << std::endl;
        ImageFlowData imageFlowData(image);
        imageFlowData.init(ImageFlowData::DilationErosion, 5);

        ProjectionFlow PF(imageFlowData);
        PF(imageOut);

        displayImage(windowName, saveImage(imageOut, outputFolder, "solution"));
        cvWaitKey(0);

        image = imageOut;
    }
}