#include <string>

#include <boost/filesystem/path.hpp>

#include <opencv/highgui.h>
#include <opencv/cv.hpp>

#include <DGtal/io/readers/GenericReader.h>

#include "../typeDefs.h"
#include "../imageProc.h"
#include "../io.h"
#include "ProjectionFlow.h"
#include "AdaptativeScaleFlow.h"

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

int main() {
    unsigned int gluedCurveLength = 5;
//    std::string originalImagePath = "../images/flow-evolution/single_square_no_padding.pgm";
    std::string originalImagePath = "../images/flow-evolution/straight-line-bruit-small.pgm";
    std::string outputFolder = "../output/projectionFlow/linel/line-DE-MC-5";


    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(originalImagePath);
    SegCut::Image2D imageOut = image;

    std::string windowName = "adaptativeScaleFlow";
    cvNamedWindow(windowName.c_str(), CV_WINDOW_AUTOSIZE);


    displayImage(windowName, originalImagePath);
    ImageProc::closing(image, image, 1);
    displayImage(windowName, IO::saveImage(image,
                                           outputFolder,
                                           "closing")
    );

    int borderWidth = 5;

    Domain newDomain(image.domain().lowerBound(),
                     image.domain().upperBound() + DGtal::Z2i::Point(2 * borderWidth, 2 * borderWidth)
    );

    SegCut::Image2D paddedImage(newDomain);
    ImageProc::createBorder(paddedImage, image, borderWidth);
    image = paddedImage;

    AdaptativeScaleFlow adaptativeFlow(gluedCurveLength);
    for(int i=0;i<10;++i)
    {
        adaptativeFlow(image,imageOut);
        std::cout << i << " - Cut Value: " << adaptativeFlow.cutValue() << std::endl;

        IO::displayImage(windowName,IO::saveImage(imageOut,outputFolder,"iteration-" + std::to_string(i)));
        image = imageOut;
    }

}