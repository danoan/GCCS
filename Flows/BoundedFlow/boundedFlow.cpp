#include <string>

#include <boost/filesystem/path.hpp>

#include <opencv/highgui.h>
#include <opencv/cv.hpp>

#include <DGtal/io/readers/GenericReader.h>

#include "../../utils/typeDefs.h"
#include "../../ImageProc/imageProc.h"
#include "../../utils/io.h"
#include "../ProjectionFlow/ProjectionFlow.h"
#include "../AdaptativeScaleFlow/AdaptativeScaleFlow.h"
#include "BoundedFlow.h"
#include "../PreprocessImage.h"

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;

    bool iteractive = false;
    std::string windowName = "boundedFlow";
};

void segmentImage(std::string originalImagePath,
                  std::string outputFolder,
                  int gluedCurveLength,
                  int maxIterations)
{
    ImageData ID(originalImagePath);
    SegCut::Image2D image = ID.preprocessedImage;
    SegCut::Image2D imageOut = image;

    std::string preprocessedFilepath = IO::saveImage(image,outputFolder,"preprocessing");

    if(Development::iteractive) {
        IO::displayImage(Development::windowName, originalImagePath);
        IO::displayImage(Development::windowName, preprocessedFilepath);
    }


    BoundedFlow boundedFlow(gluedCurveLength);

    for(int i=0;i<maxIterations;++i)
    {
        boundedFlow(image,imageOut);
        std::cout << i << " - Cut Value: " << boundedFlow.cutValue() << std::endl;

        std::string outputFilepath = IO::saveImage(imageOut,outputFolder,"iteration-" + std::to_string(i));
        if(Development::iteractive)
            IO::displayImage(Development::windowName,outputFilepath);

        image = imageOut;
    }
}

int main()
{
    cvNamedWindow(Development::windowName.c_str(), CV_WINDOW_AUTOSIZE);

    int gluedCurveLength = 5;
    std::string outputFolder = "../output/boundedFlow";
    std::string datasetFolder = "../images/segSet";


    typedef boost::filesystem::path path;
    typedef boost::filesystem::directory_iterator directory_iterator;

    path p( datasetFolder );
    for(directory_iterator it(p);it!=directory_iterator{};++it)
    {
        if( boost::filesystem::is_regular_file(*it) )
        {
            std::string filename = it->path().stem().generic_string();
            std::cout << "Segmentation of image:" << filename << std::endl;

            try {
                segmentImage(it->path().generic_string(), outputFolder + "/" + filename, gluedCurveLength, 5);
            }catch (Exception ex)
            {
                std::cout << "Segmentation could not be finished." << std::endl;
                std::cout << ex.what() << std::endl;
            }
        }
    }


}