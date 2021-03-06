#include <string>

#include <boost/filesystem/path.hpp>

#include <opencv/highgui.h>
#include <opencv/cv.hpp>

#include <DGtal/io/readers/GenericReader.h>

#include "../../utils/typeDefs.h"
#include "../../ImageProc/imageProc.h"
#include "../../utils/io.h"
#include "../ProjectionFlow/ProjectionFlow.h"
#include "AdaptativeScaleFlow.h"
#include "../PreprocessImage.h"

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;

    bool iteractive = false;
    std::string windowName = "adaptativeScaleFlow";
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


    AdaptativeScaleFlow adaptativeFlow(gluedCurveLength);
    for(int i=0;i<maxIterations;++i)
    {
        adaptativeFlow(image,imageOut);
        std::cout << i << " - Cut Value: " << adaptativeFlow.cutValue() << std::endl;

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
    std::string outputFolder = "../output/adaptativeScaleFlow";
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
            }catch(...)
            {
                std::cout << "Segmentation could not be finished." << std::endl;
            }
        }
    }


}