

#include <string>
#include <DGtal/io/readers/GenericReader.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "../../utils/typeDefs.h"
#include "../../ImageProc/imageProc.h"
#include "ProjectionFlow.h"
#include "../PreprocessImage.h"

using namespace DGtal;

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;

    bool iteractive = false;
    std::string windowName = "projectionFlow";
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



    for(int i=0;i<maxIterations;++i)
    {
        ImageFlowData imageFlowData(image);
        imageFlowData.init(ImageFlowData::DilationErosion, gluedCurveLength);

        ProjectionFlow PF(imageFlowData);
        PF(imageOut);

        std::cout << i << " - Cut Value: " << PF.cutValue() << std::endl;

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
    std::string outputFolder = "../output/projectionFlow";
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
