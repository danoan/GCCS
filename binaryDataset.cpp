#include <opencv/highgui.h>
#include <opencv/cv.hpp>

#include "utils/utils.h"
#include "utils/io.h"
#include "Flows/RefundFlow/RefundFlow.h"

using namespace UtilsTypes;

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};

namespace Configuration
{
    int gluedCurveLength = 5;
    std::string outputFolder = "";
    std::string inputFolder = "";
    ImageFlowData::FlowMode flowMode = ImageFlowData::DilationErosion;
}

void displayImage(std::string windowName,
                  std::string imagePath)
{
    cv::imshow( windowName,
                cv::imread(imagePath,CV_8U)
    );
    cv::waitKey(0);
}

bool noImprovement(double* currentEnergyValue)
{
    return ( fabs(currentEnergyValue[0]-currentEnergyValue[1])<0.00001 &&
            fabs(currentEnergyValue[1]-currentEnergyValue[2])<0.00001);
}

void preprocess(Image2D& image)
{
    Image2D outImage(image.domain());
    ImageProc::invertColors(outImage,outImage);

    ImageProc::resize(image,outImage);
    image = outImage;

    ImageProc::invertColors(image,image);
    ImageProc::closing(image,image,1);

    int borderWidth = 5;

    Domain newDomain( image.domain().lowerBound(),
                      image.domain().upperBound() + DGtal::Z2i::Point(2*borderWidth,2*borderWidth)
    );

    Image2D paddedImage(newDomain);
    ImageProc::createBorder(paddedImage,image,borderWidth);
    image = paddedImage;
}

void regularizeImage(std::string originalImagePath,
                     std::string outputFolder,
                     int maxIterations=50)
{
    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(originalImagePath);

    preprocess(image);
    IO::saveImage(image,outputFolder,"original");
    SegCut::Image2D imageOut = image;

    double currentEnergyValue[3];
    for(int i=0;i<maxIterations;++i)
    {
        ImageFlowData imageFlowData(image);

        imageFlowData.init(Configuration::flowMode,
                           Configuration::gluedCurveLength);


        RefundFlow refundFlow(imageFlowData);
        currentEnergyValue[i%3] = refundFlow.run(i);

        imageOut = refundFlow.outputImage();
        image=imageOut;

        IO::saveImage(imageOut,
                      outputFolder,
                      std::to_string(i)
        );

        std::cout << "(" << i << ") Energy Value: " << currentEnergyValue[i%3]
                  << std::endl;

        if( noImprovement(currentEnergyValue) ) break;
    }
}

int main()
{
    Configuration::flowMode = ImageFlowData::DilationErosion;
    Configuration::outputFolder = "../output/refundFlow/binary-set";
    Configuration::inputFolder = "../images/binary_images";

    IO::VectorOfPath vectorOfFiles;
    IO::listFiles(vectorOfFiles,Configuration::inputFolder);

    for(auto itp= vectorOfFiles.begin()+1;itp!=vectorOfFiles.end();++itp)
    {
        IO::Path p= *itp;
        std::cout << "Start Regularization Of Image: " << p.filename() << std::endl;
        regularizeImage(p.string(),
                        Configuration::outputFolder + "/" + p.filename().string() );
        std::cout << "End of Regularization" << std::endl;
    }

    return 0;
}