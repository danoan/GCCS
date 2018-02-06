#ifndef SEGCUT_UTILS
#define SEGCUT_UTILS

#include <iostream>

#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"

namespace Patch{
    extern bool cross_element;
};

namespace ImageProcTypes{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
    typedef DGtal::functors::SimpleThresholdForegroundPredicate<Image2D> ThreshPredicate;
};


void resize(ImageProcTypes::Image2D &input,ImageProcTypes::Image2D &out);

void computeBoundaryCurve(ImageProcTypes::Curve &boundCurve,
                          ImageProcTypes::KSpace &KImage,
                          const ImageProcTypes::Image2D &image,
                          unsigned int thresh_value=0);

void computeBoundaryCurve(ImageProcTypes::Curve &boundCurve,
                          ImageProcTypes::KSpace &KImage,
                          const ImageProcTypes::Image2D &image,
                          unsigned int thresh_value,
                          ImageProcTypes::Z2i::SCell imageBel);

void computeBoundaryCurve(ImageProcTypes::Curve &boundCurve,
                          ImageProcTypes::KSpace &KImage,
                          const ImageProcTypes::Image2D &image,
                          const ImageProcTypes::Image2D &mask,
                          unsigned int thresh_value=0);

void dilate(ImageProcTypes::Image2D &newImage, const ImageProcTypes::Image2D &inputImage, const int &dilation_size);
void dilate(ImageProcTypes::Image2D &newImage, const std::string &filepath, const int &dilation_size);

void dilateWithFilters(ImageProcTypes::Image2D &newImage, const std::string &filepath, const int &dilation_size);

void dilateWithMorphology(ImageProcTypes::Image2D &newImage, const std::string &filepath, const int &dilation_size);
void dilateWithMorphology(ImageProcTypes::Image2D &newImage, const ImageProcTypes::Image2D &inputImage, const int &dilation_size);
void dilateWithMorphology(ImageProcTypes::Image2D &newImage, const cv::Mat &src, const int &dilation_size);

void erode(ImageProcTypes::Image2D &newImage, const ImageProcTypes::Image2D &inputImage, const int &dilation_size);

void fromImage2DToMat(const ImageProcTypes::Image2D &dgtalImg,cv::Mat &cvImg);
void fromMatToImage2D(const cv::Mat &cvImg, ImageProcTypes::Image2D &dgtalImg,int shift=0);




#endif