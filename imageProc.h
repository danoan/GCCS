#ifndef SEGCUT_UTILS
#define SEGCUT_UTILS

#include <iostream>

#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"

namespace Development{
    extern bool crossElement;
}

namespace ImageProc {
    typedef DGtal::Z2i::Curve Curve;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::SurfelAdjacency<KSpace::dimension> SurfelAdjacency;
    typedef DGtal::Surfaces<KSpace> Surfaces;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
    typedef DGtal::functors::SimpleThresholdForegroundPredicate<Image2D> ThreshPredicate;


    void resize(Image2D &input, Image2D &out);


    void computeBoundaryCurve(const Image2D& image,
                              Curve& boundCurve,
                              unsigned int thresh_value);

    void computeBoundaryCurve(const Image2D& image,
                              Curve& boundCurve,
                              unsigned int thresh_value,
                              SCell imageBel);

    void dilate(Image2D &newImage, const Image2D &inputImage, const int &dilation_size);

    void dilate(Image2D &newImage, const std::string &filepath, const int &dilation_size);

    void dilateWithMorphology(Image2D &newImage, const std::string &filepath, const int &dilation_size);

    void dilateWithMorphology(Image2D &newImage, const Image2D &inputImage,
                              const int &dilation_size);

    void dilateWithMorphology(Image2D &newImage, const cv::Mat &src, const int &dilation_size);

    void erode(Image2D &newImage, const Image2D &inputImage, const int &dilation_size);

    void closing(Image2D& newImage, const Image2D& inputImage, const int& element_size );

    void fromImage2DToMat(const Image2D &dgtalImg, cv::Mat &cvImg);

    void fromMatToImage2D(const cv::Mat &cvImg, Image2D &dgtalImg, int shift = 0);
}


#endif