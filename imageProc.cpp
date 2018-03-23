#include "imageProc.h"
#include "DGtal/io/writers/GenericWriter.h"

void ImageProc::fromMatToImage2D(const cv::Mat &cvImg, Image2D &dgtalImg, int shift){
    int ubY = cvImg.rows-1;
    for(int i=0;i<cvImg.rows;i++){
        for(int j=0;j<cvImg.cols;j++){
            unsigned char v( cvImg.at<unsigned char>(i,j) );
            dgtalImg.setValue( Point(j+shift,ubY-i+shift), v );
        }
    }
}

void ImageProc::fromImage2DToMat(const Image2D &dgtalImg, cv::Mat &cvImg){
    int ubY = dgtalImg.domain().upperBound()[1];

    for(auto it=dgtalImg.domain().begin();it!=dgtalImg.domain().end();++it){
        Point p = *it;
        unsigned char v( dgtalImg(*it) );
        cvImg.at<unsigned char>( (ubY- p[1]),p[0] ) = v;
    }
}

void ImageProc::dilateWithMorphology(Image2D &newImage, const std::string &filepath, const int &dilation_size){
    cv::Mat src = cv::imread(filepath.c_str(),cv::IMREAD_GRAYSCALE);
    dilateWithMorphology(newImage,src,dilation_size);
}

void ImageProc::dilateWithMorphology(Image2D &newImage, const Image2D &inputImage, const int &dilation_size){
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvSrc(r,c,CV_8UC1);
    fromImage2DToMat(inputImage,cvSrc);
    dilateWithMorphology(newImage,cvSrc,dilation_size);
}

void ImageProc::dilateWithMorphology(Image2D &newImage, const cv::Mat &src, const int &dilation_size){
    cv::Mat dilation_dst;
    int dilation_type;

    if(Development::crossElement){
        dilation_type = cv::MORPH_CROSS;
    }else{
        dilation_type = cv::MORPH_RECT;
    }


    cv::Mat element = cv::getStructuringElement( dilation_type,
                                                 cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
                                                 cv::Point( dilation_size, dilation_size ) );

    cv::dilate( src, dilation_dst, element );

    fromMatToImage2D(dilation_dst,newImage);
}

void ImageProc::dilate(Image2D &newImage, const std::string &filepath, const int &dilation_size){
    dilateWithMorphology(newImage,filepath,dilation_size);
}

void ImageProc::dilate(Image2D &newImage, const Image2D &inputImage, const int &dilation_size){
    dilateWithMorphology(newImage,inputImage,dilation_size);
}

void ImageProc::erode(Image2D &newImage, const Image2D &inputImage, const int &erosion_size){
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvSrc(r,c,CV_8UC1);
    fromImage2DToMat(inputImage,cvSrc);

    cv::Mat dilation_dst;
    int erosion_type;

    if(Development::crossElement){
        erosion_type = cv::MORPH_CROSS;
    }else{
        erosion_type = cv::MORPH_RECT;
    }

    cv::Mat element = cv::getStructuringElement( erosion_type,
                                                 cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                                                 cv::Point( erosion_size, erosion_size ) );

    cv::erode( cvSrc, dilation_dst, element );

    fromMatToImage2D(dilation_dst,newImage);
}

void ImageProc::opening(Image2D& newImage, const Image2D& inputImage, const int& element_size )
{
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvSrc(r,c,CV_8UC1);
    fromImage2DToMat(inputImage,cvSrc);

    cv::Mat erosion_dst;
    cv::Mat opening_dst;
    int element_type;

    if(Development::crossElement){
        element_type = cv::MORPH_CROSS;
    }else{
        element_type = cv::MORPH_RECT;
    }

    cv::Mat element = cv::getStructuringElement( element_type,
                                                 cv::Size( 2*element_size + 1, 2*element_size+1 ),
                                                 cv::Point( element_size, element_size ) );

    cv::erode( cvSrc, erosion_dst, element );
    cv::dilate( erosion_dst,opening_dst, element );

    fromMatToImage2D(opening_dst,newImage);
}

void ImageProc::closing(Image2D& newImage, const Image2D& inputImage, const int& element_size )
{
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvSrc(r,c,CV_8UC1);
    fromImage2DToMat(inputImage,cvSrc);

    cv::Mat dilation_dst;
    cv::Mat closing_dst;
    int closing_type;

    if(Development::crossElement){
        closing_type = cv::MORPH_CROSS;
    }else{
        closing_type = cv::MORPH_RECT;
    }

    cv::Mat element = cv::getStructuringElement( closing_type,
                                                 cv::Size( 2*element_size + 1, 2*element_size+1 ),
                                                 cv::Point( element_size, element_size ) );

    cv::dilate( cvSrc, dilation_dst, element );
    cv::erode( dilation_dst,closing_dst, element );

    fromMatToImage2D(closing_dst,newImage);
}


void ImageProc::resize(Image2D &input,Image2D &out, double factor)
{
    int rIn = input.domain().upperBound()[1] + 1;
    int cIn = input.domain().upperBound()[0] + 1;

    int rOut = out.domain().upperBound()[1] + 1;
    int cOut = out.domain().upperBound()[0] + 1;

    cv::Mat cvInput(rIn,cIn,CV_8UC1);
    cv::Mat cvOut(rOut,cOut,CV_8UC1);
    fromImage2DToMat(input,cvInput);
    cv::resize(cvInput,cvOut,cvOut.size(),factor,factor,cv::INTER_NEAREST);
    fromMatToImage2D(cvOut,out);
}

void ImageProc::computeBoundaryCurve(const Image2D &image,
                                     Curve &boundCurve,
                                     unsigned int thresh_value)
{
    Domain imageDomain = image.domain();
    KSpace KImage;

    KImage.init(imageDomain.lowerBound(),imageDomain.upperBound(),true);

    ThreshPredicate imagePredicate (image,thresh_value);
    SCell imageBel = Surfaces::findABel(KImage, imagePredicate, 10000);

    SurfelAdjacency SAdj(true);

    std::vector<SCell> boundarySCells;
    Surfaces::track2DBoundary(boundarySCells,
                              KImage,
                              SAdj,
                              imagePredicate,
                              imageBel);

    boundCurve.initFromSCellsVector(boundarySCells);

    eliminateLoops(boundCurve,KImage,boundCurve);
}

void ImageProc::computeBoundaryCurve(const Image2D& image,
                                     Curve& boundCurve,
                                     unsigned int thresh_value,
                                     SCell imageBel)
{
    Domain imageDomain = image.domain();
    KSpace KImage;

    KImage.init(imageDomain.lowerBound(),imageDomain.upperBound(),true);

    ThreshPredicate imagePredicate (image,thresh_value);

    SurfelAdjacency SAdj(true);

    std::vector<SCell> boundarySCells;
    Surfaces::track2DBoundary(boundarySCells,
                              KImage,
                              SAdj,
                              imagePredicate,
                              imageBel);

    boundCurve.initFromSCellsVector(boundarySCells);

}

void ImageProc::invertColors(Image2D& outputImage, Image2D& inputImage)
{
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvIn(r,c,CV_8UC1);

    fromImage2DToMat(inputImage,cvIn);
    auto begin = cvIn.begin<unsigned char>();
    auto end = cvIn.end<unsigned char>();
    for(auto it = begin;it!=end;++it)
    {
        *it = *it==255?0:255;
    }
    fromMatToImage2D(cvIn,outputImage);
}

void ImageProc::createBorder(Image2D& outputImage, Image2D& inputImage, int borderWidth)
{
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvIn(r,c,CV_8UC1);

    fromImage2DToMat(inputImage,cvIn);
    cv::copyMakeBorder(cvIn,cvIn,borderWidth,borderWidth,borderWidth,borderWidth,CV_HAL_BORDER_CONSTANT,0);

    fromMatToImage2D(cvIn,outputImage);
}