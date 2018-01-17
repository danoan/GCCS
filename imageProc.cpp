#include "imageProc.h"

void fromMatToImage2D(const cv::Mat &cvImg, ImageProcTypes::Image2D &dgtalImg){
    int ubY = cvImg.rows-1;
    for(int i=0;i<cvImg.rows;i++){
        for(int j=0;j<cvImg.cols;j++){
            unsigned char v( cvImg.at<unsigned char>(i,j) );
            dgtalImg.setValue( ImageProcTypes::Point(j,ubY-i), v );
        }
    }
}

void fromImage2DToMat(const ImageProcTypes::Image2D &dgtalImg, cv::Mat &cvImg){
    int ubY = dgtalImg.domain().upperBound()[1];

    for(auto it=dgtalImg.domain().begin();it!=dgtalImg.domain().end();++it){
        ImageProcTypes::Point p = *it;
        unsigned char v( dgtalImg(*it) );
        cvImg.at<unsigned char>( (ubY- p[1]),p[0] ) = v;
    }
}

void dilateWithMorphology(ImageProcTypes::Image2D &newImage, const std::string &filepath, const int &dilation_size){
    cv::Mat src = cv::imread(filepath.c_str(),cv::IMREAD_GRAYSCALE);
    dilateWithMorphology(newImage,src,dilation_size);
}

void dilateWithMorphology(ImageProcTypes::Image2D &newImage, const ImageProcTypes::Image2D &inputImage, const int &dilation_size){
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvSrc(r,c,CV_8UC1);
    fromImage2DToMat(inputImage,cvSrc);
    dilateWithMorphology(newImage,cvSrc,dilation_size);
}

void dilateWithMorphology(ImageProcTypes::Image2D &newImage, const cv::Mat &src, const int &dilation_size){
    cv::Mat dilation_dst;
    int dilation_type = cv::MORPH_RECT;

    cv::Mat element = cv::getStructuringElement( dilation_type,
                                                 cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
                                                 cv::Point( dilation_size, dilation_size ) );

    dilate( src, dilation_dst, element );

    fromMatToImage2D(dilation_dst,newImage);
}

void dilateWithFilters(ImageProcTypes::Image2D &newImage, const std::string &filepath, const int &dilation_size){
    cv::Mat src,binary_src;
    int dilation_type = cv::MORPH_RECT;

    src = cv::imread(filepath.c_str(),cv::IMREAD_GRAYSCALE);
    cv::threshold(src,binary_src,0,1,cv::THRESH_BINARY);

    //Bottom
    cv::Mat ker_bottom = (cv::Mat_<char>(3,3) << 0,1,0,
            0,1,0,
            0,0,0);

    //Top
    cv::Mat ker_top = (cv::Mat_<char>(3,3) << 0,0,0,
            0,1,0,
            0,1,0);

    //Left
    cv::Mat ker_left = (cv::Mat_<char>(3,3) << 0,0,0,
            0,1,1,
            0,0,0);

    //Right
    cv::Mat ker_right= (cv::Mat_<char>(3,3) << 0,0,0,
            1,1,0,
            0,0,0);

    cv::Mat binary_temp,binary_extended;
    cv::filter2D(binary_src,binary_temp,binary_src.depth(),ker_bottom);
    cv::filter2D(binary_temp,binary_extended,binary_src.depth(),ker_top);
    cv::filter2D(binary_extended,binary_temp,binary_src.depth(),ker_left);
    cv::filter2D(binary_temp,binary_extended,binary_src.depth(),ker_right);


    int r = binary_extended.rows;
    int c = binary_extended.cols;

    std::cout << src.rows << " ; " << r << std::endl;

    uchar pixval;
    int maxY = newImage.domain().upperBound()[1];

    for(int i=0;i<r;i++) {
        for(int j=0;j<c;j++){
            unsigned char pixval = binary_extended.at<unsigned char>(i, j);
            newImage.setValue( ImageProcTypes::Point(j,maxY-i),pixval);
        }
    }
}

void dilate(ImageProcTypes::Image2D &newImage, const std::string &filepath, const int &dilation_size){
    dilateWithMorphology(newImage,filepath,dilation_size);
}

void dilate(ImageProcTypes::Image2D &newImage, const ImageProcTypes::Image2D &inputImage, const int &dilation_size){
    dilateWithMorphology(newImage,inputImage,dilation_size);
}

void erode(ImageProcTypes::Image2D &newImage, const ImageProcTypes::Image2D &inputImage, const int &erosion_size){
    int r = inputImage.domain().upperBound()[1] + 1;
    int c = inputImage.domain().upperBound()[0] + 1;

    cv::Mat cvSrc(r,c,CV_8UC1);
    fromImage2DToMat(inputImage,cvSrc);

    cv::Mat dilation_dst;
    int erosion_type = cv::MORPH_RECT;

    cv::Mat element = cv::getStructuringElement( erosion_type,
                                                 cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                                                 cv::Point( erosion_size, erosion_size ) );

    erode( cvSrc, dilation_dst, element );

    fromMatToImage2D(dilation_dst,newImage);
}

void computeBoundaryCurve(ImageProcTypes::Curve &boundCurve,
                          ImageProcTypes::KSpace &KImage,
                          const ImageProcTypes::Image2D &image,
                          unsigned int thresh_value){
    ImageProcTypes::Domain imageDomain = image.domain();

    ImageProcTypes::SurfelAdjacency<ImageProcTypes::KSpace::dimension> SAdj(true);
    KImage.init(imageDomain.lowerBound(),imageDomain.upperBound(),true);

    ImageProcTypes::ThreshPredicate imagePredicate (image,thresh_value);
    ImageProcTypes::SCell imageBel = ImageProcTypes::Surfaces<ImageProcTypes::KSpace>::findABel(KImage, imagePredicate, 10000);

    std::vector<ImageProcTypes::SCell> boundarySCells;
    ImageProcTypes::Surfaces<ImageProcTypes::KSpace>::track2DBoundary(boundarySCells,
                                             KImage,
                                             SAdj,
                                             imagePredicate,
                                             imageBel);

    boundCurve.initFromSCellsVector(boundarySCells);


//    If gridCurve is constructed using the code below, I lost the ordering.
//    ThreshSurface threshSurf(KImage,imagePredicate,SAdj,imageBel);
//    for(ThreshSurface::SurfelConstIterator it = threshSurf.begin(); it!=threshSurf.end();it++) {
//        boundCurve.push_back(*it);
//    }
}

void computeBoundaryCurve(ImageProcTypes::Curve &boundCurve,
                          ImageProcTypes::KSpace &KImage,
                          const ImageProcTypes::Image2D &image,
                          const ImageProcTypes::Image2D &mask,
                          unsigned int thresh_value)
{

    ImageProcTypes::Domain imageDomain = image.domain();

    ImageProcTypes::SurfelAdjacency<ImageProcTypes::KSpace::dimension> SAdj(true);
    KImage.init(imageDomain.lowerBound(),imageDomain.upperBound(),true);

    ImageProcTypes::ThreshPredicate imagePredicate (mask,thresh_value);
    ImageProcTypes::SCell imageBel = ImageProcTypes::Surfaces<ImageProcTypes::KSpace>::findABel(KImage, imagePredicate, 10000);

    std::vector<ImageProcTypes::SCell> boundarySCells;
    ImageProcTypes::Surfaces<ImageProcTypes::KSpace>::track2DBoundary(boundarySCells,
                                             KImage,
                                             SAdj,
                                             imagePredicate,
                                             imageBel);

    boundCurve.initFromSCellsVector(boundarySCells);
}