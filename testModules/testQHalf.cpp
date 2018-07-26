#include "../ImageProc/ImageProc.h"

namespace Development{
    bool solveShift=false;
    bool crossElement=false;

    bool lambdaEstimator=false;
    bool pessimistEstimator=false;

    bool makeConvexArcs=false;
    bool invertGluedArcs=false;
};

int main()
{
    typedef ImageProc::ImageAsDigitalSet::Domain Domain;
    typedef ImageProc::ImageAsDigitalSet::Image2D Image2D;
    typedef ImageProc::ImageAsDigitalSet::DigitalSet DigitalSet;
    typedef ImageProc::DigitalBallIntersection::EuclideanBall EuclideanBall;
    typedef ImageProc::DigitalBallIntersection::Point Point;

    typedef ImageProc::DigitalBoundary< EightNeighborhoodPredicate<DigitalSet>  > EightBoundary;



    Image2D image = DGtal::GenericReader<Image2D>::import("../images/segSet/single_square.pgm");
    Domain domain( Point(-50,-50), Point(50,50));

    DigitalSet ds(domain);
    ImageProc::DigitalRectangle(ds,domain,Point(-25,-25),50,50);

    DigitalSet boundary(domain);
    EightBoundary(boundary,ds);

    Point center = *boundary.begin();
    unsigned long int radius = 3;
    double grid_size = 1;

    EuclideanBall eb(center, radius);
    DigitalSet db(domain);
    DGtal::Shapes<Domain>::euclideanShaper(db, eb, grid_size);



    DigitalSet dsQHalf(domain);
    ImageProc::QHalf(dsQHalf,ds,radius,center);

    std::cout << dsQHalf.size() << std::endl;

    DGtal::Board2D board;

    board <<  DGtal::CustomStyle(ds.className(),
                                 new DGtal::CustomColors(DGtal::Color::Gray,
                                                         DGtal::Color::Gray))
          <<  ds
          <<  DGtal::CustomStyle(ds.className(),
                                 new DGtal::CustomColors(DGtal::Color::Blue,
                                                         DGtal::Color::Blue))
          << boundary
          << DGtal::CustomStyle(ds.className(),
                                new DGtal::CustomColors(DGtal::Color::Red,
                                                        DGtal::Color::Red))
          << db
          << DGtal::CustomStyle(ds.className(),
                                new DGtal::CustomColors(DGtal::Color::Yellow,
                                                        DGtal::Color::Yellow))
          << dsQHalf;

    board.saveEPS("qhalf.eps");

    return 0;
}