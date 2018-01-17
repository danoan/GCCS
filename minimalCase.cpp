
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/StabbingCircleComputer.h"

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"

#include "typeDefs.h"
#include "imageProc.h"
#include "Patch/patch.h"

using namespace SegCut;

void computeCurvatureUsingPatch(std::string imgFilePath)
{
    Curve intCurvePriorGS;
    Image2D image = GenericReader<Image2D>::import(imgFilePath);
    KSpace KImage;

    computeBoundaryCurve(intCurvePriorGS,KImage,image);

    Curve intCurve;
    Patch::initializeCurveCurvatureEstimator(KImage,intCurvePriorGS,intCurve);

    DGtal::functors::SCellToIncidentPoints<KSpace> scelltopF(KImage);

    typedef ConstRangeAdapter< Curve::ConstIterator,
            DGtal::functors::SCellToIncidentPoints<KSpace>,
            std::pair<KSpace::Point,KSpace::Point> > MyRange;


    MyRange myIt(intCurve.begin(),
                 intCurve.end(),
                 scelltopF);

    MyRange::ConstCirculator begin = myIt.c();
    MyRange::ConstCirculator end = myIt.c();


    std::vector<double> estimations;
    Patch::estimationsPatchMCMSECurvature(begin,end,estimations);

    std::cout << estimations.size() << std::endl;
}

void computeCurvature(std::string imgFilePath)
{
    Curve intCurve;
    Image2D image = GenericReader<Image2D>::import(imgFilePath);
    KSpace KImage;

    computeBoundaryCurve(intCurve,KImage,image);

    DGtal::functors::SCellToIncidentPoints<KSpace> scelltopF(KImage);

    typedef ConstRangeAdapter< Curve::ConstIterator,
            DGtal::functors::SCellToIncidentPoints<KSpace>,
            std::pair<KSpace::Point,KSpace::Point> > MyRange;


    MyRange myIt(intCurve.begin(),
                 intCurve.end(),
                 scelltopF);

    MyRange::ConstCirculator begin = myIt.c();
    MyRange::ConstCirculator end = myIt.c();

    std::vector<double> estimations;
    Patch::estimationsDGtalMCMSECurvature(begin,end,estimations);

    std::cout << estimations.size() << std::endl;
}

namespace Patch
{
    bool useDGtal;
};

namespace UtilsTypes{
    std::function<double(double)> toDouble = [](double x){return x;};
};

int main()
{
    Patch::useDGtal = false;

    std::string imgFilePath = "../images/graph-weight-test/single_square.pgm";
//    std::string imgFilePath = "../images/graph-weight-test/single_triangle.pgm";
//    std::string imgFilePath = "../images/graph-weight-test/smallest_disk.pgm";
//    std::string imgFilePath = "../images/graph-weight-test/last_image.pgm";
    int i=100;
    do{
        computeCurvatureUsingPatch(imgFilePath);
        --i;
    }while(i>0);

}