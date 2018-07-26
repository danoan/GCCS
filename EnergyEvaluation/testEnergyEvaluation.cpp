#include "EnergyEvaluation.h"

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};

int main()
{
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
    Image2D image = DGtal::GenericReader<Image2D>::import("../images/img.pgm");

    EnergyEvaluation::DigitalSet ds(image.domain());
    ImageProc::ImageAsDigitalSet(ds,image);

    std::cout << EnergyEvaluation::integratedSquaredCurvature(ds) << std::endl;

    return 0;
}