#ifndef SEGBYCUT_EXPERIMENTARTIST_H
#define SEGBYCUT_EXPERIMENTARTIST_H

#include "DGtal/io/boards/Board2D.h"
#include "boost/filesystem.hpp"

#include "typeDefs.h"
#include "imageProc.h"
#include "Patch/patch.h"
#include "utils.h"

namespace ArtistTypes{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    using namespace SegCut;
}

class Artist
{
public:
    ArtistTypes::KSpace KImage;
    ArtistTypes::Board2D board;
    bool squaredCurvature;

    Artist(ArtistTypes::KSpace& KImage,
           ArtistTypes::Board2D& board):KImage(KImage),board(board),squaredCurvature(false){}

    void setOptional(bool squaredCurvature=false){this->squaredCurvature =squaredCurvature;}

    void drawCurvatureMap(ArtistTypes::Curve& myCurve);

    void drawCurvatureMap(ArtistTypes::Curve& myCurve,
                          double& cmin,
                          double& cmax);


    void drawGluedCurvatureMap(ArtistTypes::SCellGluedCurveIterator itb,
                               ArtistTypes::SCellGluedCurveIterator ite);

    void drawGluedCurvatureMap(ArtistTypes::SCellGluedCurveIterator itb,
                               ArtistTypes::SCellGluedCurveIterator ite,
                               double& cmin,
                               double& cmax);


    void drawCurvesAndConnectionsCurvatureMap(std::string imgFilePath,
                                              std::string outputFilePath);



    void drawAllGluedCurves(std::string imgFilePath,
                            std::string outputFolder);


    void drawTangentMap(ArtistTypes::Curve& myCurve,
                        double& cmin,
                        double& cmax);

    void drawMaximalStabbingCircles(ArtistTypes::Curve& myCurve);
};


#endif //SEGBYCUT_EXPERIMENTARTIST_H
