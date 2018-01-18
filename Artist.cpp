#include "Artist.h"

void Artist::drawCurvatureMap(ArtistTypes::Curve& myCurve)
{
    double cmax=-100;
    double cmin=100;
    drawCurvatureMap(myCurve,cmin,cmax);
}

void Artist::drawCurvatureMap(ArtistTypes::Curve& myCurve,
                              double& cmin,
                              double& cmax)
{
    std::vector<double> estimations;
    curvatureEstimatorsGridCurve(myCurve.begin(),
                                 myCurve.end(),
                                 KImage,
                                 estimations);

    /*Max and Min value computation for Color Map*/
    if(squaredCurvature) {
        updateToSquared(estimations.begin(), estimations.end());
    }

    max_and_min(estimations,cmin,cmax);


    draw(estimations,
         myCurve.begin(),
         myCurve.end(),
         board,
         cmin,
         cmax);

}


void Artist::drawGluedCurvatureMap(SegCut::SCellGluedCurveIterator itb,
                                   SegCut::SCellGluedCurveIterator ite)
{
    double cmax=-100;
    double cmin=100;
    drawGluedCurvatureMap(itb,ite);
}

void Artist::drawGluedCurvatureMap(SegCut::SCellGluedCurveIterator itb,
                                   SegCut::SCellGluedCurveIterator ite,
                                   double& cmin,
                                   double& cmax)
{
    std::vector<double> estimations;
    curvatureEstimatorsGluedCurve(itb,
                                  ite,
                                  KImage,
                                  estimations);


    if(squaredCurvature){
        updateToSquared(estimations.begin(),estimations.end());
    }

    max_and_min(estimations,cmin,cmax);

    draw(estimations,
         itb,
         ite,
         board,
         cmin,
         cmax);
}


void Artist::drawCurvesAndConnectionsCurvatureMap(std::string imgFilePath,
                                                  std::string outputFilePath)
{
    Curve intCurvePriorGS,extCurvePriorGS;
    KSpace KImage;
    
    setCurves(imgFilePath,intCurvePriorGS,extCurvePriorGS);
    setKImage(imgFilePath,KImage);
    
    Curve intCurve,extCurve;
    Patch::initializeCurveCurvatureEstimator(KImage,
                                             intCurvePriorGS,
                                             intCurve);
    
    Patch::initializeCurveCurvatureEstimator(KImage,
                                             extCurvePriorGS,
                                             extCurve);    
    
    
    ArtistTypes::ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);


    unsigned int gluedCurveLength = 10;
    ArtistTypes::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    ArtistTypes::GluedCurveSetRange gluedCurveSetRange( seedRange.begin(),
                                                                seedRange.end(),
                                                                stgcF);
    double cmax=-100;
    double cmin=100;

    for(int i=0;i<2;++i){
        for(auto it=gluedCurveSetRange.begin();it!=gluedCurveSetRange.end();++it){
            drawGluedCurvatureMap(it->first,it->second,cmin,cmax );
        }

        drawCurvatureMap(intCurve,cmin,cmax);
        drawCurvatureMap(extCurve,cmin,cmax);
    }


    boost::filesystem::path p(outputFilePath.c_str());
    p.remove_filename();
    boost::filesystem::create_directories(p);

    board.saveEPS(outputFilePath.c_str());
};

void Artist::drawAllGluedCurves(std::string imgFilePath,
                                std::string outputFolder) {

    Curve intCurve, extCurve;

    setCurves(imgFilePath, intCurve, extCurve);
    ArtistTypes::ConnectorSeedRangeType seedRange = getSeedRange(KImage, intCurve, extCurve);


    unsigned int gluedCurveLength = 10;
    ArtistTypes::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    ArtistTypes::GluedCurveSetRange gluedCurveSetRange(seedRange.begin(),
                                                       seedRange.end(),
                                                       stgcF);
    double cmax = -100;
    double cmin = 100;

    drawCurvatureMap(intCurve, cmax, cmin);
    drawCurvatureMap(extCurve, cmax, cmin);
    for (auto it = gluedCurveSetRange.begin(); it != gluedCurveSetRange.end(); ++it) {
        drawGluedCurvatureMap(it->first, it->second, cmin, cmax);
    }


    int i = 0;
    std::string outputFilePath;
    boost::filesystem::create_directories(outputFolder);
    for (auto it = gluedCurveSetRange.begin(); it != gluedCurveSetRange.end(); ++it) {
        board.clear(DGtal::Color::White);

        drawGluedCurvatureMap(it->first, it->second, cmin, cmax);

        outputFilePath = outputFolder + "/" + std::to_string(i) + ".eps";
        board.saveEPS(outputFilePath.c_str());
        ++i;
    }


}

void Artist::drawTangentMap(Curve& myCurve,
                            double& cmin,
                            double& cmax)
{
    std::vector<ArtistTypes::TangentVector> estimations;

    tangentEstimatorsGridCurve(myCurve.begin(),
                               myCurve.end(),
                               KImage,
                               estimations);

    std::function< double( ArtistTypes::TangentVector ) > toDouble = []( ArtistTypes::TangentVector tv){return tv[0];};
    max_and_min(estimations,cmin,cmax,toDouble);


    draw(estimations,
         myCurve.begin(),
         myCurve.end(),
         board,
         cmin,
         cmax,
         toDouble);
}

void Artist::drawMaximalStabbingCircles(Curve& myCurve)
{
    typedef typename Curve::IncidentPointsRange::ConstCirculator AdapterCirculator;
    typedef StabbingCircleComputer<AdapterCirculator> SegmentComputer;
    typedef CurvatureFromDCAEstimator<SegmentComputer, false> SCFunctor;

    Curve::IncidentPointsRange intRange = myCurve.getIncidentPointsRange();

    typedef DGtal::SaturatedSegmentation<SegmentComputer> Segmentation;
    Segmentation seg( intRange.c(), intRange.c(), SegmentComputer() );


    board << DGtal::SetMode(SegmentComputer().className(), "Annulus");
    int n = 0;
    for (Segmentation::SegmentComputerIterator it = seg.begin(),
                 itEnd = seg.end(); it != itEnd; ++it, ++n)
    {
        board << (*it);
    }
}