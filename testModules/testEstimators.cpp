
#include "../utils/typeDefs.h"
#include "../utils/utils.h"
#include "../utils/Artist.h"

using namespace DGtal;
using namespace DGtal::Z2i;

using namespace SegCut;

/*Tests if the generated GluedCurves are connected*/
void testConnectdeness(std::string imgFilePath)
{
    KSpace KImage;
    setKImage(imgFilePath,KImage);

    Curve intCurve,extCurve;
    setCurves(imgFilePath,intCurve,extCurve);

    ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);
    unsigned int gluedCurveLength = 10;
    SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    GluedCurveSetRange gcsRange( seedRange.begin(),
                                         seedRange.end(),
                                         stgcF);


    DGtal::functors::SCellToIncidentPoints<KSpace> myFunctor(KImage);
    int gluedCurveNumber =0;
    for(auto itgc=gcsRange.begin();itgc!=gcsRange.end();++itgc ){
//        std::cout << "Glued Curve Start " << gluedCurveNumber << std::endl;

        SCellGluedCurveIterator begin = itgc->first;
        SCellGluedCurveIterator end = itgc->second;

        if(begin.connectorType()==makeConvex) continue;

        GluedCurveIncidentPointsRange gcipRange(begin,
                                                end,
                                                myFunctor );

        GridCurve<KSpace> gc;
        gc.initFromSCellsRange(begin,end);
        std::cout << gc << std::endl;
        ++gluedCurveNumber;
    }

    Board2D board;
    board << intCurve;
    board << extCurve;
    board.save("original-dilation.eps");


}

/*Evaluates curvature for each generated GluedCurve*/
void testCurvatureEvaluation(std::string imgFilePath)
{
    KSpace KImage;
    setKImage(imgFilePath,KImage);

    Curve intCurve,extCurve;
    setCurves(imgFilePath,intCurve,extCurve);



    ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);
    unsigned int gluedCurveLength = 10;
    SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);


    DGtal::functors::SCellToIncidentPoints<KSpace> myFunctor(KImage);
    int gluedCurveNumber=0;
    for(auto itgc=gcsRange.begin();itgc!=gcsRange.end();++itgc ){
//        std::cout << "Glued Curve Start " << gluedCurveNumber << std::endl;

        SCellGluedCurveIterator begin = itgc->first;
        SCellGluedCurveIterator end = itgc->second;
        
        if(begin.connectorType()==makeConvex) continue;

        GluedCurveIncidentPointsRange gcipRange(begin,
                                                end,
                                                myFunctor);

        std::vector<double> estimations;

        curvatureEstimatorsGluedCurve(begin,
                                      end,
                                      KImage,
                                      estimations);

        {
            auto it = gcipRange.begin();
            int i =0;
            do {
                ++it;
                ++i;
            } while (it != gcipRange.end());
        }

        ++gluedCurveNumber;
    }
}


void runTests(const std::string& imgPath,
              const std::string& outputFolder)
{
    testConnectdeness(imgPath);
    testCurvatureEvaluation(imgPath);



    KSpace KImage;
    setKImage(imgPath,KImage);

    Board2D board;
    Artist EA(KImage,board);

    EA.setOptional(true);

    Curve intCurve,extCurve;
    setCurves(imgPath,intCurve,extCurve);

    /*All Glued Curves*/
    EA.board.clear(DGtal::Color::White);
    EA.drawAllGluedCurves(imgPath,outputFolder+"/gluedCurves");


    /*Complete Curvature Map*/
    EA.board.clear(DGtal::Color::White);
    std::string outputFilePath = outputFolder + "/completeCurvatureMap.eps";
    EA.drawCurvesAndConnectionsCurvatureMap(imgPath,outputFilePath);


    /*Curvature Maps*/
    EA.board.clear(DGtal::Color::White);
    outputFilePath = outputFolder + "/curvatureMap.eps";
    {
        int i=2;
        double cmax=-100;
        double cmin=100;
        do{
            EA.drawCurvatureMap(intCurve,cmin,cmax);
            EA.drawCurvatureMap(extCurve,cmin,cmax);
            --i;
        }while(i>0);
    }
    EA.board.saveEPS(outputFilePath.c_str());



    /*Tangent Maps*/
    EA.board.clear(DGtal::Color::White);
    outputFilePath = outputFolder + "/tangentMap.eps";
    {
        int i=2;
        double cmax=-100;
        double cmin=100;
        do{
            EA.drawTangentMap(intCurve,cmin,cmax);
            EA.drawTangentMap(extCurve,cmin,cmax);
            --i;
        }while(i>0);
    }
    EA.board.saveEPS(outputFilePath.c_str());

    /*Maximal Stabbing Circles*/
    EA.board.clear(DGtal::Color::White);
    EA.drawMaximalStabbingCircles(intCurve);
    outputFilePath = outputFolder + "/maximalStabbingCircles.eps";
    EA.board.save(outputFilePath .c_str());

}

void testSequence()
{
    std::string outputRootPath = "../output/testModules/testEstimators";
    if(Development::solveShift)
        outputRootPath += "/solve-shift+average";
    else
        outputRootPath += "/average";
    
//    runTests("../images/graph-weight-test/last_image.pgm",
//             outputRootPath + "/last_image");
//
//    runTests("../images/graph-weight-test/single_square.pgm",
//                 outputRootPath + "/square");
//
//    runTests("../images/graph-weight-test/single_triangle.pgm",
//             outputRootPath + "/triangle");
//
//    runTests("../images/graph-weight-test/smallest_disk.pgm",
//             outputRootPath + "/disk");

    runTests("../images/flow-evolution/test.pgm",
             outputRootPath + "/test");
}

namespace Development{
    bool solveShift;
    bool crossElement;
    bool lambdaEstimator;
    bool pessimistEstimator=false;

    bool makeConvexArcs;
    bool invertGluedArcs;
};


int main()
{
    Development::crossElement = false;

    Development::makeConvexArcs = false;
    Development::invertGluedArcs = false;

    Development::lambdaEstimator = false;

//    Development::solveShift = true;
//    testSequence();

    Development::solveShift = false;
    testSequence();

}

