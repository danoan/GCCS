
#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "utils.h"
#include "FlowGraph/weightSettings.h"

#include "Artist.h"
#include "FlowGraph/FlowGraphBuilder.h"
#include "FlowGraph/ImageFlowDataDebug.h"
#include "FlowGraph/FlowGraphDebug.h"

using namespace UtilsTypes;


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

        SCellGluedCurveIterator begin = itgc->first;
        SCellGluedCurveIterator end = itgc->second;

        if(begin.connectorType()==makeConvex){
            std::cout << "make-convex-glued-curve" << std::endl;
        }

        GluedCurveIncidentPointsRange gcipRange(begin,
                                                end,
                                                myFunctor );

        GridCurve<KSpace> gc;
        gc.initFromSCellsRange(begin,end);
        ++gluedCurveNumber;
    }

}

void computeFlow(SegCut::Image2D& image,
                 unsigned int gluedCurveLength,
                 std::map<Z2i::SCell,double>& weightMap,
                 std::string outputFolder,
                 std::string suffix)
{
    ImageFlowData imageFlowData(image);

    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,
                       gluedCurveLength);


    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        setGridCurveWeight(it->curve,
                           imageFlowData.getKSpace(),
                           weightMap
        );
    }

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        setGluedCurveWeight( it->gcsRangeBegin(),
                             it->gcsRangeEnd(),
                             imageFlowData.getKSpace(),
                             gluedCurveLength,
                             weightMap);
    }

    FlowGraphBuilder fgb(imageFlowData);
    fgb(weightMap);

    FlowGraphDebug flowGraphDebug(fgb);
    flowGraphDebug.drawCutGraph(outputFolder,suffix);
}

void drawCurvatureMaps(Image2D& image,
                       std::map<Z2i::SCell,double>& weightMap,
                       int gluedCurveLength,
                       std::string outputFolder,
                       std::string suffix)
{
    ImageFlowData imageFlowData(image);
    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,gluedCurveLength);

    ImageFlowDataDebug imageFlowDataDebug(imageFlowData);

    imageFlowDataDebug.drawNoConnectionsCurvatureMap(weightMap,
                                                     outputFolder,
                                                     suffix);

    imageFlowDataDebug.drawInteriorConnectionsCurvatureMap(weightMap,
                                                           outputFolder,
                                                           suffix);

    imageFlowDataDebug.drawExteriorConnectionsCurvatureMap(weightMap,
                                                           outputFolder,
                                                           suffix);

    imageFlowDataDebug.drawMakeConvexConnectionsCurvatureMap(weightMap,
                                                             outputFolder,
                                                             suffix);

}

void updateImage(std::vector<Z2i::Point>& coordPixelsSourceSide,
                 std::vector<Z2i::SCell>& pixelsInTheGraph,
                 Image2D& out,
                 std::string& imageOutputPath)
{
    KSpace KImage;
    KImage.init(out.domain().lowerBound(),out.domain().upperBound(),true);
    for(auto it=coordPixelsSourceSide.begin();it!=coordPixelsSourceSide.end();++it)
    {
        out.setValue(*it,255);
    }


    /*Fill the holes*/
    for(auto it=pixelsInTheGraph.begin();it!=pixelsInTheGraph.end();++it)
    {
        Z2i::SCell pixel = *it;
        Z2i::Point p = KImage.sCoords(pixel);
        Z2i::SCells N = KImage.sNeighborhood(pixel);

        unsigned char v = out(p);
        unsigned char nv;
        bool hole=true;
        for(int j=1;j<5;++j){
            nv = out(KImage.sCoords(N[j]));
            if(v==nv){
                hole=false;
                break;
            }
        }
        if(hole){
            out.setValue(p,nv);
        }
    }


    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    GenericWriter<Image2D>::exportFile(imageOutputPath.c_str(),out);

}

void particularCurve(std::string outputFolder)
{
    Curve c;
    KSpace KImage;
    KImage.init( Z2i::Point(0,0),Z2i::Point(50,50),true);

    c.push_back( KImage.sCell( Z2i::Point(75,76),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(74,75),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(73,74),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(72,73),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(71,72),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(70,71),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(69,70),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(68,69),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(68,67),KSpace::POS) );

    c.push_back( KImage.sCell( Z2i::Point(67,66),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(66,67),KSpace::NEG) );
    c.push_back( KImage.sCell( Z2i::Point(65,68),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(64,67),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(64,65),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(64,63),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(63,62),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(62,61),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(61,60),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(60,59),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(60,57),KSpace::POS) );
    c.push_back( KImage.sCell( Z2i::Point(60,55),KSpace::POS) );

    c.isClosed();
    c.isValid();

    Board2D board;

    std::vector<double> estimations;
    curvatureEstimatorsGridCurve(c.begin(),c.end(),KImage,estimations,false);
    updateToSquared(estimations.begin(),estimations.end());

    double cmin=100;
    double cmax=-100;
    max_and_min(estimations,cmin,cmax);
    draw(estimations,c.begin(),c.end(),board,cmin,cmax);

    board.saveEPS( (outputFolder + "/particularCurve.eps").c_str() );

}

namespace Development{
    bool solveShift;
    bool crossElement;

    bool makeConvexArcs;
    bool invertGluedArcs;
};

int main(){
    Development::solveShift = false;
    Development::crossElement = false;

    Development::makeConvexArcs = false;
    Development::invertGluedArcs = false;

    unsigned int gluedCurveLength = 7;

    std::string imgPath = "../images/flow-evolution/out20.pgm";
    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(imgPath);

    std::string outputFolder = "../output/issue/out20";


    KSpace KImage;
    setKImage("../images/flow-evolution/out20.pgm",KImage);

    Board2D board;
    Artist EA(KImage,board);

    EA.setOptional(true);

    Curve intCurve,extCurve;
    setCurves(imgPath,intCurve,extCurve);

    testConnectdeness(imgPath);

    /*All Glued Curves*/
    EA.board.clear(DGtal::Color::White);
    EA.drawAllGluedCurves(imgPath,outputFolder+"/gluedCurves");


    std::map<Z2i::SCell,double> weightMap;
    computeFlow(image,
                gluedCurveLength,
                weightMap,
                outputFolder,
                "");


    drawCurvatureMaps(image,
                      weightMap,
                      gluedCurveLength,
                      outputFolder,
                      std::to_string(1));

    particularCurve(outputFolder);



    return 0;
}