#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "utils/utils.h"
#include "FlowGraph/weightSettings.h"

#include "FlowGraph/FlowGraphBuilder.h"
#include "FlowGraph/ImageFlowDataDebug.h"
#include "FlowGraph/FlowGraphDebug.h"
#include "utils/Artist.h"
#include "ImageProc/ImageProc.h"

using namespace UtilsTypes;

namespace Development{
    bool solveShift=false;
    bool crossElement=false;

    bool lambdaEstimator=false;
    bool pessimistEstimator=false;

    bool makeConvexArcs=false;
    bool invertGluedArcs=false;
};

void cellular_grid_space(std::string outputFolder)
{
    Board2D board;
    Domain domain( KSpace::Point(0,0), KSpace::Point(5,5));
    KSpace KImage;
    KImage.init(domain.lowerBound(),domain.upperBound(),true);

    KSpace::SCell pixel = KImage.sCell( KSpace::Point(5,5),true);
    KSpace::SCell linel = KImage.sCell( KSpace::Point(5,4),true);
    KSpace::SCell pointel = KImage.sCell( KSpace::Point(4,4),true);

    board << domain
          << DGtal::CustomStyle(pixel.className(), new DGtal::CustomColors(DGtal::Color::Blue, DGtal::Color::Blue))
          << pixel
          << DGtal::CustomStyle(pixel.className(), new DGtal::CustomColors(DGtal::Color::Green, DGtal::Color::Green))
          << linel
          << DGtal::CustomStyle(pixel.className(), new DGtal::CustomColors(DGtal::Color::Red, DGtal::Color::Red))
          << pointel;

    board.saveEPS( (outputFolder + "/cellular_grid_space.eps").c_str() );
}

void grid_curve(std::string outputFolder)
{
    Image2D imageOriginal = GenericReader<Image2D>::import("../images/flow-evolution/small_square.pgm");
    Domain domain = imageOriginal.domain();
    Board2D board;

    DigitalSet digitalImage(domain);
    ImageProc2::ImageAsDigitalSet(digitalImage,imageOriginal);

    Curve boundaryCurve;
    ImageProc::computeBoundaryCurve(imageOriginal,boundaryCurve,100);

    board << domain
          << DGtal::CustomStyle(boundaryCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Green, DGtal::Color::Green))
          << boundaryCurve;

    board.saveEPS( (outputFolder + "/grid_curve.eps").c_str() );
}

void glued_curve(std::string outputFolder)
{
    std::string imgFilePath = "../images/flow-evolution/small_square.pgm";
    Image2D imageOriginal = GenericReader<Image2D>::import(imgFilePath);
    Domain domain = imageOriginal.domain();
    Curve intCurve, extCurve;

    KSpace KImage;
    KImage.init(domain.lowerBound(),domain.upperBound(),true);

    setCurves(imgFilePath, intCurve, extCurve);
    ArtistTypes::ConnectorSeedRangeType seedRange = getSeedRange(KImage, intCurve, extCurve);

    Board2D board;

    unsigned int gluedCurveLength = 5;
    ArtistTypes::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    ArtistTypes::GluedCurveSetRange gluedCurveSetRange(seedRange.begin(),
                                                       seedRange.end(),
                                                       stgcF);
    board << domain
          << DGtal::CustomStyle(intCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Green, DGtal::Color::Green))
          << intCurve
          << DGtal::CustomStyle(intCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Blue, DGtal::Color::Blue))
          << extCurve;

    board << DGtal::CustomStyle(intCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Red, DGtal::Color::Red));
    GluedCurveSetRange::ConstIterator it = gluedCurveSetRange.begin();

    SegCut::SCellGluedCurveIterator itBegin = it->first;
    SegCut::SCellGluedCurveIterator itEnd = it->second;
    for(auto it = itBegin;it!=itEnd;++it)
    {
        board << *it;
    }

    board.saveEPS( (outputFolder + "/glued_curve.eps").c_str() );
}


void dual_graph(std::string outputFolder)
{
    std::string imgFilePath = "../images/flow-evolution/small_square.pgm";
    Image2D imageOriginal = GenericReader<Image2D>::import(imgFilePath);
    Domain domain = imageOriginal.domain();
    Curve intCurve, extCurve;

    KSpace KImage;
    KImage.init(domain.lowerBound(),domain.upperBound(),true);

    setCurves(imgFilePath, intCurve, extCurve);
    ArtistTypes::ConnectorSeedRangeType seedRange = getSeedRange(KImage, intCurve, extCurve);

    Board2D board;

    unsigned int gluedCurveLength = 5;
    ArtistTypes::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    ArtistTypes::GluedCurveSetRange gluedCurveSetRange(seedRange.begin(),
                                                       seedRange.end(),
                                                       stgcF);
    board << domain
          << DGtal::CustomStyle(intCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Green, DGtal::Color::Green))
          << intCurve
          << DGtal::CustomStyle(intCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Blue, DGtal::Color::Blue))
          << extCurve;

    board << DGtal::CustomStyle(intCurve.begin()->className(), new DGtal::CustomColors(DGtal::Color::Red, DGtal::Color::Red));
    for(GluedCurveSetRange::ConstIterator it = gluedCurveSetRange.begin();it!=gluedCurveSetRange.end();++it)
    {
        board << *it->first.connectorsBegin();
    }

    board.saveEPS( (outputFolder + "/primal_graph.eps").c_str() );



    ImageFlowData imf(imageOriginal);
    imf.init(ImageFlowData::DilationOnly,gluedCurveLength);

    FlowGraphBuilder::LinelWeightMap weightMap;
    setArcsWeight(imf,weightMap);

    FlowGraph fg;
    FlowGraphBuilder fgb(fg,imf,weightMap);
    FlowGraphDebug fgd(fg);
    fgd.drawFlowGraph(outputFolder,"dual_graph_intext");

}


int main()
{
    std::string outputFolder = "../output/reportFigures";
    boost::filesystem::path p2(outputFolder.c_str());
    boost::filesystem::create_directories(p2);

    cellular_grid_space(outputFolder);
    grid_curve(outputFolder);
    glued_curve(outputFolder);
    dual_graph(outputFolder);

    /*
    Image2D imageOriginal = GenericReader<Image2D>::import("../images/comedic-mi-parcours/line/line-artifacts.pgm");
    Image2D imageIteration = GenericReader<Image2D>::import("../images/comedic-mi-parcours/line/8.pgm");

    int gluedCurveLenght = 5;

    ImageFlowData imageFlowData(imageOriginal);
    imageFlowData.init(ImageFlowData::DilationOnly,gluedCurveLenght);

    std::string outputFolder = "../output/comedic-mi-parcours";


    Board2D board;
    board << imageFlowData.getMostInnerCurve();
    board << imageFlowData.getMostOuterCurve();

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        for(auto itR=it->gcsRangeBegin();itR!=it->gcsRangeEnd();++itR)
        {
            board << *itR->first.connectorsBegin();
        }
    }
    board.save( (outputFolder+"/base-grid.eps").c_str() );





    std::map<DGtal::Z2i::SCell,double> weightMap;
    setArcsWeight(imageFlowData,weightMap);
    FlowGraph fg;
    FlowGraphBuilder fgb(fg,imageFlowData,weightMap);

    FlowGraphDebug fgd(fg);
    ListDigraph::ArcMap<bool> sourceArcs(fg.graph(),false);
    ListDigraph::ArcMap<bool> targetArcs(fg.graph(),false);
    ListDigraph::ArcMap<bool> notTerminalArcs(fg.graph(),true);
    FlowGraphQuery::filterArcs(fg,sourceArcs,FlowGraph::ArcType::SourceArc,true);
    FlowGraphQuery::filterArcs(fg,sourceArcs,FlowGraph::ArcType::TargetArc,true);

    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> fc(fg.graph(),notTerminalArcs);
    fc-sourceArcs-targetArcs;

    fgd.highlightArcs(fc(),outputFolder,"flow-graph");


    ListDigraph::ArcMap<bool> cutFilter(fg.graph(),false);
    FlowGraphQuery::sourceComponent(fg,cutFilter);
    ListDigraph::ArcMap<int> arcColors(fg.graph(),0);
    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(cutFilter[a])
        {
            arcColors[a]=1;
        }
    }

    ListDigraph::NodeMap<bool> allNodes(fg.graph(),true);
    fgd.highlightArcs(allNodes,fc(),arcColors,outputFolder,"highlight-cut-Arcs");


    ImageFlowData imf2(imageIteration);
    imf2.init(ImageFlowData::DilationOnly,gluedCurveLenght);

    board.clear();
    board << imf2.getMostInnerCurve();
    board.save("../output/comedic-mi-parcours/iteration.eps");


    {
        ImageFlowData imf(imageOriginal);
        imf.init(ImageFlowData::DilationOnly,5);

        Artist EA(imf.getKSpace(), board);
        //Maximal Stabbing Circles
        EA.board.clear(DGtal::Color::White);
        EA.drawMaximalStabbingCircles(imf.getMostInnerCurve());
        std::string outputFilePath = outputFolder + "/Original-MaximalStabbingCircles.eps";
        EA.board.save(outputFilePath.c_str());
    }

    {
        ImageFlowData imf(imageIteration);
        imf.init(ImageFlowData::DilationOnly,5);

        Artist EA(imf.getKSpace(), board);
        //Maximal Stabbing Circles
        EA.board.clear(DGtal::Color::White);
        EA.drawMaximalStabbingCircles(imf.getMostInnerCurve());
        std::string outputFilePath = outputFolder + "/Iteration-MaximalStabbingCircles.eps";
        EA.board.save(outputFilePath.c_str());
    }
    */

    return 0;
}