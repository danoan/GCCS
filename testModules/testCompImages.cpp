#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;


#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "../utils/utils.h"
#include "../FlowGraph/weightSettings.h"

#include "../FlowGraph/FlowGraphBuilder.h"
#include "../FlowGraph/FlowGraphQuery.h"
#include "../FlowGraph/FlowGraphDebug.h"
#include "../utils/Artist.h"

using namespace UtilsTypes;

double computeEnergyValue(Image2D& image,
                          ImageFlowData& model,
                          std::string name)
{
    ImageFlowData imf(image);
    imf.init(model.getFlowMode(),model.getGluedCurveLength());

    std::map<Z2i::SCell,double> weightMap;
    setArcsWeight(imf,weightMap);


    FlowGraph f;
    FlowGraphBuilder fgb(f,imf,weightMap);


    ListDigraph::ArcMap<bool> internalArcs(f.graph(),false);
    FlowGraphQuery::filterArcs(f,
                               internalArcs,
                               FlowGraph::ArcType::InternalCurveArc,
                               true);


    double s = 0;
    for(ListDigraph::ArcIt a(f.graph());a!=INVALID;++a)
    {
        if(internalArcs[a]){
            s+= f.weight(a);
        }
    }


    return s;

}

void constructSCellsVectorFromFilter(std::vector<Z2i::SCell>& scells,
                                     FlowGraph& fg,
                                     ListDigraph::ArcMap<bool>& arcFilter)
{
    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(arcFilter[a])
        {
            if(fg.arcType(a)==FlowGraph::IntExtGluedArc ||
               fg.arcType(a)==FlowGraph::ExtIntGluedArc ||
               fg.arcType(a)==FlowGraph::InternalCurveArc ||
               fg.arcType(a)==FlowGraph::ExternalCurveArc )
            {
                scells.push_back(fg.scell(a));
            }
        }
    }
}


void internalLinels(Image2D& image,
                    std::vector<Z2i::SCell>& internalLinels,
                    std::map<Z2i::SCell,double>& weightMap)
{
    ImageFlowData imf(image);
    imf.init(ImageFlowData::DilationOnly,5);

    setArcsWeight(imf,weightMap);


    FlowGraph f;
    FlowGraphBuilder fgb(f,imf,weightMap);


    ListDigraph::ArcMap<bool> internalArcs(f.graph(),false);
    FlowGraphQuery::filterArcs(f,
                               internalArcs,
                               FlowGraph::ArcType::InternalCurveArc,
                               true);

    constructSCellsVectorFromFilter(internalLinels,f,internalArcs);
}

void drawCurvatureMaps(Image2D& comp1,
                       Image2D& comp2,
                       std::string outputFolder)
{
    std::vector<Z2i::SCell> internalLinelsComp1;
    std::map<Z2i::SCell,double> weightMapComp1;
    internalLinels(comp1,internalLinelsComp1,weightMapComp1);

    std::vector<Z2i::SCell> internalLinelsComp2;
    std::map<Z2i::SCell,double> weightMapComp2;
    internalLinels(comp2,internalLinelsComp2,weightMapComp2);

    std::function<double(Z2i::SCell)> fnTempToDoubleComp1 = [weightMapComp1](Z2i::SCell s){ return weightMapComp1.at(s); };
    std::function<double(Z2i::SCell)> fnTempToDoubleComp2 = [weightMapComp2](Z2i::SCell s){ return weightMapComp2.at(s); };
    double cmin=100;double cmax=-100;
    max_and_min(internalLinelsComp1,cmin,cmax,fnTempToDoubleComp1);
    max_and_min(internalLinelsComp2,cmin,cmax,fnTempToDoubleComp2);


    Board2D board;
    draw(weightMapComp1,internalLinelsComp1,board,cmin,cmax);
    board.saveEPS( std::string( (outputFolder + "/comp1.eps") ).c_str() );

    board.clear();

    draw(weightMapComp2,internalLinelsComp2,board,cmin,cmax);
    board.saveEPS( std::string( (outputFolder + "/comp2.eps") ).c_str() );
}

void drawStabbingCircles(Image2D& image,
                         std::string outputFolder,
                         std::string suffix)
{
    ImageFlowData imf(image);
    imf.init(ImageFlowData::DilationOnly,5);

    Board2D board;
    Artist EA(imf.getKSpace(),board);

    EA.drawMaximalStabbingCircles( imf.getMostInnerCurve() );
    EA.board.save( std::string(outputFolder + "/" + suffix + ".eps").c_str() );

}

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};

void evaluate_energy_results()
{
    typedef boost::filesystem::path path;
    typedef boost::filesystem::directory_iterator directory_iterator;

    std::string datasetFolder = "../images/energy-digital-circle";

    path p( datasetFolder );
    int i =0;
    for(directory_iterator it(p);it!=directory_iterator{};++it)
    {
            if (boost::filesystem::is_regular_file(*it))
            {
                std::string filename = it->path().stem().generic_string();

                SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(it->path().generic_string());
                ImageFlowData imf(image);
                imf.init(ImageFlowData::DilationOnly, 5);
                double energyValue = computeEnergyValue(image, imf, std::to_string(i));

                std::cout << filename << ":" << energyValue << std::endl;
                ++i;
            }
    }
}

int main()
{
    evaluate_energy_results();
    return 0;

    SegCut::Image2D comp1 = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/comp1.pgm");
    SegCut::Image2D comp2 = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/comp2.pgm");

    int gluedCurveLength = 5;
    std::string outputFolder = "../output/testModules/testCompImages";

    boost::filesystem::path p2(outputFolder.c_str());
    boost::filesystem::create_directories(p2);

    ImageFlowData imf1(comp1);
    imf1.init(ImageFlowData::DilationOnly,gluedCurveLength);
    double energyComp1 = computeEnergyValue(comp1,imf1,"1");

    ImageFlowData imf2(comp2);
    imf2.init(ImageFlowData::DilationOnly,gluedCurveLength);
    double energyComp2 = computeEnergyValue(comp2,imf2,"1");

    std::cout <<  energyComp1 << std::endl;
    std::cout << energyComp2 << std::endl;

    drawCurvatureMaps(comp1,comp2,outputFolder);
    drawStabbingCircles(comp1,outputFolder,"stabbingComp1");
    drawStabbingCircles(comp2,outputFolder,"stabbingComp2");

    drawCurvatureMaps(comp1,comp2,outputFolder);
    drawStabbingCircles(comp1,outputFolder,"stabbingComp1");
    drawStabbingCircles(comp2,outputFolder,"stabbingComp2");


    return 0;
}