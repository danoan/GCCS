#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;


#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include <boost/filesystem/path.hpp>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "utils.h"
#include "weightSettings.h"

#include "ImageFlowData.h"
#include "FlowGraphBuilder.h"
#include "ImageFlowDataDebug.h"
#include "FlowGraphDebug.h"


using namespace UtilsTypes;



double computeEnergyValue(Image2D& image,
                          ImageFlowData& model)
{
    ImageFlowData imf(image);
    imf.init(model.getFlowMode(),model.getGluedCurveLength());



    Flow f(imf);
    std::map<Z2i::SCell,double>& weightMap = f.getWeightMap();

    ListDigraph::ArcMap<bool> internalArcs(f.graph(),false);
    FlowGraphQuery::filterArcs(f.graphBuilder(),
                               internalArcs,
                               FlowGraphBuilder::ArcType::InternalCurveArc,
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


void updateAndCompare(Flow& f1,
                      Image2D& out,
                      ListDigraph::ArcMap<bool>& arcFilter,
                      ListDigraph::ArcMap<double>& arcDiff)
{
    f1.updateImage(out);

    ImageFlowData imageFlowData(out);
    f1.initImageFlowDataAsInFlow(imageFlowData);

    Flow f2(imageFlowData);

    ListDigraph::ArcMap<int> arcsMatch(f1.graph());
    Flow::matchFlows(arcsMatch,f1,f2);

    for(ListDigraph::ArcIt arcF1(f1.graph());arcF1!=INVALID;++arcF1)
    {
        if(arcFilter[arcF1])
        {
            ListDigraph::Arc arcF2 = f2.graph().arcFromId( arcsMatch[arcF1] );
            arcDiff[arcF1] = f2.weight(arcF2) - f1.weight(arcF1);
        }
    }

}

void computeFlow(SegCut::Image2D& image,
                 unsigned int gluedCurveLength,
                 Image2D& out,
                 std::string outputFolder,
                 int mainIteration)
{
    ImageFlowData imageFlowData(image);

    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,
                       gluedCurveLength);

    Flow f1(imageFlowData);
    double previousCutValue = f1.cutValue();
    double currentCutValue = previousCutValue*10;
    int iteration = 1;

    Image2D partialImage = f1.baseImage();
    Image2D previousPartialImage = f1.baseImage();
    while( true ) {

        FlowGraphDebug fgd(f1);

        ListDigraph::ArcMap<bool> detourArcs(f1.graph(),false);
        f1.detourArcsFilter(detourArcs);

        ListDigraph::ArcMap<double> arcDiff(f1.graph());

        updateAndCompare(f1,
                         partialImage,
                         detourArcs,
                         arcDiff);

        if(!f1.hasChanges(partialImage,previousPartialImage))
            break;

        previousPartialImage = partialImage;

        double s;
        int n;
        for(auto it = f1.detourArcsMapBegin();it!=f1.detourArcsMapEnd();++it)
        {
            FlowGraphQuery::ArcPair key = it->first;
            const std::vector<ListDigraph::Arc> &values = it->second;

            s=0;
            n=0;
            for(auto vit=values.begin();vit!=values.end();++vit)
            {
                s+= arcDiff[*vit];
                ++n;
            }

            s/=(2*n);

            f1.addRefundArcs(key.first,key.second,s);
        }


        std::string suffix = "RefundArcs-" + std::to_string(mainIteration) + "-" + std::to_string(iteration);
        fgd.drawCutGraph(outputFolder,suffix);
        fgd.drawRefundGraph(outputFolder,suffix);
        fgd.highlightArcs(detourArcs,
                          outputFolder,
                          "DetourArcs-" + std::to_string(mainIteration) + "-" + std::to_string(iteration) );


        previousCutValue = currentCutValue;
        currentCutValue = f1.cutValue();

        iteration++;
    };

    std::cout<< "iterations: " << iteration << std::endl;

    f1.updateImage(out);
    double actualEnergyValue = computeEnergyValue(out,imageFlowData);

    std::cout << "Cut Value: " << currentCutValue << std::endl;
    std::cout << "Energy Value: " << actualEnergyValue << std::endl;


}


void saveImage(Image2D& out,
               std::string outputFolder,
               std::string suffix)
{

    std::string imageOutputPath = outputFolder + "/" + suffix + ".pgm";

    boost::filesystem::path p2(imageOutputPath.c_str());
    p2.remove_filename();
    boost::filesystem::create_directories(p2);

    GenericWriter<Image2D>::exportFile(imageOutputPath.c_str(),out);

}

void drawCurvatureMaps(Image2D& image,
                       int gluedCurveLength,
                       std::string outputFolder,
                       std::string suffix)
{
    ImageFlowData imageFlowData(image);
    imageFlowData.init(ImageFlowData::FlowMode::DilationOnly,gluedCurveLength);

    ImageFlowDataDebug imageFlowDataDebug(imageFlowData);

    Flow flow(imageFlowData);
    std::map<Z2i::SCell,double>& weightMap = flow.getWeightMap();

    imageFlowDataDebug.drawNoConnectionsCurvatureMap(weightMap,
                                                     outputFolder,
                                                     suffix);

//    imageFlowDataDebug.drawInteriorConnectionsCurvatureMap(weightMap,
//                                                           outputFolder,
//                                                           suffix);
//
//    imageFlowDataDebug.drawExteriorConnectionsCurvatureMap(weightMap,
//                                                           outputFolder,
//                                                           suffix);

//    imageFlowDataDebug.drawMakeConvexConnectionsCurvatureMap(weightMap,
//                                                             outputFolder,
//                                                             suffix);

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

    unsigned int gluedCurveLength = 20;

    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/single_square.pgm");
    SegCut::Image2D imageOut = image;

    std::string outputFolder = "../output/refundFlow/square/square-20";

    for(int i=0;i<200;++i)
    {
        computeFlow(image,
                    gluedCurveLength,
                    imageOut,
                    outputFolder,
                    i);


        drawCurvatureMaps(image,
                          gluedCurveLength,
                          outputFolder,
                          std::to_string(i));


        saveImage(imageOut,
                  outputFolder,
                  std::to_string(i));

        image=imageOut;

        std::cout << "OK " << i << std::endl;
    }

    return 0;
}