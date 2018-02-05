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
#include "FlowGraphBuilder.h"

using namespace UtilsTypes;

void tangentWeight(Curve::ConstIterator begin,
                   Curve::ConstIterator end,
                   KSpace& KImage,
                   std::vector< TangentVector >& estimationsTangent,
                   std::vector< double >& tangentWeightVector)
{
    auto it = begin;
    int i =0;
    do{
        UtilsTypes::KSpace::Point pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
        UtilsTypes::KSpace::Point pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );

        UtilsTypes::KSpace::Point scellVector = pTarget-pSource;

        tangentWeightVector.push_back( fabs( estimationsTangent[i].dot(scellVector) ) );
        ++it;
        ++i;
    }while(it!=end);
}

void setGridCurveWeight(Curve curvePriorGS,
                        KSpace& KImage,
                        std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> curvatureEstimations;
    curvatureEstimatorsGridCurve(curvePriorGS.begin(),
                                 curvePriorGS.end(),
                                 KImage,
                                 curvatureEstimations);


    updateToSquared(curvatureEstimations.begin(),curvatureEstimations.end());


    std::vector<TangentVector> tangentEstimations;
    tangentEstimatorsGridCurve(curvePriorGS.begin(),
                               curvePriorGS.end(),
                               KImage,
                               tangentEstimations);


    std::vector<double> tangentWeightVector;
    tangentWeight(curvePriorGS.begin(),
                  curvePriorGS.end(),
                  KImage,
                  tangentEstimations,
                  tangentWeightVector);


    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] = curvatureEstimations[i];
            ++i;
        }
    }

    {
        int i =0;
        for(auto it=curvePriorGS.begin();it!=curvePriorGS.end();++it){
            weightMap[*it] *=tangentWeightVector[i];
            ++i;
        }
    }

}

void setGluedCurveWeight(GluedCurveSetRange gcRange,
                         KSpace& KImage,
                         unsigned int gluedCurveLength,
                         std::map<Z2i::SCell,double>& weightMap)
{
    std::vector<double> estimationsCurvature;
    curvatureEstimatorsConnections(gcRange.begin(),gcRange.end(),KImage,gluedCurveLength,estimationsCurvature);

    updateToSquared(estimationsCurvature.begin(),estimationsCurvature.end());

    {
        int i = 0;
        for (GluedCurveIteratorPair it = gcRange.begin(); it != gcRange.end(); ++it) {

            auto itC = it->first.connectorsBegin();
            do{
                weightMap[*itC] = estimationsCurvature[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);

        }
    }

    std::vector<TangentVector> estimationsTangent;
    tangentEstimatorsConnections(gcRange.begin(),gcRange.end(),KImage,gluedCurveLength,estimationsTangent);


    std::vector<double> tangentWeightVector;
    {
        GluedCurveIteratorPair it = gcRange.begin();
        int i = 0;
        for (GluedCurveIteratorPair it = gcRange.begin(); it != gcRange.end(); ++it) {

            auto itC = it->first.connectorsBegin();
            do {

                KSpace::SCell linel = *itC;

                UtilsTypes::KSpace::Point pTarget = KImage.sCoords(KImage.sDirectIncident(linel, *KImage.sDirs(linel)));
                UtilsTypes::KSpace::Point pSource = KImage.sCoords(KImage.sIndirectIncident(linel, *KImage.sDirs(linel)));

                UtilsTypes::KSpace::Point scellVector = pTarget - pSource;

                tangentWeightVector.push_back(fabs(estimationsTangent[i].dot(scellVector)));

                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);

        }

    }

    {
        int i = 0;
        for (GluedCurveIteratorPair it = gcRange.begin(); it != gcRange.end(); ++it) {
            auto itC = it->first.connectorsBegin();
            do {
                weightMap[*itC]*= tangentWeightVector[i];
                ++i;
                if(itC==it->first.connectorsEnd()) break;
                ++itC;
            }while(true);
        }
    }

}


void prepareFlowGraph(SegCut::Image2D& mask,
                      unsigned int gluedCurveLength,
                      std::map<Z2i::SCell,double>& weightMap,
                      FlowGraphBuilder** fgb,
                      bool toDilate) //it should be an empty pointer
{
    Curve intCurvePriorGS,extCurvePriorGS;
    KSpace KImage;

    SegCut::Image2D dilatedImage(mask.domain());

    if(toDilate){
        dilate(dilatedImage,mask,1);
        computeBoundaryCurve(intCurvePriorGS,KImage,mask,100);
        computeBoundaryCurve(extCurvePriorGS,KImage,dilatedImage);

//        gluedCurveLength = intCurvePriorGS.size()/2.0;
    }else{
        erode(dilatedImage,mask,1);
        computeBoundaryCurve(extCurvePriorGS,KImage,mask,100);
        computeBoundaryCurve(intCurvePriorGS,KImage,dilatedImage);

//        gluedCurveLength = extCurvePriorGS.size()/2.0;
    }

    Board2D board;
    board << intCurvePriorGS;
    board << extCurvePriorGS;
    board.saveEPS("int.eps");

    std::vector< Curve > curvesVector = { intCurvePriorGS,extCurvePriorGS};

    for(auto it=curvesVector.begin();it!=curvesVector.end();++it){
        setGridCurveWeight(*it,
                           KImage,
                           weightMap
        );
    }


    SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);


    ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurvePriorGS,extCurvePriorGS);
    GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);
    setGluedCurveWeight(gcsRange,KImage,gluedCurveLength,weightMap);


    *fgb = new FlowGraphBuilder(curvesVector,
                                KImage,
                                gluedCurveLength);

    (*fgb)->addPair(0,1);
}


typedef SubDigraph< ListDigraph,ListDigraph::NodeMap<bool> > MySubGraph;
double drawCutUpdateImage(FlowGraphBuilder* fgb,
                          std::map<Z2i::SCell,double>& weightMap,
                          MySubGraph** sg, //must be an empty pointer
                          std::string& cutOutputPath,
                          Image2D& out,
                          std::string& imageOutputPath)
{
    fgb->operator()(weightMap);
    fgb->draw();



    Preflow<ListDigraph,ListDigraph::ArcMap<double> > flow = fgb->preparePreFlow();

    flow.run();
    double v = flow.flowValue();
    std::cout << v << std::endl;

    ListDigraph::NodeMap<bool> node_filter(fgb->graph(),false);
    ListDigraph::ArcMap<bool> arc_filter(fgb->graph(),true);
    for(ListDigraph::NodeIt n(fgb->graph());n!=INVALID;++n){
        if( flow.minCut(n) ){
            node_filter[n] = true;
        }else{
            node_filter[n] = false;
        }
    }
    node_filter[fgb->source()]=false;



    *sg = new MySubGraph(fgb->graph(),node_filter,arc_filter);

    Palette palette;

    boost::filesystem::path p1(cutOutputPath.c_str());
    p1.remove_filename();
    boost::filesystem::create_directories(p1);

    graphToEps(**sg,cutOutputPath.c_str())
        .coords(fgb->coordsMap())
        .nodeScale(.0005)
        .arcWidthScale(0.0005)
        .drawArrows()
        .arrowWidth(0.25)
        .run();


    KSpace KImage;
    KImage.init(out.domain().lowerBound(),out.domain().upperBound(),true);
    for(MySubGraph::NodeIt n(**sg);n!=INVALID;++n){
        Z2i::SCell pixel = fgb->pixelsMap()[n];
        Z2i::Point p = KImage.sCoords(pixel);
        out.setValue(p,255);
    }


    /*Fill the holes*/
    for(ListDigraph::NodeIt n(fgb->graph());n!=INVALID;++n){
        Z2i::SCell pixel = fgb->pixelsMap()[n];
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

    
    return v;
}

void drawCurvatureMaps(Image2D& image,
                       std::map<Z2i::SCell,double>& weightMap,
                       std::string outputFolder,
                       int iteration) {
    Board2D board;
    Curve intCurve, extCurve;
    KSpace KImage;

    SegCut::Image2D dilatedImage(image.domain());

    dilate(dilatedImage, image, 1);
    computeBoundaryCurve(intCurve, KImage, image, 100);
    computeBoundaryCurve(extCurve, KImage, dilatedImage);

    std::vector<Z2i::SCell> intConnection;
    std::vector<Z2i::SCell> extConnection;
    std::vector<Z2i::SCell> makeConvexConnection;


    ConnectorSeedRangeType seedRange = getSeedRange(KImage, intCurve, extCurve);
    SeedToGluedCurveRangeFunctor stgcF(10);
    GluedCurveSetRange gcsRange(seedRange.begin(),
                                seedRange.end(),
                                stgcF);

    for (GluedCurveIteratorPair it = gcsRange.begin(); it != gcsRange.end(); ++it) {
        ConnectorType ct = it->first.connectorType();

        switch (ct) {

            case internToExtern:
                intConnection.push_back(it->first.linkSurfel());
                break;
            case externToIntern:
                extConnection.push_back(it->first.linkSurfel());
                break;
            case makeConvex:
                auto itC = it->first.connectorsBegin();
                do{
                    makeConvexConnection.push_back(*itC);
                    if(itC==it->first.connectorsEnd()) break;
                    ++itC;
                }while(true);
        }

    }

    double cmin = 100;
    double cmax = -100;
    for (int i = 0; i < 2; i++) {
        drawCurvatureMap(intCurve.begin(),
                         intCurve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );

        drawCurvatureMap(extCurve.begin(),
                         extCurve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );

        drawCurvatureMap(intConnection.begin(),
                         intConnection.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );
    }

    std::string intConnsOutputFilePath = outputFolder +
                                         "/intconnCurvatureMap" + std::to_string(iteration) + ".eps";
    board.save(intConnsOutputFilePath.c_str());

    board.clear();

    cmin = 100;
    cmax = -100;
    for (int i = 0; i < 2; i++) {
        drawCurvatureMap(intCurve.begin(),
                         intCurve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );

        drawCurvatureMap(extCurve.begin(),
                         extCurve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );

        drawCurvatureMap(extConnection.begin(),
                         extConnection.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );
    }

    std::string extConnsOutputFilePath = outputFolder +
                                         "/extconnCurvatureMap" + std::to_string(iteration) + ".eps";

    board.save(extConnsOutputFilePath.c_str());


    board.clear();

    cmin = 100;
    cmax = -100;
    for (int i = 0; i < 2; i++) {
        drawCurvatureMap(intCurve.begin(),
                         intCurve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );

        drawCurvatureMap(extCurve.begin(),
                         extCurve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap
        );
    }

    std::string noConnsOutputFilePath = outputFolder +
                                         "/noConnCurvatureMap" + std::to_string(iteration) + ".eps";

    board.save(noConnsOutputFilePath.c_str());


    board.clear();

    cmin = 100;
    cmax = -100;
    if (makeConvexConnection.size() > 0){
        for (int i = 0; i < 2; i++) {
            drawCurvatureMap(makeConvexConnection.begin(),
                             makeConvexConnection.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap
            );

            drawCurvatureMap(intCurve.begin(),
                             intCurve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap
            );

            drawCurvatureMap(extCurve.begin(),
                             extCurve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap
            );

        }
    }

    std::string makeConvexConnsOutputFilePath = outputFolder +
                                                "/makeConvexConnCurvatureMap" + std::to_string(iteration) + ".eps";

    board.save(makeConvexConnsOutputFilePath.c_str());
}


namespace Patch{
    bool solveShift;
    bool cross_element;
};

namespace UtilsTypes
{
    std::function< double(double) > toDouble = [](double x){return x;};
};

int main(){
    Patch::solveShift = false;
    Patch::cross_element = false;

    unsigned int gluedCurveLength = 7;

    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import("../images/flow-evolution/single_square.pgm");

    std::string outImageFolder = "output/images/flow-evolution/square-7";
    std::string cutOutputPath;
    std::string imageOutputPath;

    FlowGraphBuilder* fgb;
    MySubGraph* sg;
    for(int i=0;i<200;++i)
    {
        std::map<Z2i::SCell,double> weightMap;
        prepareFlowGraph(image,
                         gluedCurveLength,
                         weightMap,
                         &fgb,
                         true);

        drawCurvatureMaps(image,
                          weightMap,
                          outImageFolder,
                          i);

        cutOutputPath = outImageFolder + "/cutGraph" + std::to_string(i) + ".eps";
        imageOutputPath = outImageFolder + "/out" + std::to_string(i+1) + ".pgm";


        if(false){
            Image2D temp(image.domain());
            erode(temp,image,1);
            image = temp;
        }

        drawCutUpdateImage(fgb,
                           weightMap,
                           &sg,
                           cutOutputPath,
                           image,
                           imageOutputPath
        );


//        if(i%30==29){
//            Image2D rImage(image.domain());
//            resize(image,rImage);
//            image = rImage;
//        }



        std::cout << "OK " << i << std::endl;

        delete fgb;
        delete sg;
    }

    return 0;
}