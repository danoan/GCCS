
#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "utils.h"
#include "issueUtils.h"
#include "Artist.h"
#include "FlowGraphBuilder.h"

using namespace UtilsTypes;


void getCutFilter(ListDigraph& fg,
                  Preflow<ListDigraph,ListDigraph::ArcMap<double> >& flow,
                  ListDigraph::NodeMap<bool>& nodesInTheCut,
                  ListDigraph::ArcMap<bool>& arcsInTheCut)
{

    for(ListDigraph::ArcIt a(fg);a!=INVALID;++a){
        ListDigraph::Node s =fg.source(a);
        ListDigraph::Node t =fg.target(a);

        if( (nodesInTheCut[s] && !nodesInTheCut[t])  )
        {
            arcsInTheCut[a] = true;
        }
    }

}

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
                      FlowGraphBuilder** fgb,
                      std::map<Z2i::SCell,double>& weightMap) //it should be an empty pointer
{
    Curve intCurvePriorGS,extCurvePriorGS;
    KSpace KImage;

    SegCut::Image2D dilatedImage(mask.domain());

    dilate(dilatedImage, mask, 1);
    computeBoundaryCurve(intCurvePriorGS, KImage, mask, 100);
    computeBoundaryCurve(extCurvePriorGS, KImage, dilatedImage);


    setGridCurveWeight(intCurvePriorGS,
                       KImage,
                       weightMap);

    setGridCurveWeight(extCurvePriorGS,
                       KImage,
                       weightMap);


    ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurvePriorGS,extCurvePriorGS);
    SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);

    std::cout << "INT CURVE" <<std::endl;
    for(auto it=intCurvePriorGS.begin();it!=intCurvePriorGS.end();++it){
        std::cout << *it << "::" << weightMap[*it] << std::endl;
    }
    std::cout << "EXT CURVE" <<std::endl;
    for(auto it=extCurvePriorGS.begin();it!=extCurvePriorGS.end();++it){
        std::cout << *it << "::" << weightMap[*it] << std::endl;
    }


    setGluedCurveWeight(gcsRange,KImage,gluedCurveLength,weightMap);

//    Issue::updateGluedWeightUsingFP(seedRange,KImage,weightMap);


    std::vector< Curve > curvesVector = { intCurvePriorGS,extCurvePriorGS };

    *fgb = new FlowGraphBuilder(curvesVector,
                                KImage,
                                gluedCurveLength);

    (*fgb)->addPair(0,1);
}

void drawCurvatureMaps(Image2D& image,
                       std::map<Z2i::SCell,double>& weightMap,
                       std::string outputFolder,
                       int iteration)
{
    boost::filesystem::path p1(outputFolder.c_str());
    boost::filesystem::create_directories(p1);

    Board2D board;
    Curve intCurve,extCurve;
    KSpace KImage;

    SegCut::Image2D dilatedImage(image.domain());

    dilate(dilatedImage,image,1);
    computeBoundaryCurve(intCurve,KImage,image,100);
    computeBoundaryCurve(extCurve,KImage,dilatedImage);

//    erode(dilatedImage,image,1);
//    computeBoundaryCurve(extCurve,KImage,image,100);
//    computeBoundaryCurve(intCurve,KImage,dilatedImage);


    std::vector<Z2i::SCell> intConnection;
    std::vector<Z2i::SCell> extConnection;
    std::vector<Z2i::SCell> makeConvexConnection;


    ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);
    SeedToGluedCurveRangeFunctor stgcF(10);
    GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);

    for(GluedCurveIteratorPair it=gcsRange.begin();it!=gcsRange.end();++it){
        ConnectorType ct = it->first.connectorType();

        switch(ct) {

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

    double cmin=100;
    double cmax=-100;
    for(int i=0;i<2;i++) {

        drawCurvatureMap(intConnection.begin(),
                         intConnection.end(),
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

    std::string intConnsOutputFilePath = outputFolder +
                                         "/intconnCurvatureMap" + std::to_string(iteration) + ".eps";
    board.save(intConnsOutputFilePath.c_str());

    board.clear();

    cmin=100;
    cmax=-100;
    for(int i=0;i<2;i++) {
        drawCurvatureMap(extConnection.begin(),
                         extConnection.end(),
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

    std::string extConnsOutputFilePath = outputFolder +
                                         "/extconnCurvatureMap" + std::to_string(iteration) + ".eps";

    board.save(extConnsOutputFilePath.c_str());



    board.clear();

    cmin=100;
    cmax=-100;
    if (makeConvexConnection.size() > 0) {
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

typedef SubDigraph< ListDigraph,ListDigraph::NodeMap<bool> > MySubGraph;
double drawCutUpdateImage(FlowGraphBuilder* fgb,
                          std::map<Z2i::SCell,double>& weightMap,
                          MySubGraph** sg, //must be an empty pointer
                          std::string& cutOutputPath,
                          Image2D& out,
                          std::string& imageOutputPath) {
    fgb->operator()(weightMap);
    fgb->draw();


    Preflow<ListDigraph, ListDigraph::ArcMap<double> > flow = fgb->preparePreFlow();

    flow.run();
    double v = flow.flowValue();
    std::cout << "flow value::" << v << std::endl;

    ListDigraph::NodeMap<bool> node_filter(fgb->graph(), false);
    ListDigraph::ArcMap<bool> arc_filter(fgb->graph(), true);
    for (ListDigraph::NodeIt n(fgb->graph()); n != INVALID; ++n) {
        if (flow.minCut(n)) {
            node_filter[n] = true;
        } else {
            node_filter[n] = false;
        }
    }


    KSpace KImage;
    Image2D tempMask = GenericReader<Image2D>::import("../images/flow-evolution/out6-mask.pgm");

    KImage.init(tempMask.domain().lowerBound(), tempMask.domain().upperBound(), true);
    Curve tempCurveMask;

    {
        KSpace::SCell pixelModel = KImage.sCell(DGtal::Z2i::RealVector(1, 1), true);
        KSpace::SCell imagePixel = KImage.sCell(DGtal::Z2i::RealVector(241, 265), pixelModel);
        KSpace::SCell imageBel = KImage.sLowerIncident(imagePixel)[0];

        computeBoundaryCurve(tempCurveMask, KImage, tempMask, 0, imageBel);
    }


    UnsignedSCellComparison uscmp;
    std::set<KSpace::SCell,UnsignedSCellComparison> myComp(uscmp);
    std::set<KSpace::SCell> pixelsInTheMask;
    for(auto it =tempCurveMask.begin();it!=tempCurveMask.end();++it){
        KSpace::SCell s = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
        KSpace::SCell t = KImage.sIndirectIncident(*it,KImage.sOrthDir(*it));
        pixelsInTheMask.insert(s);
        pixelsInTheMask.insert(t);
    }


    ListDigraph::ArcMap<bool> arcsInTheCut(fgb->graph(),false);
    getCutFilter(fgb->graph(),flow,node_filter,arcsInTheCut);

    Palette palette;
    ListDigraph::ArcMap<int> arcColors(fgb->graph());
    double stemp = 0;
    for(ListDigraph::ArcIt a(fgb->graph());a!=INVALID;++a){
        if(arcsInTheCut[a]){
            arcColors[a] = 2;
            ListDigraph::Node s,t;
            s = fgb->graph().source(a);
            t = fgb->graph().target(a);

            KSpace::SCell scellS = fgb->pixelsMap()[s];
            KSpace::SCell scellT = fgb->pixelsMap()[t];

            if( pixelsInTheMask.find(scellS)!= pixelsInTheMask.end() ||
                pixelsInTheMask.find(scellT)!=pixelsInTheMask.end() )
            {
                arcColors[a] = 1;
                std::cout << fgb->coordsMap()[s] << "," << fgb->coordsMap()[t] << std::endl;
            }

            stemp += fgb->getEdgeWeight()[a];
        }
    }


    std::vector< std::pair< KSpace::Vector, KSpace::Vector> > toRemove =
            { std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(491,539), KSpace::Vector(489,539) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(491,539), KSpace::Vector(491,541) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(491,539), KSpace::Vector(493,539) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(495,539), KSpace::Vector(493,539) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(495,539), KSpace::Vector(495,541) ),
            };


    std::vector< std::pair< KSpace::Vector, KSpace::Vector> > toAdd =
            { std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(491,535), KSpace::Vector(491,537) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(493,537), KSpace::Vector(491,537) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(495,537), KSpace::Vector(495,539) ),
              std::pair< KSpace::Vector,KSpace::Vector >( KSpace::Vector(497,539), KSpace::Vector(495,539) ),
            };



    std::cout << "Sum edges on the cut filter:: " << stemp << std::endl;


    double sumToRemove=0;
    double sumToAdd=0;
    for(ListDigraph::ArcIt a(fgb->graph());a!=INVALID;++a){
        ListDigraph::Node s,t;
        s = fgb->graph().source(a);
        t = fgb->graph().target(a);

        KSpace::SCell scellS = fgb->pixelsMap()[s];
        KSpace::SCell scellT = fgb->pixelsMap()[t];

        for(auto it=toRemove.begin();it!=toRemove.end();it++){
            if( KImage.sKCoords(scellS) == it->first &&
                KImage.sKCoords(scellT) == it->second )
            {
                std::cout << "Remove Edge " << it->first << it->second << " weight::" << fgb->getEdgeWeight()[a] << std::endl;
                sumToRemove += fgb->getEdgeWeight()[a];
                break;
            }
        }

        for(auto it=toAdd.begin();it!=toAdd.end();it++){
            if( KImage.sKCoords(scellS) == it->first &&
                KImage.sKCoords(scellT) == it->second )
            {
                sumToAdd += fgb->getEdgeWeight()[a];
                break;
            }
        }

    }

    std::cout << "Sum to Remove::" << sumToRemove << std::endl;
    std::cout << "Sum to Add::" << sumToAdd << std::endl;

    graphToEps(fgb->graph(),"cutEdges.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .arcColors(composeMap(palette,arcColors))
            .coords(fgb->coordsMap())
            .run();


    *sg = new MySubGraph(fgb->graph(),node_filter,arc_filter);

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



    return v;
}

void particularCurve()
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

    board << c;
    board.saveEPS("particularCurve.eps");

    std::vector<double> estimations;
    curvatureEstimatorsGridCurve(c.begin(),c.end(),KImage,estimations,false);
    updateToSquared(estimations.begin(),estimations.end());

    double cmin=100;
    double cmax=-100;
    max_and_min(estimations,cmin,cmax);
    draw(estimations,c.begin(),c.end(),board,cmin,cmax);

    board.saveEPS("particularCurve.eps");

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

    std::string imgPath = "../images/flow-evolution/out20.pgm";
    SegCut::Image2D image = GenericReader<SegCut::Image2D>::import(imgPath);

    std::string outImageFolder = "issueout/out20";
    std::string cutOutputPath = "issueout/out20/cut.eps";
    std::string imageOutputPath = "issueout/out20/image.pgm";

    FlowGraphBuilder* fgb;
    MySubGraph* sg;

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
    EA.drawAllGluedCurves(imgPath,outImageFolder+"/gluedCurves");

    Issue::drawFaithfulPolygon(image,outImageFolder);

    std::map<Z2i::SCell,double> weightMap;
    prepareFlowGraph(image,
                     gluedCurveLength,
                     &fgb,
                     weightMap);

    drawCurvatureMaps(image,weightMap,outImageFolder,1);
    particularCurve();

    drawCutUpdateImage(fgb,
                       weightMap,
                       &sg,
                       cutOutputPath,
                       image,
                       imageOutputPath
    );

    return 0;
}