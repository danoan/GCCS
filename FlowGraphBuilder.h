#ifndef SEGBYCUT_FLOWGRAPHBUILDER_H
#define SEGBYCUT_FLOWGRAPHBUILDER_H

#include <ctime>
#include <cstdlib>

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;

#include "typeDefs.h"
#include "utils.h"

class FlowGraphBuilder{

public:
    typedef dim2::Point<int> LemonPoint;

    FlowGraphBuilder( Curve intCurve,
                      Curve extCurve,
                      KSpace KImage,
                      unsigned int gluedCurveLength,
                      double fidelityWeight=0.0,
                      double curvatureWeight=1.0):intCurve(intCurve),
                                                  extCurve(extCurve),
                                                  gluedCurveLength(gluedCurveLength),
                                                  KImage(KImage),
                                                  fFactor(fidelityWeight),wFactor(curvatureWeight),
                                                  edgeWeight(fg),
                                                  pixelMap(fg),
                                                  coords(fg)
    {
        srand(time(NULL));
        sourceNode = fg.addNode();
        targetNode = fg.addNode();

    };

    void operator()(std::map<Z2i::SCell,double>& weightMap);
    
    void draw();
    
    ListDigraph& graph(){return fg;};
    ListDigraph::Node& source(){return sourceNode;}
    ListDigraph::Node& target(){return targetNode;}

    ListDigraph::NodeMap< LemonPoint >& coordsMap(){return coords;}
    ListDigraph::NodeMap<Z2i::SCell>& pixelsMap(){return pixelMap;}
    
    
    Preflow <ListDigraph,ListDigraph::ArcMap<double> > preparePreFlow(){
        Preflow <ListDigraph,ListDigraph::ArcMap<double> > flow(fg,edgeWeight,sourceNode,targetNode);
        return flow;
    }

private:
    void randomizeCurve(std::vector<Z2i::SCell>::const_iterator cBegin,
                        std::vector<Z2i::SCell>::const_iterator cEnd,
                        int size,
                        std::vector<Z2i::SCell>& scellsRand);

    void createEdgeFromLinel(Curve::SCell& linel,
                             KSpace& KImage,
                             std::map<Z2i::SCell,double>& weightMap,
                             ListDigraph::Node& sourcePixelNode,
                             ListDigraph::Node& targetPixelNode);

    void createGridCurveEdges(Curve::ConstIterator curveBegin,
                              Curve::ConstIterator curveEnd,
                              ListDigraph::Node terminalToConnect,
                              bool toTerminal,
                              std::map<Z2i::SCell,double>& weightMap);

    void createGluedCurveEdges(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                               SegCut::GluedCurveIteratorPair gluedRangeEnd,
                               std::map<Z2i::SCell,double>& weightMap);


private:
    int xavg=0;
    int yavg=0;

    unsigned int gluedCurveLength;

    KSpace KImage;

    Curve intCurve;
    Curve extCurve;

    double fFactor;
    double wFactor;

    ListDigraph fg;

    std::map<KSpace::Point,ListDigraph::Node> coordToNode;
    ListDigraph::ArcMap<double> edgeWeight;
    ListDigraph::NodeMap<Z2i::SCell> pixelMap;
    ListDigraph::NodeMap< LemonPoint > coords;

    ListDigraph::Node sourceNode;
    ListDigraph::Node targetNode;
};

#endif //SEGBYCUT_FLOWGRAPHBUILDER_H
