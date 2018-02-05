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

namespace Stats{
    extern int gluedEdges;
    extern int sourceEdges;
    extern int targetEdges;
    extern int externalEdges;
    extern int internalEdges;

    extern int makeConvexEdges;
    extern int escapeEdges;

    extern double extMakeConvexCut;

    typedef SubDigraph< ListDigraph,ListDigraph::NodeMap<bool> > MySubGraph;
}

class UnsignedSCellComparison{
public:
    bool operator()(const SCell& s1, const SCell& s2) const{
        return s1.preCell().coordinates < s2.preCell().coordinates;
    }
};

class FlowGraphBuilder{

public:
    typedef dim2::Point<int> LemonPoint;

    FlowGraphBuilder( std::vector< Curve > curvesVector,
                      KSpace KImage,
                      unsigned int gluedCurveLength,
                      double fidelityWeight=0.0,
                      double curvatureWeight=1.0):curvesVector(curvesVector),
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

    void addPair(int i1,int i2){ curvesPairs.push_back( std::pair<int,int>(i1,i2) ); }

    void operator()(std::map<Z2i::SCell,double>& weightMap);
    void operator()(std::map<Z2i::SCell,double>& weightMap, bool multiplePair);
    
    void draw();

    Stats::MySubGraph curveSubgraph(Curve& curve,
                                    ListDigraph::NodeMap<bool>& node_filter,
                                    ListDigraph::ArcMap<bool>& arc_filter );

    
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
                             std::map<Z2i::SCell,double>& weightMap);

    int createGridCurveEdges(Curve::ConstIterator curveBegin,
                              Curve::ConstIterator curveEnd,
                              std::map<Z2i::SCell,double>& weightMap);

    void createGluedCurveEdges(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                               SegCut::GluedCurveIteratorPair gluedRangeEnd,
                               std::map<Z2i::SCell,double>& weightMap,
                               std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes);

    int createEscapeEdges(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                          Curve& fromCurve,
                          Curve& toCurve);


    int createSourceEdges(Curve& erodedCurve,
                          std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes);

    void connectNodeToPixel(ListDigraph::Node& sourceNode,
                            KSpace::SCell pTarget,
                            ListDigraph::Arc& a);

    void connectPixelToNode(KSpace::SCell pSource,
                            ListDigraph::Node& targetNode,
                            ListDigraph::Arc& a);

    void createNodeFromPixel(KSpace::SCell pixel,
                             ListDigraph::Node& node);

    void createEdgeFromPixels(KSpace::SCell pSource,
                              KSpace::SCell pTarget,
                              ListDigraph::Arc& a);

    int explore(KSpace::SCell candidate,
                 std::set<KSpace::SCell,UnsignedSCellComparison>& borderPoints,
                 std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes);


    int createTargetEdgesFromGluedSegments(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                           SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                           SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                           std::map<Z2i::SCell,double>& weightMap);

    int createTargetEdgesFromExteriorCurve(Curve& extCurve,
                                           std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                           std::map<Z2i::SCell,double>& weightMap);

    template<typename TType>
    int createFromIteratorsQueue(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                 std::queue<TType> intervals,std::map<Z2i::SCell,
                                 double>& weightMap)
    {

        int c=0;
        TType start,end;
        KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
        while(!intervals.empty()){
            start = intervals.front();
            intervals.pop();
            end = intervals.front();
            intervals.pop();

            TType itC = start;
            do{
                Stats::extMakeConvexCut+= weightMap[*itC];
                KSpace::SCell indirectIncidentPixel = KImage.sIndirectIncident(*itC,KImage.sOrthDir(*itC));
                KSpace::SCell directIncidentPixel = KImage.sDirectIncident(*itC,KImage.sOrthDir(*itC));

                visitedNodes.insert(directIncidentPixel);

                ListDigraph::Arc a;
                connectPixelToNode(indirectIncidentPixel,targetNode,a);
                edgeWeight[a] = 10;

                ++c;

                if(itC==end) break;
                ++itC;
            }while(true);

        }

        return c;

    }


private:
    int xavg=0;
    int yavg=0;

    unsigned int gluedCurveLength;

    KSpace KImage;

    double fFactor;
    double wFactor;

    ListDigraph fg;

    std::vector<Curve> curvesVector;
    std::vector< std::pair<int,int> > curvesPairs;

    std::map<KSpace::Point,ListDigraph::Node> coordToNode;
    ListDigraph::ArcMap<double> edgeWeight;
    ListDigraph::NodeMap<Z2i::SCell> pixelMap;
    ListDigraph::NodeMap< LemonPoint > coords;

    ListDigraph::Node sourceNode;
    ListDigraph::Node targetNode;
};

#endif //SEGBYCUT_FLOWGRAPHBUILDER_H
