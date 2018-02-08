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
#include "ImageFlowData.h"

namespace Development{
      extern bool invertGluedArcs;  
};

class UnsignedSCellComparison{
public:
    bool operator()(const SCell& s1, const SCell& s2) const{
        return s1.preCell().coordinates < s2.preCell().coordinates;
    }
};

class FlowGraphBuilder{

public:
    friend class FlowGraphDebug;

    enum ArcType{
        SourceArc,TargetArc,GluedArc,CurveArc,EscapeArc,MakeConvexArc,
    };

    typedef dim2::Point<int> LemonPoint;
    typedef SubDigraph< ListDigraph,ListDigraph::NodeMap<bool> > SubGraph;
    typedef Preflow <ListDigraph,ListDigraph::ArcMap<double> > Flow;

    typedef std::set<KSpace::SCell,UnsignedSCellComparison> UnsignedSCellSet;

    typedef std::map<KSpace::Point,ListDigraph::Node> CoordToNodeMap;

    typedef ListDigraph::NodeMap<Z2i::SCell> NodeToPixelMap;
    typedef ListDigraph::NodeMap< LemonPoint > NodeToCoordMap;

    typedef ListDigraph::ArcMap<double> ArcWeightMap;
    typedef ListDigraph::ArcMap<ArcType> ArcTypeMap;

    typedef std::map<Z2i::SCell,double> LinelWeightMap;

    FlowGraphBuilder( ImageFlowData& imageFlowData):imageFlowData(imageFlowData),
                                                    pixelMap(fg),
                                                    coords(fg),
                                                    arcWeight(fg),
                                                    arcType(fg)
    {
        srand(time(NULL));
        sourceNode = fg.addNode();
        targetNode = fg.addNode();

    };

    void initializeStats();
    void operator()(std::map<Z2i::SCell,double>& weightMap);



    ListDigraph& graph(){return fg;};
    ListDigraph::Node& source(){return sourceNode;}
    ListDigraph::Node& target(){return targetNode;}

    ListDigraph::NodeMap< LemonPoint >& coordsMap(){return coords;}
    ListDigraph::NodeMap<Z2i::SCell>& pixelsMap(){return pixelMap;}
    ListDigraph::ArcMap<double>& getEdgeWeight(){return arcWeight;}
    
    
    Flow preparePreFlow(){
        Flow flow(fg,arcWeight,sourceNode,targetNode);
        return flow;
    }

private:

    void createArcFromLinel(Curve::SCell& linel,
                            std::map<Z2i::SCell,double>& weightMap,
                            ArcType at,
                            bool invert=false);

    void createCurveArcs(Curve::ConstIterator curveBegin,
                         Curve::ConstIterator curveEnd,
                         std::map<Z2i::SCell,double>& weightMap);

    void createGluedArcs(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                         SegCut::GluedCurveIteratorPair gluedRangeEnd,
                         std::map<Z2i::SCell,double>& weightMap,
                         UnsignedSCellSet& visitedNodes);

    void createEscapeArcs(Curve& fromCurve,
                          Curve& toCurve,
                          UnsignedSCellSet& visitedNodes);


    void createSourceArcs(Curve& erodedCurve,
                          std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes);

    void connectNodeToPixel(ListDigraph::Node& sourceNode,
                            KSpace::SCell pTarget,
                            ListDigraph::Arc& a);

    void connectPixelToNode(KSpace::SCell pSource,
                            ListDigraph::Node& targetNode,
                            ListDigraph::Arc& a);

    void createNodeFromPixel(KSpace::SCell pixel,
                             ListDigraph::Node& node);

    void createArcFromPixels(KSpace::SCell pSource,
                              KSpace::SCell pTarget,
                              ListDigraph::Arc& a);

    void explore(KSpace::SCell candidate,
                 UnsignedSCellSet& forbiddenPoints,
                 UnsignedSCellSet& visitedNodes);


    void createTargetArcsFromGluedSegments(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                           SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                           std::map<Z2i::SCell,double>& weightMap,
                                           UnsignedSCellSet& visitedNodes);

    void createTargetArcsFromExteriorCurve(Curve& extCurve,
                                           UnsignedSCellSet& visitedNodes);

    template<typename TType>
    int createFromIteratorsQueue(std::queue<TType> intervals,
                                 std::map<Z2i::SCell,double>& weightMap,
                                 UnsignedSCellSet& visitedNodes);



private:
    ImageFlowData imageFlowData;

    ListDigraph fg;

    ListDigraph::Node sourceNode;
    ListDigraph::Node targetNode;

    CoordToNodeMap coordToNode;

    NodeToPixelMap pixelMap;
    NodeToCoordMap coords;

    ArcWeightMap arcWeight;
    ArcTypeMap arcType;

};

#endif //SEGBYCUT_FLOWGRAPHBUILDER_H
