#ifndef FLOW_GRAPH_H
#define FLOW_GRAPH_H

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

#include "DGtal/helpers/StdDefs.h"

#include "ImageFlowData.h"

class FlowGraph
{
public:
    friend class FlowGraphBuilder;

    enum ArcType{
        SourceArc,TargetArc,IntExtGluedArc,ExtIntGluedArc,
        InternalCurveArc,ExternalCurveArc,EscapeArc,
        MakeConvexArc,RefundArc,ConflictArc
    };

    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef lemon::ListDigraph ListDigraph;

    typedef std::map<KSpace::Point,ListDigraph::Node> CoordToNodeMap;
    typedef std::map<SCell,ListDigraph::Arc> SCellArcMap;

    typedef lemon::dim2::Point<int> LemonPoint;

    typedef ListDigraph::NodeMap<SCell> NodeToPixelMap;    //Image reconstruction from cut
    typedef ListDigraph::NodeMap< LemonPoint > NodeToCoordMap;  //Graph Drawing

    typedef ListDigraph::ArcMap<double> ArcWeightMap;   //MinCut computation
    typedef ListDigraph::ArcMap<ArcType> ArcTypeMap;

    typedef DGtal::Z2i::Curve::ConstIterator SCellIteratorType;
    typedef DGtal::Circulator<SCellIteratorType> SCellCirculator;
    typedef std::pair<SCellCirculator,SCellCirculator> CirculatorPair;

    typedef ListDigraph::ArcMap<SCell> ArcSCellMap;
    typedef ListDigraph::ArcMap<CirculatorPair> ArcCirculatorMap; //DetourArc Computation

    typedef lemon::Preflow <ListDigraph,ListDigraph::ArcMap<double> > FlowComputer;


public:
    FlowGraph():pixelMap(digraph),
                arcWeightMap(digraph),
                arcTypeMap(digraph),
                arcCirculatorMap(digraph),
                arcSCellMap(digraph),
                coords(digraph)
    {
        sourceNode = digraph.addNode();
        targetNode = digraph.addNode();
    };


    int id(const ListDigraph::Arc& a){return digraph.id(a);}
    int id(const ListDigraph::Node& n){return digraph.id(n);}

    const ListDigraph::Node source(const ListDigraph::Arc& a){ return digraph.source(a); }
    const ListDigraph::Node target(const ListDigraph::Arc& a){ return digraph.target(a); }

    SCell pixel(const ListDigraph::Node& node){return pixelMap[node];}

    ListDigraph::Arc arc(const SCell& scell){ return scellArcMap[scell]; }
    ArcType arcType(const ListDigraph::Arc& a){ return arcTypeMap[a]; }

    CirculatorPair circulatorPair(const ListDigraph::Arc& a){ return arcCirculatorMap[a]; }

    double weight(const ListDigraph::Arc& a){ return arcWeightMap[a]; }

    int getGluedCurveLength(){ return gluedCurveLength; }
    int getConsecutiveGluedPairDistance(){ return consecutiveGluedPairDistance; }
    int getDiffDistance(){ return diffDistance; }
    ImageFlowData::FlowMode getFlowMode(){ return flowGraphMode; }

    ListDigraph::Node& source(){return sourceNode;}
    ListDigraph::Node& target(){return targetNode;}

    NodeToCoordMap& coordsMap(){ return coords; };

    ListDigraph& graph(){return digraph;}

    FlowComputer prepareFlow(){ return FlowComputer(digraph,arcWeightMap,sourceNode,targetNode); }


private:
    int gluedCurveLength;
    int consecutiveGluedPairDistance;
    int diffDistance;
    ImageFlowData::FlowMode flowGraphMode;

    ListDigraph digraph;

    ListDigraph::Node sourceNode;
    ListDigraph::Node targetNode;

    CoordToNodeMap coordToNode;
    NodeToCoordMap coords;

    NodeToPixelMap pixelMap;

    ArcWeightMap arcWeightMap;
    ArcTypeMap arcTypeMap;
    ArcCirculatorMap arcCirculatorMap; //Valid entries only for Glued Arcs
    ArcSCellMap arcSCellMap;   //Valid entries only for Curve Arcs

    SCellArcMap scellArcMap;
};

#endif