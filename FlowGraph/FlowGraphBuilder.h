#ifndef SEGBYCUT_FLOWGRAPHBUILDER_H
#define SEGBYCUT_FLOWGRAPHBUILDER_H

#include <ctime>
#include <cstdlib>

#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>

#include <lemon/preflow.h>
#include <lemon/adaptors.h>

using namespace lemon;

#include "ImageFlowData.h"
#include "FlowGraph.h"

namespace Development{
      extern bool invertGluedArcs;  
};

class UnsignedSCellComparison{
public:
    bool operator()(const SCell& s1, const SCell& s2) const{
        return s1.preCell().coordinates < s2.preCell().coordinates;
    }
};

class FlowGraphBuilder
{

public:
    typedef dim2::Point<int> LemonPoint;

    typedef std::set<KSpace::SCell,UnsignedSCellComparison> UnsignedSCellSet;
    typedef std::map<Z2i::SCell,double> LinelWeightMap;

    typedef Curve::ConstIterator SCellIteratorType;
    typedef DGtal::Circulator<SCellIteratorType> SCellCirculator;
    typedef std::pair<SCellCirculator,SCellCirculator> CirculatorPair;


public:
    FlowGraphBuilder(FlowGraph& fg,
                     ImageFlowData& imageFlowData,
                     LinelWeightMap& weightMap);

    static void addArc(FlowGraph& fg,
                       ListDigraph::Node& u,
                       ListDigraph::Node& v,
                       FlowGraph::ArcType at,
                       double weight);


private:

    void addArc(ListDigraph::Arc& arc,
                ListDigraph::Node& u,
                ListDigraph::Node& v,
                FlowGraph::ArcType at,
                double weight);

    void addArc(ListDigraph::Node& u,
                ListDigraph::Node& v,
                FlowGraph::ArcType at,
                double weight);

    void createArcFromLinel(ListDigraph::Arc& a,
                            Curve::SCell& linel,
                            LinelWeightMap& weightMap,
                            FlowGraph::ArcType at,
                            bool invert);

    void createArcFromLinel(Curve::SCell& linel,
                            std::map<Z2i::SCell,double>& weightMap,
                            FlowGraph::ArcType at,
                            bool invert);

    void createArcFromPixels(ListDigraph::Arc& arc,
                             KSpace::SCell pSource,
                             KSpace::SCell pTarget,
                             FlowGraph::ArcType at,
                             double weight);

    void createArcFromPixels(KSpace::SCell pSource,
                             KSpace::SCell pTarget,
                             FlowGraph::ArcType at,
                             double weight);

    void connectNodeToPixel(ListDigraph::Node& sourceNode,
                            KSpace::SCell pTarget,
                            ListDigraph::Arc& a);

    void connectPixelToNode(KSpace::SCell pSource,
                            ListDigraph::Node& targetNode,
                            ListDigraph::Arc& a);

    void createNodeFromPixel(KSpace::SCell pixel,
                             ListDigraph::Node& node);


    void createCurveArcs(Curve::ConstIterator curveBegin,
                         Curve::ConstIterator curveEnd,
                         ImageFlowData::CurveType ct,
                         LinelWeightMap& weightMap);

    void createGluedArcs(ImageFlowData::GluedCurveIteratorPair gluedRangeBegin,
                         ImageFlowData::GluedCurveIteratorPair gluedRangeEnd,
                         LinelWeightMap& weightMap,
                         UnsignedSCellSet& visitedNodes);

    void createEscapeArcs(Curve& fromCurve,
                          Curve& toCurve,
                          UnsignedSCellSet& visitedNodes);


    void createSourceArcs(Curve& erodedCurve,
                          std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes);



    void explore(KSpace::SCell candidate,
                 UnsignedSCellSet& forbiddenPoints,
                 UnsignedSCellSet& visitedNodes);


    void createTargetArcsFromGluedSegments(ImageFlowData::GluedCurveIteratorPair gluedRangeBegin,
                                           ImageFlowData::GluedCurveIteratorPair gluedRangeEnd,
                                           LinelWeightMap& weightMap,
                                           UnsignedSCellSet& visitedNodes);

    void createTargetArcsFromExteriorCurve(Curve& extCurve,
                                           UnsignedSCellSet& visitedNodes);

    template<typename TType>
    int createFromIteratorsQueue(std::queue<TType> intervals,
                                 std::map<Z2i::SCell,double>& weightMap,
                                 UnsignedSCellSet& visitedNodes);

    void setTerminalsCoordinates();

private:
    FlowGraph& fg;
    ImageFlowData& imageFlowData;
    LinelWeightMap& weightMap;
};



#endif //SEGBYCUT_FLOWGRAPHBUILDER_H
