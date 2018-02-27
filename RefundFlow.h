#ifndef SEGBYCUT_FLOW_H
#define SEGBYCUT_FLOW_H

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "FlowGraph/ImageFlowData.h"
#include "FlowGraph/FlowGraphBuilder.h"
#include "FlowGraph/FlowGraph.h"
#include "FlowGraph/FlowGraphQuery.h"

#include "FlowGraph/weightSettings.h"


class RefundFlow
{
public:
    friend class FlowGraphDebug;

    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
    typedef FlowGraphBuilder::LinelWeightMap LinelWeightMap;

    enum IterationReturnType {Continue,Stop};

public:


public:
    RefundFlow(ImageFlowData& imageFlowData);


    void run();
    Image2D& outputImage(){ return imageOut; }




private:
    void addRefundArcs(FlowGraph& fg,
                       ListDigraph::Arc& intExtArc,
                       ListDigraph::Arc& extIntArc,
                       double weight);

    void addConflictArcs(FlowGraph& fg,
                         ListDigraph::Arc& intExtArc,
                         ListDigraph::Arc& extIntArc,
                         double weight);

    void computeCutEnergyDifference(FlowGraph& f1,
                                    Image2D& image,
                                    ListDigraph::ArcMap<bool>& arcFilter,
                                    ListDigraph::ArcMap<double>& arcDiff);

    int findConflicts(FlowGraph& fg,
                      std::set<FlowGraphQuery::ArcPair>& listOfPairs,
                      FlowGraphQuery::ArcPair& pairToTest,
                      FlowGraphQuery::ArcPair conflictedPair[]);

    int runIteration(FlowGraph& fg,
                     Image2D& partialImage);

    double cutValue(FlowGraph& fg);

    void updateImage(FlowGraph& fg,
                     Image2D& out);

    void fillHoles(FlowGraph& fg,
                   Image2D& out);

    void pixels(FlowGraph& fg,
                std::vector<SCell>& pixelsVector,
                ListDigraph::NodeMap<bool>& nodeFilter);


    static void matchFlows(ListDigraph::ArcMap<int>& mapping,
                           FlowGraph& f1,
                           FlowGraph& f2);


private:
    ImageFlowData& imageFlowData;
    Image2D imageOut;

};


#endif //SEGBYCUT_FLOW_H
