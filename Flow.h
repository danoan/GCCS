#ifndef SEGBYCUT_FLOW_H
#define SEGBYCUT_FLOW_H

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "ImageFlowData.h"
#include "FlowGraphBuilder.h"

#include "weightSettings.h"
#include "FlowGraphQuery.h"

class Flow
{
public:
    friend class FlowGraphDebug;

    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

    typedef std::pair< Point, Point > MyKey;

public:
    static void matchFlows(ListDigraph::ArcMap<int>& mapping,
                           Flow& f1,
                           Flow& f2);

public:
    Flow(ImageFlowData& imageFlowData);

    Domain domain();
    ListDigraph& graph(){ return fgb.graph(); };
    std::map<SCell,double>& getWeightMap(){ return weightMap; };


    void initImageFlowDataAsInFlow(ImageFlowData& imageFlowData);
    double weight(ListDigraph::Arc& arc){ return fgb.getEdgeWeight()[arc]; }

    void detourArcsFilter(ListDigraph::ArcMap<bool> & detourArcs);
    std::map< FlowGraphQuery::ArcPair, std::vector< ListDigraph::Arc > >::const_iterator detourArcsMapBegin(){ return fgq.detourArcsBegin(); };
    std::map< FlowGraphQuery::ArcPair, std::vector< ListDigraph::Arc > >::const_iterator detourArcsMapEnd(){ return fgq.detourArcsEnd(); };

    void addRefundArcs(ListDigraph::Arc& intExtArc,
                       ListDigraph::Arc& extIntArc,
                       double weight);

    Image2D& baseImage(){ return imageFlowData.getOriginalImage(); }
    Image2D& dilatedImage(){ return imageFlowData.getDilatedImage(); }

    double cutValue();

    void updateImage(Image2D& out);
    bool hasChanges(Image2D& im1, Image2D& im2);
    FlowGraphBuilder& graphBuilder(){return fgb;}


private:
    void setArcsWeight();
    void getPixelsFilter(ListDigraph::NodeMap<bool>& pixelsFilter);
    void pixels(std::vector<SCell>& pixelsVector,
                ListDigraph::NodeMap<bool>& nodeFilter);
    void fillHoles(Image2D& out);

    static MyKey arcKey(FlowGraphBuilder& fgb,
                        ListDigraph::Arc& arc);

private:
    ImageFlowData imageFlowData;
    FlowGraphBuilder fgb;
    FlowGraphQuery fgq;

    std::map<SCell,double> weightMap;
};


#endif //SEGBYCUT_FLOW_H
