#include "Flow.h"

Flow::MyKey Flow::arcKey(FlowGraphBuilder& fgb,
                   ListDigraph::Arc& arc)
{
    SCell source = fgb.pixelsMap()[ fgb.graph().source(arc) ];
    SCell target = fgb.pixelsMap()[ fgb.graph().target(arc) ];

    return MyKey( source.preCell().coordinates,
                  target.preCell().coordinates );
}

void Flow::matchFlows(ListDigraph::ArcMap<int>& mapping,
                      Flow &f1,
                      Flow &f2)
{
    std::map< MyKey,ListDigraph::Arc > unifiedArcMapping;
    for(ListDigraph::ArcIt arcF1(f1.graph());arcF1!=INVALID;++arcF1)
    {
        unifiedArcMapping[Flow::arcKey(f1.fgb, arcF1)] = arcF1;
    }

    for(ListDigraph::ArcIt arcF2(f2.graph());arcF2!=INVALID;++arcF2)
    {
        MyKey arcF2Key = Flow::arcKey(f2.fgb, arcF2);

        if (unifiedArcMapping.find(arcF2Key) != unifiedArcMapping.end()) {
            ListDigraph::Arc arcF1 = unifiedArcMapping[arcF2Key];
            mapping[arcF1] = f2.graph().id(arcF2);
            std::cout << f2.graph().id(arcF2) << std::endl;
        }
    }

}

Flow::Flow(ImageFlowData& imageFlowData):imageFlowData(imageFlowData),
                                         fgb(imageFlowData),
                                         fgq(fgb)
{
    setArcsWeight();
    fgb(weightMap);
}

Domain Flow::domain()
{
    KSpace& KImage = imageFlowData.getKSpace();
    return Domain(KImage.lowerBound(),KImage.upperBound());
}

void Flow::initImageFlowDataAsInFlow(ImageFlowData& newImageFlowData)
{
    newImageFlowData.init( imageFlowData.getFlowMode(),
                           imageFlowData.getGluedCurveLength()
    );
}

void Flow::detourArcsFilter(ListDigraph::ArcMap<bool>& detourArcs)
{
    for (FlowGraphQuery::DetourArcMapIterator dami = fgq.detourArcsBegin(); dami != fgq.detourArcsEnd(); ++dami) {
        const std::vector<ListDigraph::Arc> &values = dami->second;

        for (FlowGraphQuery::DetourArcIterator dai = values.begin(); dai != values.end(); ++dai) {
            detourArcs[*dai] =true;

            if( fgb.getArcType(*dai) != FlowGraphBuilder::ArcType::ExternalCurveArc) throw "Not ExternalCurveArc identified as a DetourArc.";
        }
    }
}

void Flow::addRefundArcs(ListDigraph::Arc& intExtArc,
                         ListDigraph::Arc& extIntArc,
                         double weight)
{
    ListDigraph::Node u1 = fgb.graph().source(intExtArc);
    ListDigraph::Node v1 = fgb.graph().target(extIntArc);

    ListDigraph::Node u2 = fgb.graph().source(extIntArc);
    ListDigraph::Node v2 = fgb.graph().target(intExtArc);


    fgb.addRefundArc(u1,v1,weight);
    fgb.addRefundArc(u2,v2,weight);
}

double Flow::cutValue()
{
    FlowGraphBuilder::FlowComputer flowComputer =  fgb.preparePreFlow();
    flowComputer.run();
    return flowComputer.flowValue();
}

void Flow::setArcsWeight()
{

    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        setGridCurveWeight(it->curve,
                           imageFlowData.getKSpace(),
                           weightMap
        );
    }

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        setGluedCurveWeight( it->gcsRangeBegin(),
                             it->gcsRangeEnd(),
                             imageFlowData.getKSpace(),
                             imageFlowData.getGluedCurveLength(),
                             weightMap);
    }
}


void Flow::updateImage(Image2D& updatedImage)
{
    ListDigraph::NodeMap<bool> nodeFilter(fgb.graph(),false);
    ListDigraph::ArcMap<bool> arcFilter(fgb.graph(),false);
    ListDigraph::ArcMap<bool> cutFilter(fgb.graph(),false);

    fgq.sourceComponent(nodeFilter,
                         arcFilter,
                         cutFilter);


    FlowGraphBuilder::SubGraph subgraph(graph(),
                                        nodeFilter,
                                        arcFilter);

    ListDigraph::NodeMap<bool> pixelsFilter(fgb.graph(),false);
    for(FlowGraphBuilder::SubGraph::ArcIt a(subgraph);a!=INVALID;++a)
    {
        pixelsFilter[ subgraph.source(a) ] = true;
    }
    pixelsFilter[fgb.source()] = false;
    pixelsFilter[fgb.target()] = false;


    KSpace& KImage = imageFlowData.getKSpace();
    std::vector<Z2i::SCell> pixelsInTheGraph;
    pixels(pixelsInTheGraph,
           pixelsFilter);

    for(std::vector<SCell>::const_iterator it=pixelsInTheGraph.begin();it!=pixelsInTheGraph.end();++it)
    {
        Z2i::Point p = KImage.sCoords(*it);
        updatedImage.setValue(p,255);
    }

    fillHoles(updatedImage);
}

void Flow::pixels(std::vector<SCell>& pixelsVector,
                  ListDigraph::NodeMap<bool>& nodeFilter)
{
    for(ListDigraph::NodeIt n(fgb.graph());n!=INVALID;++n)
    {
        if(nodeFilter[n]){
            Z2i::SCell pixel = fgb.pixelsMap()[n];
            pixelsVector.push_back(pixel);
        }
    }
}

void Flow::fillHoles(Image2D& out)
{
    KSpace& KImage = imageFlowData.getKSpace();
    std::vector<Z2i::SCell> pixelsInTheGraph;

    ListDigraph::NodeMap<bool> pixelsFilter(fgb.graph(),true);
    pixelsFilter[fgb.source()] = false;
    pixelsFilter[fgb.target()] = false;

    pixels(pixelsInTheGraph,
           pixelsFilter);

    for(auto it=pixelsInTheGraph.begin();it!=pixelsInTheGraph.end();++it)
    {
        Z2i::SCell pixel = *it;
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
}