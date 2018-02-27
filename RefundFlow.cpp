#include "RefundFlow.h"

void RefundFlow::matchFlows(ListDigraph::ArcMap<int>& mapping,
                            FlowGraph &f1,
                            FlowGraph &f2)
{
    std::map< FlowGraphQuery::ArcKey,ListDigraph::Arc > unifiedArcMapping;
    for(ListDigraph::ArcIt arcF1(f1.graph());arcF1!=INVALID;++arcF1)
    {
        unifiedArcMapping[FlowGraphQuery::arcKey(f1, arcF1)] = arcF1;
    }

    for(ListDigraph::ArcIt arcF2(f2.graph());arcF2!=INVALID;++arcF2)
    {
        FlowGraphQuery::ArcKey arcF2Key = FlowGraphQuery::arcKey(f2, arcF2);

        if (unifiedArcMapping.find(arcF2Key) != unifiedArcMapping.end()) {
            ListDigraph::Arc arcF1 = unifiedArcMapping[arcF2Key];
            mapping[arcF1] = f2.graph().id(arcF2);
        }
    }

}

RefundFlow::RefundFlow(ImageFlowData& imageFlowData):imageFlowData(imageFlowData),
                                                     imageOut(imageFlowData.getOriginalImage()){};





void RefundFlow::updateImage(FlowGraph& fg,
                             Image2D& updatedImage)
{

    ListDigraph::NodeMap<bool> pixelsFilter(fg.graph(),false);
    FlowGraphQuery::pixelsFilter(fg,
                                 pixelsFilter);


    KSpace& KImage = imageFlowData.getKSpace();
    std::vector<Z2i::SCell> pixelsInTheGraph;
    pixels(fg,
           pixelsInTheGraph,
           pixelsFilter);

    for(std::vector<SCell>::const_iterator it=pixelsInTheGraph.begin();it!=pixelsInTheGraph.end();++it)
    {
        Z2i::Point p = KImage.sCoords(*it);
        updatedImage.setValue(p,255);
    }

    fillHoles(fg,
              updatedImage);
}


void RefundFlow::pixels(FlowGraph& fg,
                        std::vector<SCell>& pixelsVector,
                        ListDigraph::NodeMap<bool>& nodeFilter)
{
    for(ListDigraph::NodeIt n(fg.graph());n!=INVALID;++n)
    {
        if(nodeFilter[n]){
            Z2i::SCell pixel = fg.pixel(n);
            pixelsVector.push_back(pixel);
        }
    }
}

void RefundFlow::fillHoles(FlowGraph& fg,
                           Image2D& out)
{
    KSpace& KImage = imageFlowData.getKSpace();
    std::vector<Z2i::SCell> pixelsInTheGraph;

    ListDigraph::NodeMap<bool> pixelsFilter(fg.graph(),true);
    pixelsFilter[fg.source()] = false;
    pixelsFilter[fg.target()] = false;

    pixels(fg,
           pixelsInTheGraph,
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


void RefundFlow::computeCutEnergyDifference(FlowGraph& f1,
                                            Image2D& image,
                                            ListDigraph::ArcMap<bool>& arcFilter,
                                            ListDigraph::ArcMap<double>& arcDiff)
{
    ImageFlowData imageFlowData(image);
    imageFlowData.init(f1.imageFlowData().getFlowMode(),
                       f1.imageFlowData().getGluedCurveLength());

    FlowGraphBuilder::LinelWeightMap weightMap;
    setArcsWeight(imageFlowData,weightMap);

    FlowGraph f2;
    FlowGraphBuilder fgb(f2,imageFlowData,weightMap);

    ListDigraph::ArcMap<int> arcsMatch(f1.graph(),-1);
    matchFlows(arcsMatch,f1,f2);

    for(ListDigraph::ArcIt arcF1(f1.graph());arcF1!=INVALID;++arcF1)
    {
        if(arcFilter[arcF1] && arcsMatch[arcF1]!=-1)
        {
            ListDigraph::Arc arcF2 = f2.graph().arcFromId( arcsMatch[arcF1] );
            arcDiff[arcF1] = f2.weight(arcF2) - f1.weight(arcF1);
        }
    }

}

int RefundFlow::findConflicts(FlowGraph& fg,
                              std::set<FlowGraphQuery::ArcPair>& listOfPairs,
                              FlowGraphQuery::ArcPair& pairToTest,
                              FlowGraphQuery::ArcPair conflictedPair[])
{
    ListDigraph::Arc testIntExt = pairToTest.first;
    ListDigraph::Arc testExtInt = pairToTest.second;

    int n =0;
    int conflictDistance = fg.imageFlowData().getConsecutiveGluedPairDistance();
    for(FlowGraphQuery::ArcPairIterator ait=listOfPairs.begin();ait!=listOfPairs.end();++ait)
    {
        ListDigraph::Arc intExtArc = ait->first; //intExt
        ListDigraph::Arc extIntArc = ait->second; //extInt

        int d1 = FlowGraphQuery::arcDistance(fg,testIntExt,extIntArc);
        int d2 = FlowGraphQuery::arcDistance(fg,testExtInt,intExtArc);

        if(d1<conflictDistance)
        {
            conflictedPair[n].first =testIntExt;
            conflictedPair[n].second = extIntArc;
            ++n;
        }else if(d2<conflictDistance){
            conflictedPair[n].first =testExtInt;
            conflictedPair[n].second = intExtArc;
            ++n;
        }

        if(n==2) break;

    }

    return n;

}

void RefundFlow::addConflictArcs(FlowGraph& fg,
                                 ListDigraph::Arc& a1,
                                 ListDigraph::Arc& a2,
                                 double weight)
{
    ListDigraph::Node u1 = fg.target(a1);
    ListDigraph::Node v1 = fg.source(a2);

    ListDigraph::Node u2 = fg.target(a2);
    ListDigraph::Node v2 = fg.source(a1);

    FlowGraphBuilder::addArc(fg,u1,v1,FlowGraph::ArcType::ConflictArc,weight);
    FlowGraphBuilder::addArc(fg,u2,v2,FlowGraph::ArcType::ConflictArc,weight);
}

void RefundFlow::addRefundArcs(FlowGraph& fg,
                               ListDigraph::Arc& a1,
                               ListDigraph::Arc& a2,
                               double weight)
{
    ListDigraph::Node u1 = fg.source(a1);
    ListDigraph::Node v1 = fg.target(a2);

    ListDigraph::Node u2 = fg.source(a2);
    ListDigraph::Node v2 = fg.target(a1);

    FlowGraphBuilder::addArc(fg,u1,v1,FlowGraph::ArcType::RefundArc,weight);
    FlowGraphBuilder::addArc(fg,u2,v2,FlowGraph::ArcType::RefundArc,weight);
}

int RefundFlow::runIteration(FlowGraph& fg,
                             Image2D& partialImage)
{
    FlowGraphQuery::DetourArcMap detourArcMap;
    FlowGraphQuery::detourArcMap(fg,detourArcMap);

    ListDigraph::ArcMap<bool> detourArcs(fg.graph(),false);
    FlowGraphQuery::detourArcFilter(fg,detourArcs,detourArcMap);

    ListDigraph::ArcMap<double> arcDiff(fg.graph());


    Image2D previousPartialImage = partialImage;
    updateImage(fg,partialImage);

    computeCutEnergyDifference(fg,
                               partialImage,
                               detourArcs,
                               arcDiff);

    if(partialImage==previousPartialImage)
        return Stop;


    FlowGraphQuery::ArcPairSet usedKeys;
    FlowGraphQuery::ArcPair conflictedPairs[2];
    double s;
    int n;
    for(FlowGraphQuery::DetourArcMapIterator dami = detourArcMap.begin();dami!=detourArcMap.end();++dami)
    {
        FlowGraphQuery::ArcPair key = dami->first;
        if(usedKeys.find(key)!=usedKeys.end()) continue;


        int conflicts = findConflicts(fg,
                                      usedKeys,
                                      key,
                                      conflictedPairs);
        for(int i=0;i<conflicts;++i)
        {
            addConflictArcs(fg,
                            conflictedPairs[i].first,
                            conflictedPairs[i].second,
                            10);
        }

        usedKeys.insert(key);

        const std::set<ListDigraph::Arc> &detourArcs = dami->second;

        s=0;
        n=0;
        for(auto vit=detourArcs.begin();vit!=detourArcs.end();++vit)
        {
            s+= arcDiff[*vit];
            ++n;
        }
        s/=2;

        //Flow graph arcs must not have negative weights
        if(s<0){
            addRefundArcs(fg,
                          key.first,
                          key.second,
                          0);
        }else{
            addRefundArcs(fg,
                          key.first,
                          key.second,
                          s);
        }
    }

    return Continue;
}

void RefundFlow::run()
{
    FlowGraph fg;
    LinelWeightMap weightMap;
    setArcsWeight(imageFlowData,weightMap);
    FlowGraphBuilder fgb(fg,imageFlowData,weightMap);


    double currentCutValue = cutValue(fg);
    int iteration = 1;

    Image2D partialImage = imageFlowData.getOriginalImage();
    while(runIteration(fg,partialImage)==Continue)
    {
        ++iteration;
    }

    imageOut = partialImage;
}
