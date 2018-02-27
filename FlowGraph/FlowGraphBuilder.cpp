#include "FlowGraphBuilder.h"

FlowGraphBuilder::FlowGraphBuilder(FlowGraph& fg,
                                   ImageFlowData& imageFlowData,
                                   LinelWeightMap& weightMap):imageFlowData(imageFlowData),
                                                              weightMap(weightMap),
                                                              fg(fg)
{
    fg.gluedCurveLength = imageFlowData.getGluedCurveLength();
    fg.consecutiveGluedPairDistance = imageFlowData.getConsecutiveGluedPairDistance();
    fg.diffDistance = imageFlowData.getDiffDistance();


    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        createCurveArcs(it->curve.begin(),it->curve.end(),it->curveType,weightMap);
    }


    UnsignedSCellComparison myComp;
    UnsignedSCellSet visitedNodes(myComp);

    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        createGluedArcs(it->gcsRangeBegin(),
                        it->gcsRangeEnd(),
                        weightMap,
                        visitedNodes);

        createEscapeArcs(it->intCurveData.curve,
                         it->extCurveData.curve,
                         visitedNodes);
    }

    createSourceArcs(imageFlowData.getMostInnerCurve(),
                     visitedNodes);


    for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
    {
        createTargetArcsFromGluedSegments(it->gcsRangeBegin(),
                                          it->gcsRangeEnd(),
                                          weightMap,
                                          visitedNodes);
    }

    createTargetArcsFromExteriorCurve(imageFlowData.getMostOuterCurve(),
                                      visitedNodes);

    setTerminalsCoordinates();

}


void FlowGraphBuilder::createSourceArcs(Curve& erodedCurve,
                                        UnsignedSCellSet& visitedNodes)
{
    typedef Curve::InnerPointsRange::ConstIterator IteratorType;

    KSpace& KImage = imageFlowData.getKSpace();

    KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
    for(IteratorType it=erodedCurve.getInnerPointsRange().begin();it!=erodedCurve.getInnerPointsRange().end();++it)
    {
        KSpace::SCell directIncidentPixel = KImage.sCell(*it,pixelModel);

        if(visitedNodes.find(directIncidentPixel)!=visitedNodes.end()) continue;
        visitedNodes.insert(directIncidentPixel);

        ListDigraph::Arc a;
        connectNodeToPixel(fg.sourceNode,directIncidentPixel,a);

        fg.arcWeightMap[a] = 10;
        fg.arcTypeMap[a] = FlowGraph::ArcType::SourceArc;
    }
}

void FlowGraphBuilder::createTargetArcsFromGluedSegments(ImageFlowData::GluedCurveIteratorPair gluedRangeBegin,
                                                         ImageFlowData::GluedCurveIteratorPair gluedRangeEnd,
                                                         LinelWeightMap& weightMap,
                                                         UnsignedSCellSet& visitedNodes)
{
    typedef Curve::OuterPointsRange::ConstIterator IteratorType;
    typedef ImageFlowData::SCellGluedCurveIterator::LinkIteratorType LinkIteratorType;
    typedef ImageFlowData::SCellGluedCurveIterator::CurveCirculator CurveCirculator;


    CurveCirculator start,end;
    std::queue<LinkIteratorType> connectorsInterval;

    for(ImageFlowData::GluedCurveIteratorPair gcIt = gluedRangeBegin;gcIt!=gluedRangeEnd;++gcIt){
        if(gcIt->first.connectorType()!=makeConvex) continue;
        connectorsInterval.push(gcIt->first.connectorsBegin());
        connectorsInterval.push(gcIt->first.connectorsEnd());
    }

    if(!connectorsInterval.empty()){
        createFromIteratorsQueue(connectorsInterval,weightMap,visitedNodes);
    }

}

void FlowGraphBuilder::createTargetArcsFromExteriorCurve(Curve& extCurve,
                                                         UnsignedSCellSet& visitedNodes)
{
    KSpace& KImage = imageFlowData.getKSpace();

    for(auto itC = extCurve.begin();itC!=extCurve.end();itC++){
        KSpace::SCell indirectIncidentPixel = KImage.sIndirectIncident(*itC,KImage.sOrthDir(*itC));

        if(visitedNodes.find(indirectIncidentPixel)!=visitedNodes.end()){
            continue;
        }

        visitedNodes.insert(indirectIncidentPixel);


        ListDigraph::Arc a;
        connectPixelToNode(indirectIncidentPixel,fg.targetNode,a);

        fg.arcWeightMap[a] = 10;
        fg.arcTypeMap[a] = FlowGraph::ArcType::TargetArc;
    }

}

void FlowGraphBuilder::createNodeFromPixel(KSpace::SCell pixel,
                                           ListDigraph::Node& node)
{
    Z2i::Point pixelCoord = pixel.preCell().coordinates;

    if(fg.coordToNode.find(pixelCoord)==fg.coordToNode.end())
    {
        fg.coordToNode[pixelCoord] = fg.digraph.addNode();
    }

    node = fg.coordToNode[pixelCoord];
    fg.pixelMap[ node ] = pixel;

    fg.coords[node] = LemonPoint( pixelCoord[0],
                                  pixelCoord[1]);
}

void FlowGraphBuilder::connectNodeToPixel(ListDigraph::Node& sourceNode,
                                          KSpace::SCell pTarget,
                                          ListDigraph::Arc& a)
{
    ListDigraph::Node targetNode;
    createNodeFromPixel(pTarget,targetNode);

    a = fg.digraph.addArc(sourceNode,targetNode);
}

void FlowGraphBuilder::connectPixelToNode(KSpace::SCell pSource,
                                          ListDigraph::Node& targetNode,
                                          ListDigraph::Arc& a)
{
    ListDigraph::Node sourceNode;
    createNodeFromPixel(pSource,sourceNode);

    a = fg.digraph.addArc(sourceNode,targetNode);
}

void FlowGraphBuilder::createArcFromPixels(ListDigraph::Arc& arc,
                                           KSpace::SCell pSource,
                                           KSpace::SCell pTarget,
                                           FlowGraph::ArcType at,
                                           double weight)
{
    ListDigraph::Node sourcePixelNode;
    ListDigraph::Node targetPixelNode;

    createNodeFromPixel(pSource,sourcePixelNode);
    createNodeFromPixel(pTarget,targetPixelNode);

    addArc(arc,sourcePixelNode,targetPixelNode,at,weight);
}

void FlowGraphBuilder::createArcFromPixels(KSpace::SCell pSource,
                                           KSpace::SCell pTarget,
                                           FlowGraph::ArcType at,
                                           double weight)
{
    ListDigraph::Node sourcePixelNode;
    ListDigraph::Node targetPixelNode;

    createNodeFromPixel(pSource,sourcePixelNode);
    createNodeFromPixel(pTarget,targetPixelNode);

    ListDigraph::Arc arc;
    addArc(arc,sourcePixelNode,targetPixelNode,at,weight);
}

void FlowGraphBuilder::createArcFromLinel(ListDigraph::Arc& arc,
                                          Curve::SCell& linel,
                                          LinelWeightMap& weightMap,
                                          FlowGraph::ArcType at,
                                          bool invert)
{
    KSpace& KImage = imageFlowData.getKSpace();

    Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
    Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

    if(invert){
        createArcFromPixels(arc,
                            indirectPixel,
                            directPixel,
                            at,
                            weightMap[linel]);
    }else{
        createArcFromPixels(arc,
                            directPixel,
                            indirectPixel,
                            at,
                            weightMap[linel]);
    }

    fg.arcSCellMap[arc] = linel;
    fg.scellArcMap[linel] = arc;
}

void FlowGraphBuilder::createArcFromLinel(Curve::SCell& linel,
                                          LinelWeightMap& weightMap,
                                          FlowGraph::ArcType at,
                                          bool invert)
{
    ListDigraph::Arc arc;
    createArcFromLinel(arc,
                       linel,
                       weightMap,
                       at,
                       invert);
}


void FlowGraphBuilder::createCurveArcs(Curve::ConstIterator curveBegin,
                                     Curve::ConstIterator curveEnd,
                                     ImageFlowData::CurveType ct,
                                     LinelWeightMap& weightMap)
{
    FlowGraph::ArcType at;
    if( ct==ImageFlowData::CurveType::OriginalCurve ){
        at = FlowGraph::ArcType::InternalCurveArc;
    }else{
        at = FlowGraph::ArcType::ExternalCurveArc;
    }

    for(auto it=curveBegin;it!=curveEnd;++it){
        Curve::SCell linel = *it;

        createArcFromLinel(linel,
                           weightMap,
                           at,
                           false);
    }
}

void FlowGraphBuilder::createGluedArcs(ImageFlowData::GluedCurveIteratorPair gluedRangeBegin,
                                       ImageFlowData::GluedCurveIteratorPair gluedRangeEnd,
                                       LinelWeightMap& weightMap,
                                       UnsignedSCellSet& visitedNodes)
{
    KSpace& KImage = imageFlowData.getKSpace();

    for(auto it=gluedRangeBegin;it!=gluedRangeEnd;++it){
        ImageFlowData::SCellGluedCurveIterator begin = it->first;

        auto linkIt = begin.connectorsBegin();
        do{
            Z2i::SCell linel = *linkIt;

            if(begin.connectorType()==makeConvex){
                createArcFromLinel(linel,
                                   weightMap,
                                   FlowGraph::ArcType::MakeConvexArc,
                                   false);
            }else{
                ListDigraph::Arc a;
                FlowGraph::ArcType at = begin.connectorType()==ConnectorType::internToExtern?
                                        FlowGraph::IntExtGluedArc:
                                        FlowGraph::ExtIntGluedArc;
                createArcFromLinel(a,
                                   linel,
                                   weightMap,
                                   at,
                                   Development::invertGluedArcs);


                fg.arcCirculatorMap[a] = CirculatorPair(begin.curveSegment1End(),
                                                     begin.curveSegment2Begin());
            }


            visitedNodes.insert( KImage.sIndirectIncident(linel,KImage.sOrthDir(linel)) );

            if(linkIt==begin.connectorsEnd()) break;
            linkIt++;
        }while(true);

    }

}

void FlowGraphBuilder::createEscapeArcs(Curve& fromCurve,
                                        Curve& toCurve,
                                        UnsignedSCellSet& visitedNodes)
{
    KSpace& KImage = imageFlowData.getKSpace();

    UnsignedSCellComparison myComp;
    UnsignedSCellSet forbiddenPoints(myComp);

    for(auto it=fromCurve.begin();it!=fromCurve.end();++it){
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        //I am using the sign attribute to identify inner pixels from intern or extern curve
        KImage.sSetSign(innerPixel,false);
        forbiddenPoints.insert( innerPixel );
    }


    for(auto it=toCurve.begin();it!=toCurve.end();++it)
    {
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        KImage.sSetSign(innerPixel,true);
        forbiddenPoints.insert( innerPixel );
    }


    for(Curve::ConstIterator it = fromCurve.begin();it!=fromCurve.end();++it){
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell candidate = KImage.sIndirectIncident(*it,orthDir); //OutterPixel

        if( forbiddenPoints.find(candidate)!=forbiddenPoints.end() ) continue;

        ListDigraph::Arc in,out;
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        explore(candidate,forbiddenPoints,visitedNodes);
    }



}


void FlowGraphBuilder::explore(KSpace::SCell candidate,
                               UnsignedSCellSet& borderPoints,
                               UnsignedSCellSet& visitedNodes)
{

    if(visitedNodes.find(candidate)!=visitedNodes.end()) return;

    KSpace& KImage = imageFlowData.getKSpace();

    visitedNodes.insert(candidate);
    SCells N = KImage.sProperNeighborhood(candidate);
    for(SCells::ConstIterator n=N.begin();n!=N.end();++n){

        if( borderPoints.find(*n)!=borderPoints.end() ){
            KSpace::SCell bPoint = *borderPoints.find(*n);
            if(KImage.sSign(bPoint)){
                createArcFromPixels(candidate,*n,FlowGraph::ArcType::EscapeArc,10);
                createArcFromPixels(*n,candidate,FlowGraph::ArcType::EscapeArc,10);
            }

            continue;
        }

        explore(*n,borderPoints,visitedNodes);

        createArcFromPixels(candidate,*n,FlowGraph::ArcType::EscapeArc,10);
        createArcFromPixels(*n,candidate,FlowGraph::ArcType::EscapeArc,10);
    }

}


template<typename TType>
int FlowGraphBuilder::createFromIteratorsQueue(std::queue<TType> intervals,
                                               LinelWeightMap& weightMap,
                                               UnsignedSCellSet& visitedNodes)
{

    KSpace& KImage = imageFlowData.getKSpace();

    TType start,end;
    KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
    while(!intervals.empty()){
        start = intervals.front();
        intervals.pop();
        end = intervals.front();
        intervals.pop();

        TType itC = start;
        do{

            KSpace::SCell indirectIncidentPixel = KImage.sIndirectIncident(*itC,KImage.sOrthDir(*itC));
            KSpace::SCell directIncidentPixel = KImage.sDirectIncident(*itC,KImage.sOrthDir(*itC));

            visitedNodes.insert(directIncidentPixel);

            ListDigraph::Node u;
            createNodeFromPixel(indirectIncidentPixel,u);

            addArc(u,
                   fg.targetNode,
                   FlowGraph::ArcType::TargetArc,
                   10);


            if(itC==end) break;
            ++itC;
        }while(true);

    }

}

void FlowGraphBuilder::setTerminalsCoordinates()
{
    double x=0;
    double y=0;
    double i=0;
    for(ListDigraph::NodeIt n(fg.digraph);n!=INVALID;++n)
    {
        if(n==fg.sourceNode || n==fg.targetNode) continue;

        x+=fg.coords[n][0];
        y+=fg.coords[n][1];
        i+=1;
    }
    x/=i;
    y/=i;

    fg.coords[fg.sourceNode] = LemonPoint(x,y);
    fg.coords[fg.targetNode] = LemonPoint(x+40,y+40);
}

void FlowGraphBuilder::addArc(ListDigraph::Arc& arc,
                              ListDigraph::Node& u,
                              ListDigraph::Node& v,
                              FlowGraph::ArcType at,
                              double weight)
{
    arc = fg.digraph.addArc(u,v);
    fg.arcWeightMap[arc] = weight;
    fg.arcTypeMap[arc] = at;
}

void FlowGraphBuilder::addArc(ListDigraph::Node& u,
                              ListDigraph::Node& v,
                              FlowGraph::ArcType at,
                              double weight)
{
    ListDigraph::Arc a1 = fg.digraph.addArc(u,v);
    fg.arcWeightMap[a1] = weight;
    fg.arcTypeMap[a1] = at;
}

void FlowGraphBuilder::addArc(FlowGraph& fg,
                              ListDigraph::Node& u,
                              ListDigraph::Node& v,
                              FlowGraph::ArcType at,
                              double weight)
{
    ListDigraph::Arc a1 = fg.digraph.addArc(u,v);
    fg.arcWeightMap[a1] = weight;
    fg.arcTypeMap[a1] = at;
}
