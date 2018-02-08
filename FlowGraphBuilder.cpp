#include "FlowGraphBuilder.h"


void FlowGraphBuilder::operator()(LinelWeightMap& weightMap)
{
    for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
    {
        createCurveArcs(it->curve.begin(),it->curve.end(),weightMap);
    }


    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> visitedNodes(myComp);

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


}


void FlowGraphBuilder::createSourceArcs(Curve& erodedCurve,
                                        std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
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
        connectNodeToPixel(sourceNode,directIncidentPixel,a);

        arcWeight[a] = 10;
        arcType[a] = ArcType::SourceArc;
    }
}

void FlowGraphBuilder::createTargetArcsFromGluedSegments(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                                        SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                                        std::map<Z2i::SCell,double>& weightMap,
                                                        std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    typedef Curve::OuterPointsRange::ConstIterator IteratorType;
    typedef SegCut::SCellGluedCurveIterator::LinkIteratorType LinkIteratorType;
    typedef SegCut::SCellGluedCurveIterator::CurveCirculator CurveCirculator;


    CurveCirculator start,end;
    std::queue<LinkIteratorType> connectorsInterval;

    for(SegCut::GluedCurveIteratorPair gcIt = gluedRangeBegin;gcIt!=gluedRangeEnd;++gcIt){
        if(gcIt->first.connectorType()!=makeConvex) continue;
        connectorsInterval.push(gcIt->first.connectorsBegin());
        connectorsInterval.push(gcIt->first.connectorsEnd());
    }

    if(!connectorsInterval.empty()){
        createFromIteratorsQueue(connectorsInterval,weightMap,visitedNodes);
    }

}

void FlowGraphBuilder::createTargetArcsFromExteriorCurve(Curve& extCurve,
                                                        std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    KSpace& KImage = imageFlowData.getKSpace();

    for(auto itC = extCurve.begin();itC!=extCurve.end();itC++){
        KSpace::SCell indirectIncidentPixel = KImage.sIndirectIncident(*itC,KImage.sOrthDir(*itC));

        if(visitedNodes.find(indirectIncidentPixel)!=visitedNodes.end()){
            continue;
        }

        visitedNodes.insert(indirectIncidentPixel);


        ListDigraph::Arc a;
        connectPixelToNode(indirectIncidentPixel,targetNode,a);

        arcWeight[a] = 10;
        arcType[a] = ArcType::TargetArc;
    }

}

void FlowGraphBuilder::createNodeFromPixel(KSpace::SCell pixel,
                                           ListDigraph::Node& node)
{
    Z2i::Point pixelCoord = pixel.preCell().coordinates;

    if(coordToNode.find(pixelCoord)==coordToNode.end())
    {
        coordToNode[pixelCoord] = fg.addNode();
    }

    node = coordToNode[pixelCoord];
    pixelMap[ node ] = pixel;

    coords[node] = LemonPoint( pixelCoord[0],
                               pixelCoord[1]);
}

void FlowGraphBuilder::connectNodeToPixel(ListDigraph::Node& sourceNode,
                                          KSpace::SCell pTarget,
                                          ListDigraph::Arc& a)
{
    ListDigraph::Node targetNode;
    createNodeFromPixel(pTarget,targetNode);

    a = fg.addArc(sourceNode,targetNode);
}

void FlowGraphBuilder::connectPixelToNode(KSpace::SCell pSource,
                                          ListDigraph::Node& targetNode,
                                          ListDigraph::Arc& a)
{
    ListDigraph::Node sourceNode;
    createNodeFromPixel(pSource,sourceNode);

    a = fg.addArc(sourceNode,targetNode);
}

void FlowGraphBuilder::createArcFromPixels(KSpace::SCell pSource,
                                           KSpace::SCell pTarget,
                                           ListDigraph::Arc& a)
{
    ListDigraph::Node sourcePixelNode;
    ListDigraph::Node targetPixelNode;

    createNodeFromPixel(pSource,sourcePixelNode);
    createNodeFromPixel(pTarget,targetPixelNode);

    a = fg.addArc(sourcePixelNode,targetPixelNode);
}

void FlowGraphBuilder::createArcFromLinel(Curve::SCell& linel,
                                          std::map<Z2i::SCell,double>& weightMap,
                                          ArcType at,
                                          bool invert)
{
    KSpace& KImage = imageFlowData.getKSpace();

    Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
    Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

    ListDigraph::Arc a;
    if(invert){
        createArcFromPixels(indirectPixel,
                            directPixel,
                            a);
    }else{
        createArcFromPixels(directPixel,
                            indirectPixel,
                            a);
    }


    arcWeight[a] = weightMap[linel];
    arcType[a] = at;
}

void FlowGraphBuilder::createCurveArcs(Curve::ConstIterator curveBegin,
                                     Curve::ConstIterator curveEnd,
                                     std::map<Z2i::SCell,double>& weightMap)
{
    for(auto it=curveBegin;it!=curveEnd;++it){
        Curve::SCell linel = *it;

        createArcFromLinel(linel,
                           weightMap,
                           ArcType::CurveArc);
    }
}

void FlowGraphBuilder::createGluedArcs(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                       SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                       std::map<Z2i::SCell,double>& weightMap,
                                       std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    KSpace& KImage = imageFlowData.getKSpace();

    for(auto it=gluedRangeBegin;it!=gluedRangeEnd;++it){
        SegCut::SCellGluedCurveIterator begin = it->first;

        auto linkIt = begin.connectorsBegin();
        do{
            Z2i::SCell linel = *linkIt;

            if(begin.connectorType()==makeConvex){
                createArcFromLinel(linel,
                                   weightMap,
                                   ArcType::MakeConvexArc,
                                   false);
            }else{
                createArcFromLinel(linel,
                                   weightMap,
                                   ArcType::GluedArc,
                                   Development::invertGluedArcs);
            }


            visitedNodes.insert( KImage.sIndirectIncident(linel,KImage.sOrthDir(linel)) );

            if(linkIt==begin.connectorsEnd()) break;
            linkIt++;
        }while(true);

    }

}

void FlowGraphBuilder::createEscapeArcs(Curve& fromCurve,
                                       Curve& toCurve,
                                       std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    KSpace& KImage = imageFlowData.getKSpace();

    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> forbiddenPoints(myComp);

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
                               std::set<KSpace::SCell,UnsignedSCellComparison>& borderPoints,
                               std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{

    if(visitedNodes.find(candidate)!=visitedNodes.end()) return;

    KSpace& KImage = imageFlowData.getKSpace();

    visitedNodes.insert(candidate);
    SCells N = KImage.sProperNeighborhood(candidate);
    for(SCells::ConstIterator n=N.begin();n!=N.end();++n){

        if( borderPoints.find(*n)!=borderPoints.end() ){
            KSpace::SCell bPoint = *borderPoints.find(*n);
            if(KImage.sSign(bPoint)){
                ListDigraph::Arc a1,a2;
                createArcFromPixels(candidate,*n,a1);
                createArcFromPixels(*n,candidate,a2);

                arcWeight[a1]=10;
                arcType[a1] = ArcType::EscapeArc;

                arcWeight[a2]=10;
                arcType[a2] = ArcType::EscapeArc;

            }

            continue;
        }

        explore(*n,borderPoints,visitedNodes);

        ListDigraph::Arc a1,a2;
        createArcFromPixels(candidate,*n,a1);
        createArcFromPixels(*n,candidate,a2);

        arcWeight[a1]=10;
        arcType[a1] = ArcType::EscapeArc;

        arcWeight[a2]=10;
        arcType[a2] = ArcType::EscapeArc;
    }

}


template<typename TType>
int FlowGraphBuilder::createFromIteratorsQueue(std::queue<TType> intervals,std::map<Z2i::SCell, double>& weightMap,
                                               std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
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

            ListDigraph::Arc a;
            connectPixelToNode(indirectIncidentPixel,targetNode,a);

            arcWeight[a] = 10; //This line didn't exist before refactoring
            arcType[a] = ArcType::TargetArc;

            if(itC==end) break;
            ++itC;
        }while(true);

    }

}