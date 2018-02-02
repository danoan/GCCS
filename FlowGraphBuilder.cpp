#include "FlowGraphBuilder.h"

namespace Stats{
    int gluedEdges;
    int sourceEdges;
    int targetEdges;
    int externalEdges;
    int internalEdges;

    int makeConvexEdges;
    int escapeEdges;

    double extMakeConvexCut;
}

void FlowGraphBuilder::randomizeCurve(std::vector<Z2i::SCell>::const_iterator cBegin,
                                      std::vector<Z2i::SCell>::const_iterator cEnd,
                                      int size,
                                      std::vector<Z2i::SCell>& scellsRand)
{
    DGtal::Circulator<std::vector<Z2i::SCell>::const_iterator> C(cBegin,
                                                                 cBegin,
                                                                 cEnd);

    int inc = rand()%size;
    while(inc>0){
        ++C;
        --inc;
    }

    DGtal::Circulator<std::vector<Z2i::SCell>::const_iterator> it = C;
    do{
        scellsRand.push_back(*it);
        ++it;
    }while(it!=C);
}

Stats::MySubGraph FlowGraphBuilder::externalCurveSubgraph(ListDigraph::NodeMap<bool>& node_filter,
                                                          ListDigraph::ArcMap<bool>& arc_filter )
{
    node_filter[targetNode] = true;
    for(auto it=extCurve.begin();it!=extCurve.end();++it){
        KSpace::SCell pxDir = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
        KSpace::SCell pxInd = KImage.sIndirectIncident(*it,KImage.sOrthDir(*it));

        ListDigraph::Node n1 = coordToNode[pxDir.preCell().coordinates];
        ListDigraph::Node n2 = coordToNode[pxInd.preCell().coordinates];

        node_filter[n1] = true;
        node_filter[n2] = true;
    }

    Stats::MySubGraph msg(fg,node_filter,arc_filter);
    return msg;
}



Stats::MySubGraph FlowGraphBuilder::particularCurveSubgraph(ListDigraph::NodeMap<bool>& node_filter,
                                                          ListDigraph::ArcMap<bool>& arc_filter )
{
    node_filter[targetNode] = true;
    for(auto it=extCurve.begin();it!=extCurve.end();++it){
        KSpace::SCell pxDir = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
        KSpace::SCell pxInd = KImage.sIndirectIncident(*it,KImage.sOrthDir(*it));

        if(pxDir.preCell().coordinates[1]<200) continue;
        if(pxDir.preCell().coordinates[0]>150) continue;

        ListDigraph::Node n1 = coordToNode[pxDir.preCell().coordinates];
        ListDigraph::Node n2 = coordToNode[pxInd.preCell().coordinates];

        node_filter[n1] = true;
        node_filter[n2] = true;
    }

    Stats::MySubGraph msg(fg,node_filter,arc_filter);
    return msg;
}

Stats::MySubGraph FlowGraphBuilder::internalCurveSubgraph(ListDigraph::NodeMap<bool>& node_filter,
                                                          ListDigraph::ArcMap<bool>& arc_filter )
{
    node_filter[sourceNode] = true;
    for(auto it=intCurve.begin();it!=intCurve.end();++it){
        KSpace::SCell px = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
        ListDigraph::Node n = coordToNode[px.preCell().coordinates];
        node_filter[n] = true;
    }


    Stats::MySubGraph msg(fg,node_filter,arc_filter);
    return msg;
}


void FlowGraphBuilder::draw()
{
    graphToEps(fg,"flowGraph.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .coords(coords)
            .run();

    ListDigraph::NodeMap<bool> isg_node_filter(fg,false);
    ListDigraph::ArcMap<bool> isg_arc_filter(fg,true);

    Stats::MySubGraph isg = internalCurveSubgraph(isg_node_filter,isg_arc_filter);

    graphToEps(isg,"internalGraph.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .coords(coords)
            .run();

    ListDigraph::NodeMap<bool> esg_node_filter(fg,false);
    ListDigraph::ArcMap<bool> esg_arc_filter(fg,true);

    Stats::MySubGraph esg = externalCurveSubgraph(esg_node_filter,esg_arc_filter);

    graphToEps(esg,"externalGraph.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .coords(coords)
            .run();



    ListDigraph::NodeMap<bool> psg_node_filter(fg,false);
    ListDigraph::ArcMap<bool> psg_arc_filter(fg,true);

    Stats::MySubGraph psg = particularCurveSubgraph(psg_node_filter,psg_arc_filter);

    graphToEps(psg,"particularGraph.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .coords(coords)
            .run();
}

void FlowGraphBuilder::operator()(std::map<Z2i::SCell,double>& weightMap)
{
    Stats::internalEdges=0;
    Stats::externalEdges=0;
    Stats::gluedEdges=0;

    Stats::escapeEdges=0;
    Stats::makeConvexEdges=0;

    Stats::sourceEdges=0;
    Stats::targetEdges=0;

    Stats::extMakeConvexCut=0;



    createGridCurveEdges(intCurve.begin(),intCurve.end(),weightMap);
    Stats::internalEdges = Stats::externalEdges;

    Stats::externalEdges = 0;
    createGridCurveEdges(extCurve.begin(),extCurve.end(),weightMap);


    SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);
    SegCut::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);


    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> visitedNodes(myComp);

    createGluedCurveEdges(gcsRange.begin(),gcsRange.end(),weightMap,visitedNodes);
    createEscapeEdges(visitedNodes);

    createSourceEdges(visitedNodes);
    createTargetEdgesFromGluedSegments(gcsRange.begin(),gcsRange.end(),weightMap);
    createTargetEdgesFromExteriorCurve(visitedNodes,weightMap);


    std::cout << "External Curve Size::" << extCurve.size() << std::endl;
    std::cout << "External Edges::" << Stats::externalEdges << std::endl;


    std::cout << "Internal Curve Size::" << intCurve.size() << std::endl;
    std::cout << "Internal Edges::" << Stats::internalEdges << std::endl;


    std::cout << "Glued Edges::" << Stats::gluedEdges << std::endl;
    std::cout << "Make Convex Edges::" << Stats::makeConvexEdges << std::endl;
    std::cout << "Escape Edges::" << Stats::escapeEdges << std::endl;

    std::cout << "Source Edges::" << Stats::sourceEdges << std::endl;
    std::cout << "Target Edges::" << Stats::targetEdges << std::endl;

    std::cout << "Total Edges::" << countArcs(fg) << std::endl;


    std::cout << "External Curve + Make-Convex Cut Value::" << Stats::extMakeConvexCut << std::endl;

    
    xavg=(int) xavg/(countNodes(fg));
    yavg=(int) yavg/(countNodes(fg));

    coords[sourceNode] = LemonPoint( xavg,yavg);
    coords[targetNode] = LemonPoint( xavg-40,yavg+40);
    
}

void FlowGraphBuilder::createSourceEdges(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    typedef Curve::InnerPointsRange::ConstIterator IteratorType;

    KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
    for(IteratorType it=intCurve.getInnerPointsRange().begin();it!=intCurve.getInnerPointsRange().end();++it)
    {
        KSpace::SCell directIncidentPixel = KImage.sCell(*it,pixelModel);

        if(visitedNodes.find(directIncidentPixel)!=visitedNodes.end()) continue;
        visitedNodes.insert(directIncidentPixel);

        ListDigraph::Arc a;
        connectNodeToPixel(sourceNode,directIncidentPixel,a);

        edgeWeight[a] = 10;

        Stats::sourceEdges++;
    }
}

bool FlowGraphBuilder::createTargetEdgesFromGluedSegments(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                                          SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                                          std::map<Z2i::SCell,double>& weightMap)
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
        createFromIteratorsQueue(connectorsInterval,weightMap);
        return true;
    }

}

void FlowGraphBuilder::createTargetEdgesFromExteriorCurve(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                                          std::map<Z2i::SCell,double>& weightMap)
{
    for(auto itC = extCurve.begin();itC!=extCurve.end();itC++){
        KSpace::SCell indirectIncidentPixel = KImage.sIndirectIncident(*itC,KImage.sOrthDir(*itC));
        if(visitedNodes.find(indirectIncidentPixel)!=visitedNodes.end()) continue;

        visitedNodes.insert(indirectIncidentPixel);

        Stats::extMakeConvexCut += weightMap[*itC];
        Stats::targetEdges++;

        ListDigraph::Arc a;
        connectPixelToNode(indirectIncidentPixel,targetNode,a);
        edgeWeight[a] = 10;
    }

}

void FlowGraphBuilder::createNodeFromPixel(KSpace::SCell pixel,
                                           ListDigraph::Node& node)
{
    Z2i::Point pixelCoord = pixel.preCell().coordinates;

    if(coordToNode.find(pixelCoord)==coordToNode.end())
    {
        coordToNode[pixelCoord] = fg.addNode();

        xavg+= pixelCoord[0];
        yavg+= pixelCoord[1];
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

void FlowGraphBuilder::createEdgeFromPixels(KSpace::SCell pSource,
                                            KSpace::SCell pTarget,
                                            ListDigraph::Arc& a)
{
    ListDigraph::Node sourcePixelNode;
    ListDigraph::Node targetPixelNode;

    createNodeFromPixel(pSource,sourcePixelNode);
    createNodeFromPixel(pTarget,targetPixelNode);

    a = fg.addArc(sourcePixelNode,targetPixelNode);
}

void FlowGraphBuilder::createEdgeFromLinel(Curve::SCell& linel,
                                           KSpace& KImage,
                                           std::map<Z2i::SCell,double>& weightMap)
{
    Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
    Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

    ListDigraph::Arc a;
    createEdgeFromPixels(directPixel,
                         indirectPixel,
                         a);

    edgeWeight[a] = wFactor*weightMap[linel];
}

void FlowGraphBuilder::createGridCurveEdges(Curve::ConstIterator curveBegin,
                                            Curve::ConstIterator curveEnd,
                                            std::map<Z2i::SCell,double>& weightMap)
{
    for(auto it=curveBegin;it!=curveEnd;++it){
        Curve::SCell linel = *it;

        createEdgeFromLinel(linel,
                            KImage,
                            weightMap);

        Stats::externalEdges++;
    }
}

void FlowGraphBuilder::createGluedCurveEdges(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                             SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                             std::map<Z2i::SCell,double>& weightMap,
                                             std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    for(auto it=gluedRangeBegin;it!=gluedRangeEnd;++it){
        SegCut::SCellGluedCurveIterator begin = it->first;

        auto linkIt = begin.connectorsBegin();
        do{
            Z2i::SCell linel = *linkIt;
            if(begin.connectorType()==makeConvex){
                Stats::makeConvexEdges++;
            }else{
                Stats::gluedEdges++;
            }

            createEdgeFromLinel(linel,
                                KImage,
                                weightMap);

            visitedNodes.insert( KImage.sIndirectIncident(linel,KImage.sOrthDir(linel)) );

            if(linkIt==begin.connectorsEnd()) break;
            linkIt++;
        }while(true);

    }

}

void FlowGraphBuilder::createEscapeEdges(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes) {

    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> forbiddenPoints(myComp);
    for(auto it=intCurve.begin();it!=intCurve.end();++it)
    {
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        //I am using the sign attribute to identify inner pixels from intern or extern curve
        KImage.sSetSign(innerPixel,false);
        forbiddenPoints.insert( innerPixel );
    }


    for(auto it=extCurve.begin();it!=extCurve.end();++it)
    {
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        KImage.sSetSign(innerPixel,true);
        forbiddenPoints.insert( innerPixel );
    }


    for(Curve::ConstIterator it = intCurve.begin();it!=intCurve.end();++it){
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

    visitedNodes.insert(candidate);
    SCells N = KImage.sProperNeighborhood(candidate);
    for(SCells::ConstIterator n=N.begin();n!=N.end();++n){

        if( borderPoints.find(*n)!=borderPoints.end() ){
            KSpace::SCell bPoint = *borderPoints.find(*n);
            if(KImage.sSign(bPoint)){
                ListDigraph::Arc a1,a2;
                createEdgeFromPixels(candidate,*n,a1);
                createEdgeFromPixels(*n,candidate,a2);

                edgeWeight[a1]=10;
                edgeWeight[a2]=10;

                Stats::escapeEdges+=2;
            }

            continue;
        }

        explore(*n,borderPoints,visitedNodes);

        ListDigraph::Arc a1,a2;
        createEdgeFromPixels(candidate,*n,a1);
        createEdgeFromPixels(*n,candidate,a2);

        edgeWeight[a1]=10;
        edgeWeight[a2]=10;

        Stats::escapeEdges+=2;
    }
}
