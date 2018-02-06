#include "FlowGraphBuilder.h"

namespace Stats{
    int gluedEdges;
    int sourceEdges;
    int targetEdges;
    int erodedEdges;
    int originalEdges;
    int dilatedEdges;

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

Stats::MySubGraph FlowGraphBuilder::curveSubgraph(Curve& curve,
                                                  ListDigraph::NodeMap<bool>& node_filter,
                                                  ListDigraph::ArcMap<bool>& arc_filter )
{
    for(auto it=curve.begin();it!=curve.end();++it){
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

    Curve& firstCurve = curvesVector[0];
    isg_node_filter[sourceNode]=true;
    isg_node_filter[targetNode]=false;
    Stats::MySubGraph isg = curveSubgraph(firstCurve,isg_node_filter,isg_arc_filter);

    graphToEps(isg,"internalGraph.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .coords(coords)
            .run();

    ListDigraph::NodeMap<bool> esg_node_filter(fg,false);
    ListDigraph::ArcMap<bool> esg_arc_filter(fg,true);

    Curve& lastCurve = curvesVector[curvesVector.size()-1];
    esg_node_filter[sourceNode]=false;
    esg_node_filter[targetNode]=true;
    Stats::MySubGraph esg = curveSubgraph(lastCurve,esg_node_filter,esg_arc_filter);

    graphToEps(esg,"externalGraph.eps")
            .nodeScale(.0005)
            .arcWidthScale(0.001)
            .drawArrows()
            .arrowWidth(0.1)
            .coords(coords)
            .run();
}

void FlowGraphBuilder::operator()(std::map<Z2i::SCell,double>& weightMap, bool multiplePair)
{
    Stats::erodedEdges=0;
    Stats::originalEdges=0;
    Stats::dilatedEdges=0;
    Stats::gluedEdges=0;

    Stats::escapeEdges=0;
    Stats::makeConvexEdges=0;

    Stats::sourceEdges=0;
    Stats::targetEdges=0;

    Stats::extMakeConvexCut=0;

    Stats::erodedEdges = createGridCurveEdges(curvesVector[0].begin(),curvesVector[0].end(),weightMap);
    Stats::originalEdges = createGridCurveEdges(curvesVector[1].begin(),curvesVector[1].end(),weightMap);
    Stats::dilatedEdges = createGridCurveEdges(curvesVector[2].begin(),curvesVector[2].end(),weightMap);


    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> visitedNodes(myComp);
    SegCut::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);

    for(auto itP=curvesPairs.begin();itP!=curvesPairs.end();++itP){
        Curve& fromCurve = curvesVector[itP->first];
        Curve& toCurve = curvesVector[itP->second];

        SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,fromCurve,toCurve);
        SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                             seedRange.end(),
                                             stgcF);

        createGluedCurveEdges(gcsRange.begin(),
                              gcsRange.end(),
                              weightMap,
                              visitedNodes);

        Stats::escapeEdges += createEscapeEdges(visitedNodes,
                                                fromCurve,
                                                toCurve);
    }

    Curve& firstCurve = curvesVector[0];
    Curve& beforeLastCurve = curvesVector[curvesVector.size()-2];
    Curve& lastCurve = curvesVector[curvesVector.size()-1];

    Stats::sourceEdges = createSourceEdges(firstCurve,visitedNodes);



    SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,beforeLastCurve,lastCurve);
    SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                         seedRange.end(),
                                         stgcF);

    Stats::targetEdges += createTargetEdgesFromGluedSegments(visitedNodes,gcsRange.begin(),gcsRange.end(),weightMap);
    Stats::targetEdges += createTargetEdgesFromExteriorCurve(lastCurve,visitedNodes,weightMap);



    for(auto it=curvesVector.begin();it!=curvesVector.end();++it){
        std::cout << "Curve size::" << it->size() << std::endl;
    }

    std::cout << "Eroded Edges::" << Stats::erodedEdges << std::endl;
    std::cout << "Original Edges::" << Stats::originalEdges << std::endl;
    std::cout << "Dilated Edges::" << Stats::dilatedEdges << std::endl;



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

void FlowGraphBuilder::operator()(std::map<Z2i::SCell,double>& weightMap)
{
    Stats::erodedEdges=0;
    Stats::originalEdges=0;
    Stats::dilatedEdges=0;
    Stats::gluedEdges=0;

    Stats::escapeEdges=0;
    Stats::makeConvexEdges=0;

    Stats::sourceEdges=0;
    Stats::targetEdges=0;

    Stats::extMakeConvexCut=0;

    Stats::originalEdges = createGridCurveEdges(curvesVector[0].begin(),curvesVector[0].end(),weightMap);
    Stats::dilatedEdges = createGridCurveEdges(curvesVector[1].begin(),curvesVector[1].end(),weightMap);


    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> visitedNodes(myComp);
    SegCut::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);

    for(auto itP=curvesPairs.begin();itP!=curvesPairs.end();++itP){
        Curve& fromCurve = curvesVector[itP->first];
        Curve& toCurve = curvesVector[itP->second];

        SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,fromCurve,toCurve);
        SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                             seedRange.end(),
                                             stgcF);

        createGluedCurveEdges(gcsRange.begin(),
                              gcsRange.end(),
                              weightMap,
                              visitedNodes);

        Stats::escapeEdges += createEscapeEdges(visitedNodes,
                                                fromCurve,
                                                toCurve);
    }

    Curve& firstCurve = curvesVector[0];
    Curve& beforeLastCurve = curvesVector[curvesVector.size()-2];
    Curve& lastCurve = curvesVector[curvesVector.size()-1];

    Stats::sourceEdges = createSourceEdges(firstCurve,visitedNodes);



    SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,beforeLastCurve,lastCurve);
    SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                         seedRange.end(),
                                         stgcF);

    Stats::targetEdges += createTargetEdgesFromGluedSegments(visitedNodes,gcsRange.begin(),gcsRange.end(),weightMap);
    Stats::targetEdges += createTargetEdgesFromExteriorCurve(lastCurve,visitedNodes,weightMap);



    for(auto it=curvesVector.begin();it!=curvesVector.end();++it){
        std::cout << "Curve size::" << it->size() << std::endl;
    }

    std::cout << "Eroded Edges::" << Stats::erodedEdges << std::endl;
    std::cout << "Original Edges::" << Stats::originalEdges << std::endl;
    std::cout << "Dilated Edges::" << Stats::dilatedEdges << std::endl;



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

int FlowGraphBuilder::createSourceEdges(Curve& erodedCurve,
                                        std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    typedef Curve::InnerPointsRange::ConstIterator IteratorType;

    int c =0;
    KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
    for(IteratorType it=erodedCurve.getInnerPointsRange().begin();it!=erodedCurve.getInnerPointsRange().end();++it)
    {
        KSpace::SCell directIncidentPixel = KImage.sCell(*it,pixelModel);

        if(visitedNodes.find(directIncidentPixel)!=visitedNodes.end()) continue;
        visitedNodes.insert(directIncidentPixel);

        ListDigraph::Arc a;
        connectNodeToPixel(sourceNode,directIncidentPixel,a);

        edgeWeight[a] = 10;

        ++c;
    }

    return c;
}

int FlowGraphBuilder::createTargetEdgesFromGluedSegments(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                                         SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                                          SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                                          std::map<Z2i::SCell,double>& weightMap)
{
    typedef Curve::OuterPointsRange::ConstIterator IteratorType;
    typedef SegCut::SCellGluedCurveIterator::LinkIteratorType LinkIteratorType;
    typedef SegCut::SCellGluedCurveIterator::CurveCirculator CurveCirculator;


    int c=0;
    CurveCirculator start,end;
    std::queue<LinkIteratorType> connectorsInterval;

    for(SegCut::GluedCurveIteratorPair gcIt = gluedRangeBegin;gcIt!=gluedRangeEnd;++gcIt){
        if(gcIt->first.connectorType()!=makeConvex) continue;
        connectorsInterval.push(gcIt->first.connectorsBegin());
        connectorsInterval.push(gcIt->first.connectorsEnd());
    }

    if(!connectorsInterval.empty()){
        createFromIteratorsQueue(visitedNodes,connectorsInterval,weightMap);
    }

    return c;

}

int FlowGraphBuilder::createTargetEdgesFromExteriorCurve(Curve& extCurve,
                                                         std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                                         std::map<Z2i::SCell,double>& weightMap)
{
    int c =0;
    for(auto itC = extCurve.begin();itC!=extCurve.end();itC++){
        KSpace::SCell indirectIncidentPixel = KImage.sIndirectIncident(*itC,KImage.sOrthDir(*itC));
        Stats::extMakeConvexCut += weightMap[*itC];

        if(visitedNodes.find(indirectIncidentPixel)!=visitedNodes.end()){
            continue;
        }

        visitedNodes.insert(indirectIncidentPixel);


        ++c;

        ListDigraph::Arc a;
        connectPixelToNode(indirectIncidentPixel,targetNode,a);
        edgeWeight[a] = 10;
    }
    return c;

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
                                           std::map<Z2i::SCell,double>& weightMap,
                                           bool invert)
{
    Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
    Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

    ListDigraph::Arc a;
    if(invert){
        createEdgeFromPixels(indirectPixel,
                             directPixel,
                             a);
    }else{
        createEdgeFromPixels(directPixel,
                             indirectPixel,
                             a);
    }


    edgeWeight[a] = wFactor*weightMap[linel];
}

int FlowGraphBuilder::createGridCurveEdges(Curve::ConstIterator curveBegin,
                                            Curve::ConstIterator curveEnd,
                                            std::map<Z2i::SCell,double>& weightMap)
{
    int c =0;
    for(auto it=curveBegin;it!=curveEnd;++it){
        Curve::SCell linel = *it;

        createEdgeFromLinel(linel,
                            KImage,
                            weightMap);

        ++c;
    }

    return c;
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

            if(begin.connectorType()==makeConvex){
                createEdgeFromLinel(linel,
                                    KImage,
                                    weightMap,
                                    false);
            }else{
                createEdgeFromLinel(linel,
                                    KImage,
                                    weightMap,
                                    true);
            }


            visitedNodes.insert( KImage.sIndirectIncident(linel,KImage.sOrthDir(linel)) );

            if(linkIt==begin.connectorsEnd()) break;
            linkIt++;
        }while(true);

    }

}

int FlowGraphBuilder::createEscapeEdges(std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes,
                                         Curve& fromCurve,
                                         Curve& toCurve)
{

    UnsignedSCellComparison myComp;
    std::set<KSpace::SCell,UnsignedSCellComparison> forbiddenPoints(myComp);

    int c =0;
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

        c+=explore(candidate,forbiddenPoints,visitedNodes);
    }

    return c;

}


int FlowGraphBuilder::explore(KSpace::SCell candidate,
                               std::set<KSpace::SCell,UnsignedSCellComparison>& borderPoints,
                               std::set<KSpace::SCell,UnsignedSCellComparison>& visitedNodes)
{
    int c = 0;
    if(visitedNodes.find(candidate)!=visitedNodes.end()) return c;

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

                c+=2;
            }

            continue;
        }

        c += explore(*n,borderPoints,visitedNodes);

        ListDigraph::Arc a1,a2;
        createEdgeFromPixels(candidate,*n,a1);
        createEdgeFromPixels(*n,candidate,a2);

        edgeWeight[a1]=10;
        edgeWeight[a2]=10;

        c+=2;
    }

    return c;
}
