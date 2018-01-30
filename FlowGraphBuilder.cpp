#include "FlowGraphBuilder.h"

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


void FlowGraphBuilder::draw()
{
    graphToEps(fg,"flowGraph.eps")
            .coords(coords)
            .run();    
}

void FlowGraphBuilder::operator()(std::map<Z2i::SCell,double>& weightMap)
{
    createGridCurveEdges(intCurve.begin(),intCurve.end(),weightMap);
    createGridCurveEdges(extCurve.begin(),extCurve.end(),weightMap);


    SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);
    SegCut::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);

    createGluedCurveEdges(gcsRange.begin(),gcsRange.end(),weightMap);
    createEscapeEdges();

    createSourceEdges();
    createTargetEdges();

    
    xavg=(int) xavg/(countNodes(fg));
    yavg=(int) yavg/(countNodes(fg));

    xavg = 160;
    yavg = 160;

    coords[sourceNode] = LemonPoint( xavg,yavg);
    coords[targetNode] = LemonPoint( xavg-40,yavg-40);
    
}

void FlowGraphBuilder::createSourceEdges()
{
    typedef Curve::InnerPointsRange::ConstIterator IteratorType;

    KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
    for(IteratorType it=intCurve.getInnerPointsRange().begin();it!=intCurve.getInnerPointsRange().end();++it)
    {
        KSpace::SCell directIncidentPixel = KImage.sCell(*it,pixelModel);
        ListDigraph::Arc a;
        connectNodeToPixel(sourceNode,directIncidentPixel,a);
        edgeWeight[a] = 10;
    }
}

void FlowGraphBuilder::createTargetEdges()
{
    typedef Curve::OuterPointsRange::ConstIterator IteratorType;

    KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
    for(IteratorType it=extCurve.getOuterPointsRange().begin();it!=extCurve.getOuterPointsRange().end();++it)
    {
        KSpace::SCell indirectIncidentPixel = KImage.sCell(*it,pixelModel);
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
    }
}

void FlowGraphBuilder::createGluedCurveEdges(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                             SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                             std::map<Z2i::SCell,double>& weightMap)
{
    for(auto it=gluedRangeBegin;it!=gluedRangeEnd;++it){
        SegCut::SCellGluedCurveIterator begin = it->first;
        Z2i::SCell linel = begin.linkSurfel();

        createEdgeFromLinel(linel,
                            KImage,
                            weightMap);
    }

}

void FlowGraphBuilder::createEscapeEdges() {

    std::set<KSpace::SCell> forbiddenPoints;
    for(auto it=intCurve.begin();it!=intCurve.end();++it)
    {
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        forbiddenPoints.insert( innerPixel );
    }


    for(auto it=extCurve.begin();it!=extCurve.end();++it)
    {
        Dimension orthDir = KImage.sOrthDir(*it);
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);

        forbiddenPoints.insert( innerPixel );
    }

    Board2D board;
    board << extCurve;
    board << intCurve;
    board.saveEPS("ext.eps");

    std::set<KSpace::SCell> visitedNodes;
    int count=0;
    for(Curve::ConstIterator it = extCurve.begin();it!=extCurve.end();++it){
        Dimension orthDir = KImage.sOrthDir(*it);
        int x = KImage.sSign(*it)?1:-1;
        x = x*(1-2*orthDir);

        KSpace::SCell innerLinel = KImage.sGetAdd(*it,orthDir,x);
        KSpace::SCell candidate = KImage.sDirectIncident(innerLinel,orthDir);

        if( forbiddenPoints.find(candidate)!=forbiddenPoints.end() ) continue;

        ListDigraph::Arc in,out;
        KSpace::SCell innerPixel = KImage.sDirectIncident(*it,orthDir);
        createEdgeFromPixels(innerPixel,candidate,out);
        edgeWeight[out] = 0;

        createEdgeFromPixels(candidate,innerPixel,in);
        edgeWeight[in] = 0;

        explore(candidate,forbiddenPoints,visitedNodes);
        count++;
    }

    std::cout << count;

}


void FlowGraphBuilder::explore(KSpace::SCell candidate,
                               std::set<KSpace::SCell>& borderPoints,
                               std::set<KSpace::SCell>& visitedNodes)
{
    if(visitedNodes.find(candidate)!=visitedNodes.end()) return;

    visitedNodes.insert(candidate);
    SCells N = KImage.sProperNeighborhood(candidate);
    for(SCells::ConstIterator n=N.begin();n!=N.end();++n){
        if( borderPoints.find(*n)!=borderPoints.end() ) continue;

        explore(*n,borderPoints,visitedNodes);

        ListDigraph::Arc a;
        createEdgeFromPixels(candidate,*n,a);
        edgeWeight[a]=10;
    }
}
