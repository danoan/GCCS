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
    std::vector<Z2i::SCell> tempInt,tempExt;
    randomizeCurve(intCurve.begin(),intCurve.end(),intCurve.size(),tempInt);
    randomizeCurve(extCurve.begin(),extCurve.end(),extCurve.size(),tempExt);

    createGridCurveEdges(tempInt.begin(),tempInt.end(),sourceNode,false,weightMap);
    createGridCurveEdges(tempExt.begin(),tempExt.end(),targetNode,true,weightMap);


    SegCut::ConnectorSeedRangeType seedRange = getSeedRange(KImage,intCurve,extCurve);
    SegCut::SeedToGluedCurveRangeFunctor stgcF(gluedCurveLength);
    SegCut::GluedCurveSetRange gcsRange( seedRange.begin(),
                                 seedRange.end(),
                                 stgcF);

    createGluedCurveEdges(gcsRange.begin(),gcsRange.end(),weightMap);

    
    xavg=(int) xavg/(countNodes(fg));
    yavg=(int) yavg/(countNodes(fg));

    xavg = 160;
    yavg = 160;

    coords[sourceNode] = LemonPoint( xavg,yavg);
    coords[targetNode] = LemonPoint( xavg-40,yavg-40);
    
}

void FlowGraphBuilder::createEdgeFromLinel(Curve::SCell& linel,
                                           KSpace& KImage,
                                           std::map<Z2i::SCell,double>& weightMap,
                                           ListDigraph::Node& sourcePixelNode,
                                           ListDigraph::Node& targetPixelNode)
{
    Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
    Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

    Z2i::Point sourceCoord = directPixel.preCell().coordinates;
    Z2i::Point targetCoord = indirectPixel.preCell().coordinates;

    if(coordToNode.find(sourceCoord)==coordToNode.end())
    {
        coordToNode[sourceCoord] = fg.addNode();

        xavg+= sourceCoord[0];
        yavg+= sourceCoord[1];
    }

    if(coordToNode.find(targetCoord)==coordToNode.end())
    {
        coordToNode[targetCoord] = fg.addNode();
    }

    sourcePixelNode = coordToNode[sourceCoord];
    targetPixelNode = coordToNode[targetCoord];

    pixelMap[sourcePixelNode] = directPixel;
    pixelMap[targetPixelNode] = indirectPixel;

    coords[sourcePixelNode] = LemonPoint( sourceCoord[0],
                                          sourceCoord[1]);

    coords[targetPixelNode] = LemonPoint( targetCoord[0],
                                          targetCoord[1]);

    ListDigraph::Arc a = fg.addArc(sourcePixelNode,targetPixelNode);
    edgeWeight[a] = wFactor*weightMap[linel];
}

void FlowGraphBuilder::createGridCurveEdges(Curve::ConstIterator curveBegin,
                                            Curve::ConstIterator curveEnd,
                                            ListDigraph::Node terminalToConnect,
                                            bool toTerminal,
                                            std::map<Z2i::SCell,double>& weightMap)
{
    for(auto it=curveBegin;it!=curveEnd;++it){
        Curve::SCell linel = *it;

        ListDigraph::Node sourcePixelNode;
        ListDigraph::Node targetPixelNode;

        createEdgeFromLinel(linel,
                            KImage,
                            weightMap,
                            sourcePixelNode,
                            targetPixelNode);

        if(toTerminal){//Extern to target, for example
            ListDigraph::Arc terminalArc = fg.addArc(targetPixelNode,terminalToConnect);
            edgeWeight[terminalArc] = 10;
        }else{
            ListDigraph::Arc terminalArc = fg.addArc(terminalToConnect,sourcePixelNode);
            edgeWeight[terminalArc] = 10;
        }

    }
}

void FlowGraphBuilder::createGluedCurveEdges(SegCut::GluedCurveIteratorPair gluedRangeBegin,
                                             SegCut::GluedCurveIteratorPair gluedRangeEnd,
                                             std::map<Z2i::SCell,double>& weightMap)
{
    for(auto it=gluedRangeBegin;it!=gluedRangeEnd;++it){
        SegCut::SCellGluedCurveIterator begin = it->first;
        Z2i::SCell linel = begin.linkSurfel();

        ListDigraph::Node sourcePixelNode;
        ListDigraph::Node targetPixelNode;

        createEdgeFromLinel(linel,
                            KImage,
                            weightMap,
                            sourcePixelNode,
                            targetPixelNode);

    }

}