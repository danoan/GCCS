#include "FlowGraphQuery.h"

void FlowGraphQuery::sourceComponent(ListDigraph::ArcMap<bool>& cutFilter)
{
    ListDigraph::NodeMap<bool> allNodes(fgb.graph(),true);
    ListDigraph::ArcMap<bool> allArcs(fgb.graph(),true);

    sourceComponent(allNodes,
                    allArcs,
                    cutFilter);
}

void FlowGraphQuery::sourceComponent(ListDigraph::NodeMap<bool>& nodeFilter,
                                     ListDigraph::ArcMap<bool>& arcFilter,
                                     ListDigraph::ArcMap<bool>& cutFilter)
{
    FlowGraphBuilder::FlowComputer flow = fgb.preparePreFlow();
    flow.run();

    for(ListDigraph::ArcIt a(fgb.graph());a!=INVALID;++a)
    {
        if( flow.minCut( fgb.graph().source(a) ) )
        {
            nodeFilter[ fgb.graph().target(a) ] = true;
            nodeFilter[ fgb.graph().source(a) ] = true;
            arcFilter[a] = true;

            if(!flow.minCut( fgb.graph().target(a) ))
            {
                cutFilter[a] = true;
            }
        }
    }

    nodeFilter[fgb.source()]=false;
    nodeFilter[fgb.target()]=false;
}


void FlowGraphQuery::setGluedEdgePairsOnCut()
{
    //We first compute the sourceComponent of the graph
    ListDigraph::NodeMap<bool> nodeFilter(fgb.graph(),false);
    ListDigraph::ArcMap<bool> arcFilter(fgb.graph(),false);
    ListDigraph::ArcMap<bool> cutFilter(fgb.graph(),false);

    sourceComponent(nodeFilter,
                    arcFilter,
                    cutFilter);

    ListDigraph::ArcMap<bool> intExtGluedArcsFilter (fgb.graph(),false);
    filterArcs(intExtGluedArcsFilter,FlowGraphBuilder::ArcType::IntExtGluedArc,true);

    ListDigraph::ArcMap<bool> refundArcsFilter (fgb.graph(),false);
    filterArcs(refundArcsFilter,FlowGraphBuilder::ArcType::RefundArc,true);

    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> arcFilterComposer(fgb.graph());
    (arcFilterComposer+arcFilter)-intExtGluedArcsFilter+cutFilter-refundArcsFilter;

    FlowGraphBuilder::SubGraph sc(fgb.graph(),nodeFilter,arcFilterComposer.initialMap);

    graphToEps(sc,"vamola.eps")
            .coords(fgb.coordsMap())
            .nodeScale(.0005)
            .arcWidthScale(0.0005)
            .drawArrows()
            .arrowWidth(0.25)
            .run();


    std::stack<ListDigraph::Arc> intExtGluedArcs;
    for(FlowGraphBuilder::SubGraph::ArcIt a(sc);a!=INVALID;++a)
    {
        if(fgb.arcType[a]==FlowGraphBuilder::ArcType::IntExtGluedArc &&
           cutFilter[a])
        {
            intExtGluedArcs.push(a);
        }
    }


    std::queue<ListDigraph::Arc> arcQueue;
    std::stack<ListDigraph::Node> traverseStack;
    ListDigraph::NodeMap<bool> visitedNodes(fgb.graph(), false);
    ListDigraph::Node currentNode;
    while(!intExtGluedArcs.empty())
    {
        ListDigraph::Arc a = intExtGluedArcs.top(); intExtGluedArcs.pop();
        arcQueue.push(a);
        ListDigraph::Node firstNode = sc.source(a);

        traverseStack.push(firstNode);
        while(!traverseStack.empty()){
            currentNode = traverseStack.top(); traverseStack.pop();
            if(visitedNodes[currentNode]) continue;

            visitedNodes[currentNode] = true;

            for(FlowGraphBuilder::SubGraph::OutArcIt oai(sc, currentNode);oai!=INVALID;++oai)
            {
                traverseStack.push(sc.target(oai));


                if( fgb.arcType[oai]==FlowGraphBuilder::ArcType::ExtIntGluedArc &&
                    cutFilter[oai])
                {
                    arcQueue.push(oai);
                }
            }

        }

    }


    {
        if(arcQueue.size()%2==1) throw "arcQueue has an odd number os elements";

        while(!arcQueue.empty())
        {
            ListDigraph::Arc intToExt = arcQueue.front();
            arcQueue.pop();
            ListDigraph::Arc extToInt = arcQueue.front();
            arcQueue.pop();

            gluedEdgePairsOnCut.push_back(ArcPair(intToExt, extToInt));
        }
    }

}

FlowGraphQuery::ArcPairIterator FlowGraphQuery::gluedEdgePairsBegin()
{
    gluedEdgePairsOnCut.clear();
    setGluedEdgePairsOnCut();
    return gluedEdgePairsOnCut.begin();
}

FlowGraphQuery::ArcPairIterator FlowGraphQuery::gluedEdgePairsEnd()
{
    return gluedEdgePairsOnCut.end();
}



void FlowGraphQuery::setDetourArcs()
{
    SCellToArc staFunctor(fgb.scellArc);
    int length = 2;
    for(ArcPairIterator ait = gluedEdgePairsBegin();ait!=gluedEdgePairsEnd();++ait)
    {
        ListDigraph::Arc intToExt = ait->first;
        ListDigraph::Arc extToInt = ait->second;

        FlowGraphBuilder::CirculatorPair cIn = fgb.arcCirculator[intToExt];
        FlowGraphBuilder::CirculatorPair cEx = fgb.arcCirculator[extToInt];


        std::set<ListDigraph::Arc>& workingSet = detourArcs[*ait];


        FlowGraphBuilder::SCellCirculator beginRightExternal = cIn.second;
        insertSCellFromArc(workingSet,staFunctor,beginRightExternal,length);

        FlowGraphBuilder::SCellCirculator endLeftExternal(cEx.first);
        ++endLeftExternal;
        FlowGraphBuilder::SCellCirculator beginLeftExternal = walkCirculator(endLeftExternal,-length);

        insertSCellFromArc(workingSet,staFunctor,beginLeftExternal,length);


        FlowGraphBuilder::SCellCirculator beginLeftInternal(cEx.second);
        insertSCellFromArc(workingSet,staFunctor,beginLeftInternal,length);


        FlowGraphBuilder::SCellCirculator endRightInternal(cIn.first);
        ++endRightInternal;
        FlowGraphBuilder::SCellCirculator beginRightInternal = walkCirculator(endRightInternal,-length);

        insertSCellFromArc(workingSet,staFunctor,beginRightInternal,length);


        if(fgb.arcType[fgb.scellArc[*beginRightExternal]]!=FlowGraphBuilder::ArcType::ExternalCurveArc) throw std::runtime_error("External Arc Expected");
        if(fgb.arcType[fgb.scellArc[*beginLeftExternal]]!=FlowGraphBuilder::ArcType::ExternalCurveArc) throw std::runtime_error("External Arc Expected");
        if(fgb.arcType[fgb.scellArc[*beginRightInternal]]!=FlowGraphBuilder::ArcType::InternalCurveArc) throw std::runtime_error("Internal Arc Expected");
        if(fgb.arcType[fgb.scellArc[*beginLeftInternal]]!=FlowGraphBuilder::ArcType::InternalCurveArc) throw std::runtime_error("Internal Arc");

    }
}

//Deprecated
void FlowGraphQuery::setDetourArcsAsExternalArcs()
{
    for(ArcPairIterator ait = gluedEdgePairsBegin();ait!=gluedEdgePairsEnd();++ait)
    {
        ListDigraph::Arc intToExt = ait->first;
        ListDigraph::Arc extToInt = ait->second;



        ListDigraph::Node currentNode;
        ListDigraph::NodeMap<bool> visitedNodes(fgb.graph(),false);

        std::stack<ListDigraph::Node> traverseStack;
        traverseStack.push((fgb.graph().source(intToExt)));

        std::set<ListDigraph::Arc>& workingSet = detourArcs[*ait];
        while(!traverseStack.empty())
        {
            currentNode = traverseStack.top(); traverseStack.pop();
            if(visitedNodes[currentNode]) continue;
            visitedNodes[currentNode]=true;

            for(ListDigraph::OutArcIt oai(fgb.graph(),currentNode);oai!=INVALID;++oai)
            {
                if(oai==intToExt) continue;
                if(oai==extToInt) continue;

                traverseStack.push( fgb.graph().target(oai) );
                if(fgb.arcType[oai]==FlowGraphBuilder::ArcType::ExternalCurveArc)
                {
                    workingSet.insert(oai);
                }
            }
        }

    }
}

FlowGraphQuery::DetourArcMapIterator FlowGraphQuery::detourArcsBegin()
{
    detourArcs.clear();
    setDetourArcs();
    return detourArcs.begin();
}


FlowGraphQuery::DetourArcMapIterator FlowGraphQuery::detourArcsEnd()
{
    return detourArcs.end();
}



































































