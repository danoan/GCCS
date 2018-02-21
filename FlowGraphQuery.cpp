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

    FlowGraphBuilder::SubGraph sc(fgb.graph(),nodeFilter,arcFilter);


    ListDigraph::Node firstNode;
    {
        //First node must be an internal node.

        ListDigraph::Arc firstArc;
        sc.first(firstArc);

        while(fgb.arcType[firstArc]!=FlowGraphBuilder::ArcType::InternalCurveArc){
            sc.next(firstArc);
        }

        firstNode = sc.source(firstArc);
    }


    {
        ListDigraph::NodeMap<bool> visitedNodes(fgb.graph(), false);

        std::queue<ListDigraph::Arc> arcQueue;
        std::stack<ListDigraph::Node> traverseStack;

        ListDigraph::Node currentNode;

        traverseStack.push(firstNode);
        while(!traverseStack.empty()){
            currentNode = traverseStack.top(); traverseStack.pop();
            if(visitedNodes[currentNode]) continue;

            visitedNodes[currentNode] = true;

            for(ListDigraph::OutArcIt oai(fgb.graph(), currentNode);oai!=INVALID;++oai)
            {
                traverseStack.push(sc.target(oai));


                if(fgb.arcType[oai]==FlowGraphBuilder::ArcType::GluedArc &&
                   cutFilter[oai])
                {
                    arcQueue.push(oai);
                }
            }

        }

        /*
         * The graph needs to be traversed one node at time in order to guarantee the correct order.
         * The first gluedEdge added in arcQueue will necessarily be a extToInt edge.
         *
         * A glued edge pair points to different directions on the graph, leading to a shifting on
         * the arcQueue. The first element will actually match with the last one.
         */

        ListDigraph::Arc a = arcQueue.front();
        arcQueue.push(a);
        arcQueue.pop();

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
    for(ArcPairIterator ait = gluedEdgePairsBegin();ait!=gluedEdgePairsEnd();++ait)
    {
        ListDigraph::Arc intToExt = ait->first;
        ListDigraph::Arc extToInt = ait->second;



        ListDigraph::Node currentNode;
        ListDigraph::NodeMap<bool> visitedNodes(fgb.graph(),false);

        std::stack<ListDigraph::Node> traverseStack;
        traverseStack.push((fgb.graph().source(intToExt)));

        std::vector<ListDigraph::Arc>& workingVector = detourArcs[*ait];
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
                    workingVector.push_back(oai);
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



































































