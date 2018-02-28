#include "FlowGraphQuery.h"

void FlowGraphQuery::sourceComponent(FlowGraph& fg,
                                     ListDigraph::ArcMap<bool>& cutFilter)
{
    ListDigraph::NodeMap<bool> allNodes(fg.graph(),true);
    ListDigraph::ArcMap<bool> allArcs(fg.graph(),true);

    FlowGraphQuery::sourceComponent(fg,
                                    allNodes,
                                    allArcs,
                                    cutFilter);
}

void FlowGraphQuery::sourceComponent(FlowGraph& fg,
                                     ListDigraph::NodeMap<bool>& nodeFilter,
                                     ListDigraph::ArcMap<bool>& arcFilter,
                                     ListDigraph::ArcMap<bool>& cutFilter)
{
    FlowGraph::FlowComputer flow = fg.prepareFlow();
    flow.run();

    for(ListDigraph::ArcIt a(fg.graph());a!=lemon::INVALID;++a)
    {
        if( flow.minCut( fg.graph().source(a) ) )
        {
            nodeFilter[ fg.graph().target(a) ] = true;
            nodeFilter[ fg.graph().source(a) ] = true;
            arcFilter[a] = true;

            if(!flow.minCut( fg.graph().target(a) ))
            {
                cutFilter[a] = true;
            }
        }
    }

    nodeFilter[fg.source()]=false;
    nodeFilter[fg.target()]=false;
}


void FlowGraphQuery::selectArcsInCut(FlowGraph& fg,
                                     ListDigraph::ArcMap<bool>& selectedArcs,
                                     FlowGraph::ArcType at)
{
    ListDigraph::ArcMap<bool> cutFilter(fg.graph(),false);
    sourceComponent(fg,
                    cutFilter);

    ListDigraph::ArcMap<bool> arcOfType (fg.graph(),false);
    filterArcs(fg,
               arcOfType,
               at,
               true);

    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> filterComposer(fg.graph(),selectedArcs);
    (filterComposer+cutFilter)*arcOfType;

}

void FlowGraphQuery::coreArcs(FlowGraph& fg,
                              ListDigraph::ArcMap<bool>& coreArcs)
{
    ListDigraph::ArcMap<bool> sourceArcs(fg.graph());
    filterArcs(fg,
               sourceArcs,
               FlowGraph::ArcType::SourceArc,
               true);

    ListDigraph::ArcMap<bool> targetArcs(fg.graph());
    filterArcs(fg,
               sourceArcs,
               FlowGraph::ArcType::TargetArc,
               true);

    FilterComposer<ListDigraph::ArcMap<bool>,ListDigraph::ArcIt> filterComposer(fg.graph(),
                                                                                coreArcs);

    ListDigraph::ArcMap<bool> allArcs(fg.graph(),true);
    (filterComposer+allArcs)-sourceArcs-targetArcs;
}


void FlowGraphQuery::gluedArcPairSet(FlowGraph& fg,
                                        std::set<ArcPair>& gluedArcPairSet)
{
    ListDigraph::NodeMap<bool> allNodes(fg.graph(),true);
    ListDigraph::ArcMap<bool> intExtGluedArcsFilter(fg.graph(),false);
    selectArcsInCut(fg,
                    intExtGluedArcsFilter,
                    FlowGraph::ArcType::IntExtGluedArc);


    SubGraph intExtGluedSubgraph(fg.graph(),
                                 allNodes,
                                 intExtGluedArcsFilter);

    lemon::graphToEps(intExtGluedSubgraph,"myTest.eps")
            .coords(fg.coordsMap())
            .nodeScale(.0005)
            .arcWidthScale(0.0005)
            .drawArrows()
            .arrowWidth(0.25)
            .run();


    ListDigraph::ArcMap<bool> noSourceTargetArcs(fg.graph(),false);
    coreArcs(fg,
             noSourceTargetArcs);

    SubGraph coreSubgraph(fg.graph(),allNodes,noSourceTargetArcs);


    ListDigraph::ArcMap<bool> cutFilter(fg.graph(),false);
    sourceComponent(fg,
                    cutFilter);


    std::queue<ListDigraph::Arc> arcQueue;
    std::stack<ListDigraph::Node> traverseStack;
    ListDigraph::NodeMap<bool> visitedNodes(fg.graph(), false);
    ListDigraph::Node currentNode;
    bool pairNotFound;
    for(SubGraph::ArcIt a(intExtGluedSubgraph);a!=lemon::INVALID;++a)
    {
        arcQueue.push(a);
        ListDigraph::Node firstNode = intExtGluedSubgraph.source(a);
        visitedNodes[intExtGluedSubgraph.target(a)]=true;

        traverseStack.push(firstNode);
        pairNotFound=true;
        while(pairNotFound){
            currentNode = traverseStack.top(); traverseStack.pop();
            if(visitedNodes[currentNode]) continue;

            visitedNodes[currentNode] = true;

            for(SubGraph::OutArcIt oai(coreSubgraph, currentNode);oai!=lemon::INVALID;++oai)
            {
                traverseStack.push(coreSubgraph.target(oai));


                if( fg.arcType(oai)==FlowGraph::ArcType::ExtIntGluedArc &&
                    cutFilter[oai])
                {
                    arcQueue.push(oai);
                    pairNotFound=false;
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

            gluedArcPairSet.insert(ArcPair(intToExt, extToInt));
        }
    }

}


void FlowGraphQuery::detourArcMap(FlowGraph& fg,
                                  DetourArcMap& detourArcMap)
{
    SCellToArc staFunctor(fg);
    int length = fg.getDiffDistance();
    
    ArcPairSet gaps;
    gluedArcPairSet(fg,gaps);
    
    for(ArcPairIterator ait = gaps.begin();ait!=gaps.end();++ait)
    {
        ListDigraph::Arc intToExt = ait->first;
        ListDigraph::Arc extToInt = ait->second;

        FlowGraph::CirculatorPair cIn = fg.circulatorPair(intToExt);
        FlowGraph::CirculatorPair cEx = fg.circulatorPair(extToInt);


        std::set<ListDigraph::Arc>& workingSet = detourArcMap[*ait];


        FlowGraph::SCellCirculator beginRightExternal = cIn.second;
        insertSCellFromArc(fg,workingSet,staFunctor,beginRightExternal,length);

        FlowGraph::SCellCirculator endLeftExternal(cEx.first);
        ++endLeftExternal;
        FlowGraph::SCellCirculator beginLeftExternal = moveIterator(endLeftExternal,-length);

        insertSCellFromArc(fg,workingSet,staFunctor,beginLeftExternal,length);


        FlowGraph::SCellCirculator beginLeftInternal(cEx.second);
        insertSCellFromArc(fg,workingSet,staFunctor,beginLeftInternal,length);


        FlowGraph::SCellCirculator endRightInternal(cIn.first);
        ++endRightInternal;
        FlowGraph::SCellCirculator beginRightInternal = moveIterator(endRightInternal,-length);

        insertSCellFromArc(fg,workingSet,staFunctor,beginRightInternal,length);


        if(fg.arcType(fg.arc(*beginRightExternal))!=FlowGraph::ArcType::ExternalCurveArc) throw std::runtime_error("External Arc Expected");
        if(fg.arcType(fg.arc(*beginLeftExternal))!=FlowGraph::ArcType::ExternalCurveArc) throw std::runtime_error("External Arc Expected");
        if(fg.arcType(fg.arc(*beginRightInternal))!=FlowGraph::ArcType::InternalCurveArc) throw std::runtime_error("Internal Arc Expected");
        if(fg.arcType(fg.arc(*beginLeftInternal))!=FlowGraph::ArcType::InternalCurveArc) throw std::runtime_error("Internal Arc");

    }
}

void FlowGraphQuery::detourArcFilter(FlowGraph& fg,
                                     ListDigraph::ArcMap<bool>& arcFilter,
                                     DetourArcMap& detourArcMap)
{
    for(DetourArcMapIterator dami=detourArcMap.begin();dami!=detourArcMap.end();++dami)
    {
        DetourArcIterator begin = dami->second.begin();
        DetourArcIterator end = dami->second.end();
        for(DetourArcIterator dai=begin;dai!=end;++dai )
        {
            arcFilter[*dai] = true;
        }
    }
}


int FlowGraphQuery::arcDistance(FlowGraph& fg,
                                ListDigraph::Arc& a1,
                                ListDigraph::Arc& a2)
{
    //a1 -> extInt
    FlowGraph::SCellCirculator c1 = fg.circulatorPair(a1).first; //End of external segment

    //a1 -> intExt
    FlowGraph::SCellCirculator c2 = fg.circulatorPair(a2).second; //Begin of external segment

    int n =0;
    while(n<fg.getConsecutiveGluedPairDistance() && c1!=c2)
    {
        ++c1;
        ++n;
    }

    return n;
}

void FlowGraphQuery::pixelsFilter(FlowGraph& fg,
                                  ListDigraph::NodeMap<bool>& pixelsFilter)
{
    ListDigraph::NodeMap<bool> nodeFilter(fg.graph(),false);
    ListDigraph::ArcMap<bool> arcFilter(fg.graph(),false);
    ListDigraph::ArcMap<bool> cutFilter(fg.graph(),false);

    sourceComponent(fg,
                    nodeFilter,
                    arcFilter,
                    cutFilter);


    SubGraph subgraph(fg.graph(),
                      nodeFilter,
                      arcFilter);

    for(SubGraph::ArcIt a(subgraph);a!=lemon::INVALID;++a)
    {
        pixelsFilter[ subgraph.source(a) ] = true;
    }
    pixelsFilter[fg.source()] = false;
    pixelsFilter[fg.target()] = false;
}



double FlowGraphQuery::computeEnergyValue(FlowGraph& fg,
                                          ListDigraph::ArcMap<bool>& arcFilter)
{
    double s = 0;
    for(ListDigraph::ArcIt a(fg.graph());a!=lemon::INVALID;++a)
    {
        if(arcFilter[a]){
            s+= fg.weight(a);
        }
    }

    return s;

}

FlowGraphQuery::ArcKey FlowGraphQuery::arcKey(FlowGraph& fg,
                                              ListDigraph::Arc& arc)
{
    SCell source = fg.pixel( fg.source(arc) );
    SCell target = fg.pixel( fg.target(arc) );

    return ArcKey( source.preCell().coordinates,
                   target.preCell().coordinates );
}


void FlowGraphQuery::arcFilterConversion(ListDigraph::ArcMap<bool>& arcMap,
                                         FlowGraph& fg,
                                         ListDigraph::ArcMap<bool>& newArcMap,
                                         FlowGraph& newFg)
{
    std::set<ArcKey> keys;
    for(ListDigraph::ArcIt a(fg.graph());a!=lemon::INVALID;++a)
    {
        if(arcMap[a])
        {
            keys.insert( arcKey(fg,a) );
        }

    }

    ListDigraph::NodeMap<bool> allNodes(newFg.graph(),true);
    ListDigraph::ArcMap<bool> internalArcs(newFg.graph(),false);

    FlowGraphQuery::filterArcs(newFg,internalArcs,
                               FlowGraph::ArcType::InternalCurveArc,
                               true);

    SubGraph subgraph(newFg.graph(),
                      allNodes,
                      internalArcs);

    for(SubGraph::ArcIt a(subgraph);a!=lemon::INVALID;++a)
    {
        if(keys.find(arcKey(newFg,a))!=keys.end())
        {
            newArcMap[a] = true;
        }
    }
}

double FlowGraphQuery::cutValue(FlowGraph& fg)
{
    FlowGraph::FlowComputer flowComputer = fg.prepareFlow();
    flowComputer.run();

    return flowComputer.flowValue();
}