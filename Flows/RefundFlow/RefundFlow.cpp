#include "RefundFlow.h"
#include "../../FlowGraph/FlowGraphDebug.h"

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
        //Only internalArcs in F2 corresponds to the energy comparison
        if(f2.arcType(arcF2)!=FlowGraph::InternalCurveArc) continue;

        FlowGraphQuery::ArcKey arcF2Key = FlowGraphQuery::arcKey(f2, arcF2);

        if (unifiedArcMapping.find(arcF2Key) != unifiedArcMapping.end()) {
            ListDigraph::Arc arcF1 = unifiedArcMapping[arcF2Key];
            mapping[arcF1] = f2.graph().id(arcF2);
        }
    }

}

RefundFlow::RefundFlow(ImageFlowData& imageFlowData,int maxIterations):imageFlowData(imageFlowData),
                                                                       imageOut(imageFlowData.getOriginalImage()),
                                                                       maxIterations(maxIterations){};





void RefundFlow::updateImage(FlowGraph& fg,
                             Image2D& updatedImage)
{

    ListDigraph::NodeMap<bool> pixelsFilter(fg.graph(),false);
    FlowGraphQuery::sourceComponentNodes(fg,
                                         pixelsFilter);


    KSpace& KImage = imageFlowData.getKSpace();
    std::vector<Z2i::SCell> pixelsInTheGraph;
    pixels(fg,
           pixelsInTheGraph,
           pixelsFilter);

    updatedImage = imageFlowData.getMostInnerImage();
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
                                            ListDigraph::ArcMap<bool>& arcFilter,
                                            ListDigraph::ArcMap<double>& arcDiff)
{
    Image2D image = imageFlowData.getOriginalImage();
    updateImage(f1,image);

    GenericWriter<Image2D>::exportFile("test.pgm",image);

    ImageFlowData tempImageFlowData(image);
    tempImageFlowData.init(f1.getFlowMode(),
                       f1.getGluedCurveLength());

    FlowGraphBuilder::LinelWeightMap weightMap;
    setArcsWeight(tempImageFlowData,weightMap);

    FlowGraph f2;
    FlowGraphBuilder fgb(f2,tempImageFlowData,weightMap);

    ListDigraph::ArcMap<int> arcsMatch(f1.graph(),-1);
    matchFlows(arcsMatch,f1,f2);

    ListDigraph::ArcMap<bool> printArcs(f1.graph(),false);

    for(ListDigraph::ArcIt arcF1(f1.graph());arcF1!=INVALID;++arcF1)
    {
        if(arcFilter[arcF1] && arcsMatch[arcF1]!=-1)
        {
            printArcs[arcF1]=true;
            ListDigraph::Arc arcF2 = f2.graph().arcFromId( arcsMatch[arcF1] );
            arcDiff[arcF1] = f2.weight(arcF2) - f1.weight(arcF1);
        }
    }

    FlowGraphDebug fgd(f1);
    fgd.highlightArcs(printArcs,"cmake-build-debug","highlight");

}

int RefundFlow::findConflicts(FlowGraph& fg,
                              std::set<FlowGraphQuery::ArcPair>& listOfPairs,
                              FlowGraphQuery::ArcPair& pairToTest,
                              FlowGraphQuery::ArcPair conflictedPair[])
{
    ListDigraph::Arc testIntExt = pairToTest.first;
    ListDigraph::Arc testExtInt = pairToTest.second;

    int n =0;
    int conflictDistance = fg.getConsecutiveGluedPairDistance();
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

void RefundFlow::checkCodeConsistence(FlowGraph& fg,Image2D& resultingImage)
{
    std::set<ListDigraph::Arc> globalDetourArcSet;
    FlowGraphQuery::globalDetourArcSet(fg,globalDetourArcSet);

    ListDigraph::ArcMap<bool> intExtCut(fg.graph(),false);
    FlowGraphQuery::selectArcsInCut(fg,intExtCut,FlowGraph::ArcType::IntExtGluedArc);

    ListDigraph::ArcMap<bool> extIntCut(fg.graph(),false);
    FlowGraphQuery::selectArcsInCut(fg,extIntCut,FlowGraph::ArcType::ExtIntGluedArc);

    ListDigraph::ArcMap<bool> globalDetourArcFilter(fg.graph(),false);
    for(auto it=globalDetourArcSet.begin();it!=globalDetourArcSet.end();++it)
    {
        globalDetourArcFilter[*it] = true;
    }

    ImageFlowData imf(resultingImage);
    imf.init(imageFlowData.getFlowMode(),imageFlowData.getGluedCurveLength());

    FlowGraphBuilder::LinelWeightMap tempWeightMap;
    setArcsWeight(imf,tempWeightMap);

    FlowGraph tempFG;
    FlowGraphBuilder fgb(tempFG,imf,tempWeightMap);

    ListDigraph::ArcMap<bool> matchFilter(fg.graph(),false);
    ListDigraph::ArcMap<int> matches(fg.graph(),-1);
    matchFlows(matches,fg,tempFG);

    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(matches[a]!=-1)
        {
            matchFilter[a] = true;
        }
    }

    FilterComposer<ListDigraph::ArcMap<bool>, ListDigraph::ArcIt> validDetourArcs(fg.graph(), globalDetourArcFilter);
    validDetourArcs*matchFilter;
    validDetourArcs+intExtCut+extIntCut;


    double cutDetourArcsSum = FlowGraphQuery::computeEnergyValue(fg,validDetourArcs.initialMap);


    ListDigraph::ArcMap<bool> tempArcFilter(tempFG.graph(),false);
    FlowGraphQuery::arcFilterConversion(validDetourArcs.initialMap,fg,tempArcFilter,tempFG);

    double solutionDetourArcsSum = FlowGraphQuery::computeEnergyValue(tempFG,tempArcFilter);

    double refundArcsSum = solutionDetourArcsSum - cutDetourArcsSum;

    double cutValue = FlowGraphQuery::cutValue(fg);
    double resultingImageEnergyValue = energyValue(resultingImage);

    std::string testResult = abs( (cutValue + refundArcsSum)-resultingImageEnergyValue ) < 0.0001?"TRUE":"FALSE";

    std::cout << "DIFF CHECK::" << cutDetourArcsSum << "::" << solutionDetourArcsSum <<  "::" << refundArcsSum <<std::endl;
    std::cout << testResult << "::" << cutValue+refundArcsSum << std::endl;

}

int RefundFlow::runIteration(FlowGraph& fg,
                             Image2D& partialImage,
                             FlowGraphQuery::ArcPairSet& usedKeys,
                             int iteration)
{
    FlowGraphQuery::DetourArcMap detourArcMap;
    FlowGraphQuery::detourArcMap(fg,detourArcMap);

    ListDigraph::ArcMap<bool> detourArcs(fg.graph(),false);
    FlowGraphQuery::detourArcFilter(fg,detourArcs,detourArcMap);

    ListDigraph::ArcMap<double> arcDiff(fg.graph(),0);


    computeCutEnergyDifference(fg,
                               detourArcs,
                               arcDiff);


    FlowGraphQuery::ArcPair conflictedPairs[2];
    double s;
    int n;
    double totalDiff=0;
    for(FlowGraphQuery::DetourArcMapIterator dami = detourArcMap.begin();dami!=detourArcMap.end();++dami)
    {
        FlowGraphQuery::ArcPair key = dami->first;
        if(usedKeys.find(key)!=usedKeys.end()) continue;


//        int conflicts = findConflicts(fg,
//                                      usedKeys,
//                                      key,
//                                      conflictedPairs);
//        for(int i=0;i<conflicts;++i)
//        {
//            addConflictArcs(fg,
//                            conflictedPairs[i].first,
//                            conflictedPairs[i].second,
//                            10);
//        }

        usedKeys.insert(key);

        const std::set<ListDigraph::Arc> &detourArcs = dami->second;

        s=0;
        n=0;
        for(auto vit=detourArcs.begin();vit!=detourArcs.end();++vit)
        {
            s+= arcDiff[*vit];
            ++n;
        }
        totalDiff+=s;
        s/=2;

        std::cout << "ArcDiff:" << s << " in " << n << " arcs" << std::endl;

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

double RefundFlow::energyValue(Image2D& image)
{
    ImageFlowData tempImageFlowData(image);
    tempImageFlowData.init(imageFlowData.getFlowMode(),imageFlowData.getGluedCurveLength());

    FlowGraphBuilder::LinelWeightMap weightMap;
    setArcsWeight(tempImageFlowData,weightMap);

    FlowGraph fg;
    FlowGraphBuilder fgb(fg,tempImageFlowData,weightMap);

    ListDigraph::ArcMap<bool> internalArcs(fg.graph(),false);
    FlowGraphQuery::filterArcs(fg,internalArcs,FlowGraph::ArcType::InternalCurveArc,true);
    return FlowGraphQuery::computeEnergyValue(fg,internalArcs);
}

void RefundFlow::constructSCellsVectorFromFilter(std::vector<SCell>& scells,
                                                 FlowGraph& fg,
                                                 ListDigraph::ArcMap<bool>& arcFilter)
{
    for(ListDigraph::ArcIt a(fg.graph());a!=INVALID;++a)
    {
        if(arcFilter[a])
        {
            if(fg.arcType(a)==FlowGraph::IntExtGluedArc ||
               fg.arcType(a)==FlowGraph::ExtIntGluedArc ||
               fg.arcType(a)==FlowGraph::InternalCurveArc ||
               fg.arcType(a)==FlowGraph::ExternalCurveArc )
            {
                scells.push_back(fg.scell(a));
            }
        }
    }
}

void RefundFlow::debugData(FlowGraph& fg,
                           Image2D& partialImage,
                           FlowGraphBuilder::LinelWeightMap& weightMap,
                           int iteration)
{
    double currentCutValue = FlowGraphQuery::cutValue(fg);
    double currentEnergyValue = energyValue(partialImage);
    std::cout << "OK-Sub-" << iteration
              << "    Cut Value:" << currentCutValue
              << "    Energy Value:" << currentEnergyValue << std::endl;


    ImageFlowData imf(partialImage);
    imf.init(imageFlowData.getFlowMode(),imageFlowData.getGluedCurveLength());

    FlowGraphBuilder::LinelWeightMap tempWeightMap;
    setArcsWeight(imf,tempWeightMap);

    FlowGraph tempFG;
    FlowGraphBuilder fgb(tempFG,imf,tempWeightMap);


    ListDigraph::ArcMap<bool> imageArcs(tempFG.graph(),false);
    FlowGraphQuery::filterArcs(tempFG,imageArcs,FlowGraph::ArcType::InternalCurveArc,true);

    std::vector<SCell> imageSCells;
    constructSCellsVectorFromFilter(imageSCells,tempFG,imageArcs);

    Board2D board;
    double cmin=100,cmax=-100;

    std::function<double(SCell)> fnTempToDouble = [tempWeightMap](SCell s){ return tempWeightMap.at(s); };
    max_and_min(imageSCells,cmin,cmax,fnTempToDouble);
    cmin = cmin==cmax?cmax-0.000001:cmin;


    std::string suffix = std::to_string(iteration) + "-" + std::to_string(iteration);
    draw(tempWeightMap,imageSCells,board,cmin,cmax);
    board.saveEPS( ("../output/refundFlow/square/square-5/" + std::string("energy-weights-") + suffix).c_str() );



    ListDigraph::ArcMap<bool> cutResultArcs(fg.graph(),false);
    FlowGraphQuery::arcFilterConversion(imageArcs,tempFG,cutResultArcs,fg);

    std::vector<SCell> cutResultSCells;
    constructSCellsVectorFromFilter(cutResultSCells,fg,cutResultArcs);

    std::function<double(SCell)> fnToDouble = [weightMap](SCell s){ return weightMap.at(s); };
    max_and_min(cutResultSCells,cmin,cmax,fnToDouble);
    cmin = cmin==cmax?cmax-0.000001:cmin;

    draw(weightMap,cutResultSCells,board,cmin,cmax);
    board.saveEPS( ("../output/refundFlow/square/square-5/" + std::string("cutResult-weights-") + suffix).c_str() );
}

double RefundFlow::run(int mainIteration)
{
    FlowGraph fg;
    LinelWeightMap weightMap;
    setArcsWeight(imageFlowData,weightMap);
    FlowGraphBuilder fgb(fg,imageFlowData,weightMap);


    int iteration = 1;


    FlowGraphQuery::ArcPairSet usedKeys;
    Image2D partialImage = imageFlowData.getOriginalImage();
    Image2D previousImage=partialImage;

    double currentEnergyValue;
    double initialEnergyValue = energyValue(partialImage);
    while(runIteration(fg,partialImage,usedKeys,iteration)==Continue && iteration <= maxIterations)
    {
        FlowGraphDebug fgd(fg);

        previousImage = partialImage;
        updateImage(fg,partialImage);

        currentEnergyValue = energyValue(partialImage);

//        checkCodeConsistence(fg,partialImage);
//        debugData(fg,partialImage,weightMap,iteration);

        if(currentEnergyValue<initialEnergyValue)
        {
            initialEnergyValue = currentEnergyValue;
        }

        if(previousImage==partialImage) break;


        usedKeys.clear();
        ++iteration;
    }

    imageOut = partialImage;

    return energyValue(imageOut);
}
