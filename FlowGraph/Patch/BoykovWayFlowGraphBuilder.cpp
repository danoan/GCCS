#include "BoykovWayFlowGraphBuilder.h"

namespace Development {
    FlowGraphBuilder::FlowGraphBuilder(FlowGraph& fg,
                                       ImageFlowData& imageFlowData,
                                       LinelWeightMap& weightMap,
                                       double* distrFrg,
                                       double* distrBkg):imageFlowData(imageFlowData),
                                                                  weightMap(weightMap),
                                                                  fg(fg)
    {
        fg.gluedCurveLength = imageFlowData.getGluedCurveLength();
        fg.consecutiveGluedPairDistance = imageFlowData.getConsecutiveGluedPairDistance();
        fg.diffDistance = imageFlowData.getDiffDistance();
        fg.flowGraphMode = imageFlowData.getFlowMode();



        std::map<SCell,bool> superposedLinels;
        for(auto it=imageFlowData.curveDataBegin();it!=imageFlowData.curveDataEnd();++it)
        {
            identifySuperposedLinels(it->curve.begin(),it->curve.end(),superposedLinels);
        }

        createDilatedPixelsTerminalArcs(imageFlowData.getMostInnerCurve(),
                                        imageFlowData.getMostOuterCurve(),
                                        weightMap,
                                        distrFrg,
                                        distrBkg);

        createOriginalPixelsTerminalArc(imageFlowData.getMostInnerCurve(),
                                        weightMap,
                                        distrFrg,
                                        distrBkg);


        UnsignedSCellComparison myComp;
        UnsignedSCellSet visitedNodes(myComp);

        for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
        {
//            createGluedArcs(it->gcsRangeBegin(),
//                            it->gcsRangeEnd(),
//                            weightMap,
//                            visitedNodes);

            createEscapeArcs(it->intCurveData.curve,
                             it->extCurveData.curve,
                             visitedNodes,
                             superposedLinels
            );
        }


        for(auto it=imageFlowData.curvePairBegin();it!=imageFlowData.curvePairEnd();++it)
        {

            /*createTargetArcsFromGluedSegments(it->gcsRangeBegin(),
                                              it->gcsRangeEnd(),
                                              weightMap,
                                              visitedNodes);
                                              */
        }

        similarityBetweenDilated(imageFlowData.getMostOuterCurve());
        similarityBetweenCurves(imageFlowData.getMostInnerCurve());


        setTerminalsCoordinates();

    }


    void FlowGraphBuilder::createSourceArcs(Curve& erodedCurve,
                                            UnsignedSCellSet& visitedNodes,
                                            double* distrFrg,
                                            double* distrBkg)
    {
        KSpace& KImage = imageFlowData.getKSpace();

        KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);
        for(auto it=erodedCurve.begin();it!=erodedCurve.end();++it)
        {
            Z2i::SCell linel = *it;

            Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));
            Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));

            if(visitedNodes.find(directPixel)!=visitedNodes.end()) continue;
            visitedNodes.insert(directPixel);

            int iv = imageFlowData.getRefImage().operator()(KImage.sCoords(directPixel));


            //Source
            {
                ListDigraph::Arc aSource;
                connectNodeToPixel(fg.sourceNode, directPixel, aSource);
                fg.arcWeightMap[aSource] = inverseProp(iv, distrBkg) + 1.0 / weightMap[linel];
                fg.arcTypeMap[aSource] = FlowGraph::ArcType::SourceArc;
            }

            //Target
            {
                ListDigraph::Arc aTarget;
                connectPixelToNode(directPixel,fg.targetNode, aTarget);
                fg.arcWeightMap[aTarget] = inverseProp(iv, distrFrg);
                fg.arcTypeMap[aTarget] = FlowGraph::ArcType::TargetArc;
            }
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
                                                             UnsignedSCellSet& visitedNodes,
                                                             double* distrBkg)
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

            fg.arcWeightMap[a] = infWeigth;
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

    void FlowGraphBuilder::createArcFromLinel(Curve::SCell& linel,
                                              double weight,
                                              FlowGraph::ArcType at)
    {
        ListDigraph::Arc a;
        createArcFromLinel(a,linel,weight,at);
    }

    void FlowGraphBuilder::createArcFromLinel(ListDigraph::Arc& arc,
                                              Curve::SCell& linel,
                                              double weight,
                                              FlowGraph::ArcType at)
    {
        KSpace& KImage = imageFlowData.getKSpace();

        Z2i::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
        Z2i::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

        DGtal::Z2i::Point p = KImage.sCoords( KImage.sDirectIncident(linel,1) );


        createArcFromPixels(arc,
                            directPixel,
                            indirectPixel,
                            at,
                            weight);

        fg.arcSCellMap[arc] = linel;
        fg.scellArcMap[linel] = arc;
    }


    void FlowGraphBuilder::identifySuperposedLinels(Curve::ConstIterator curveBegin,
                                                    Curve::ConstIterator curveEnd,
                                                    std::map<SCell,bool>& superposedLinels)
    {

        for(auto it=curveBegin;it!=curveEnd;++it){
            if(superposedLinels.find(*it)!=superposedLinels.end())
            {
                if(superposedLinels[*it]==false)
                {
                    superposedLinels[*it] = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                superposedLinels[*it] = false;
            }

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


                KSpace::SCell directIncident = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
                KSpace::SCell indirectIncident = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));

                //Source Arc
                ListDigraph::Arc aSource;

                connectNodeToPixel(fg.sourceNode,directIncident,aSource);

                fg.arcWeightMap[aSource] = 1.0/weightMap[linel];
                fg.arcTypeMap[aSource] = FlowGraph::ArcType::SourceArc;

                //Target Arc
                ListDigraph::Arc aTarget;

                connectPixelToNode(indirectIncident,fg.targetNode,aTarget);

                fg.arcWeightMap[aTarget] = 1.0/weightMap[linel];
                fg.arcTypeMap[aTarget] = FlowGraph::ArcType::TargetArc;


                visitedNodes.insert( KImage.sIndirectIncident(linel,KImage.sOrthDir(linel)) );

                if(linkIt==begin.connectorsEnd()) break;
                linkIt++;
            }while(true);

        }

    }

    void FlowGraphBuilder::createEscapeArcs(Curve& fromCurve,
                                            Curve& toCurve,
                                            UnsignedSCellSet& visitedNodes,
                                            std::map<SCell,bool>& superposedLinels)
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

            if(superposedLinels.find(*it)!=superposedLinels.end())
            {
                if(superposedLinels[*it]==true) continue;
            }

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
                    createArcFromPixels(candidate,*n,FlowGraph::ArcType::EscapeArc,infWeigth);
                    createArcFromPixels(*n,candidate,FlowGraph::ArcType::EscapeArc,infWeigth);
                }

                continue;
            }

            explore(*n,borderPoints,visitedNodes);

            createArcFromPixels(candidate,*n,FlowGraph::ArcType::EscapeArc,infWeigth);
            createArcFromPixels(*n,candidate,FlowGraph::ArcType::EscapeArc,infWeigth);
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
                       infWeigth);


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

    double FlowGraphBuilder::dataTerm(int x, bool norm)
    {
        double nf = norm?0.64:1;
        double nt = norm?0.36:0;

        return -alpha*log(nt + nf*x);
    }

    double FlowGraphBuilder::directProp(int iv, double* distr, bool norm)
    {
        double x = 1-distr[iv];
        return dataTerm(x,norm);
    }

    double FlowGraphBuilder::inverseProp(int iv, double* distr, bool norm)
    {
        double x = distr[iv];
        return dataTerm(x,norm);
    }

    double FlowGraphBuilder::similarityTerm(int iv1, int iv2)
    {
        double niv1 = iv1/255.0;
        double niv2 = iv2/255.0;

        return beta*exp( -pow(iv1-iv2,2) );
    }

    void FlowGraphBuilder::similarityBetweenDilated(Curve& curve)
    {
        KSpace& KImage = imageFlowData.getKSpace();
        KSpace::SCell pixelModel = KImage.sCell( KSpace::Point(1,1),KSpace::POS);

        Curve::Point cPixel,nPixel;
        auto it = curve.getInnerPointsRange().begin();
        cPixel = *it;
        do
        {
            ++it;

            if(it==curve.getInnerPointsRange().end())
                it = curve.getInnerPointsRange().begin();

            nPixel = *it;

            int iv1 = imageFlowData.getRefImage().operator()(cPixel);
            int iv2 = imageFlowData.getRefImage().operator()(nPixel);

            double st = similarityTerm(iv1,iv2);


            ListDigraph::Arc arc;
            ListDigraph::Node node1 = fg.coordToNode[KImage.sKCoords( KImage.sCell(cPixel,pixelModel) )];
            ListDigraph::Node node2 = fg.coordToNode[KImage.sKCoords( KImage.sCell(nPixel,pixelModel) )];

            addArc(arc,node1,node2,FlowGraph::ArcType::SimilarityArc,st);
            addArc(arc,node2,node1,FlowGraph::ArcType::SimilarityArc,st);

            cPixel = nPixel;
        }while(it!=curve.getInnerPointsRange().begin());
    }

    void FlowGraphBuilder::similarityBetweenCurves(Curve& curve)
    {
        KSpace& KImage = imageFlowData.getKSpace();

        for(auto it=curve.begin();it!=curve.end();++it)
        {
            Curve::SCell linel = *it;

            Curve::SCell directPixel = KImage.sDirectIncident(linel,KImage.sOrthDir(linel));
            Curve::SCell indirectPixel = KImage.sIndirectIncident(linel,KImage.sOrthDir(linel));


            int iv1 = imageFlowData.getRefImage().operator()( KImage.sCoords(directPixel) );
            int iv2 = imageFlowData.getRefImage().operator()( KImage.sCoords(indirectPixel) );

            double st = similarityTerm(iv1,iv2);

            ListDigraph::Arc arc;
            ListDigraph::Node node1 = fg.coordToNode[ KImage.sKCoords(directPixel) ];
            ListDigraph::Node node2 = fg.coordToNode[ KImage.sKCoords(indirectPixel) ];

            addArc(arc,node1,node2,FlowGraph::ArcType::SimilarityArc,st);
            addArc(arc,node2,node1,FlowGraph::ArcType::SimilarityArc,st);

        }

    }

    void FlowGraphBuilder::createDilatedPixelsTerminalArcs(Curve& innerCurve,
                                                           Curve& outerCurve,
                                                           LinelWeightMap& weightMap,
                                                           double* distrFrg,
                                                           double* distrBkg)
    {
        KSpace& KImage = imageFlowData.getKSpace();

        LinelWeightMap directIncidentPixelWeight,indirectIncidentPixelWeight;
        for(auto it=outerCurve.begin();it!=outerCurve.end();++it)
        {
            Curve::SCell pixel = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
            directIncidentPixelWeight[pixel] = 0.1; //10 MAXIMUM
            indirectIncidentPixelWeight[pixel] = 0.1;
        }

        //Source Arc
        {
            for (auto it = outerCurve.begin(); it != outerCurve.end(); ++it) {
                Curve::SCell pixel = KImage.sDirectIncident(*it, KImage.sOrthDir(*it));
                directIncidentPixelWeight[pixel] += weightMap[*it];
            }

            Curve::SCell pixelModel = KImage.sCell(KSpace::Point(1, 1), KSpace::POS);
            for (auto it = outerCurve.getInnerPointsRange().begin();
                 it != outerCurve.getInnerPointsRange().end();
                 ++it) {
                Curve::SCell pixel = KImage.sCell(*it, pixelModel);

                int iv = imageFlowData.getRefImage().operator()(*it);

                ListDigraph::Arc aSource;
                connectNodeToPixel(fg.sourceNode, pixel, aSource);

                fg.arcWeightMap[aSource] = inverseProp(iv, distrBkg) + 1.0 / directIncidentPixelWeight[pixel];
                fg.arcTypeMap[aSource] = FlowGraph::ArcType::SourceArc;

            }
        }

        //Target Arc
        {
            for (auto it = innerCurve.begin(); it != innerCurve.end(); ++it) {
                Curve::SCell pixel = KImage.sIndirectIncident(*it, KImage.sOrthDir(*it));
                indirectIncidentPixelWeight[pixel] += weightMap[*it];
            }

            Curve::SCell pixelModel = KImage.sCell(KSpace::Point(1, 1), KSpace::POS);
            for (auto it = outerCurve.getInnerPointsRange().begin();
                 it != outerCurve.getInnerPointsRange().end();
                 ++it) {
                Curve::SCell pixel = KImage.sCell(*it, pixelModel);

                int iv = imageFlowData.getRefImage().operator()(*it);

                ListDigraph::Arc aTarget;
                connectPixelToNode(pixel,fg.targetNode, aTarget);

                fg.arcWeightMap[aTarget] = inverseProp(iv, distrFrg) + 1.0 / indirectIncidentPixelWeight[pixel];
                fg.arcTypeMap[aTarget] = FlowGraph::ArcType::SourceArc;
            }
        }
    }

    void FlowGraphBuilder::createOriginalPixelsTerminalArc(Curve& innerCurve,
                                                           LinelWeightMap& weightMap,
                                                           double* distrFrg,
                                                           double* distrBkg)
    {
        KSpace& KImage = imageFlowData.getKSpace();

        LinelWeightMap directIncidentPixelWeight;
        for(auto it=innerCurve.begin();it!=innerCurve.end();++it)
        {
            Curve::SCell pixel = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
            directIncidentPixelWeight[pixel] = 0.1;
        }

        //Source Arc
        {
            for (auto it = innerCurve.begin(); it != innerCurve.end(); ++it) {
                Curve::SCell pixel = KImage.sDirectIncident(*it, KImage.sOrthDir(*it));
                directIncidentPixelWeight[pixel] += weightMap[*it];
            }

            Curve::SCell pixelModel = KImage.sCell(KSpace::Point(1, 1), KSpace::POS);
            for (auto it = innerCurve.getInnerPointsRange().begin();
                 it != innerCurve.getInnerPointsRange().end();
                 ++it) {
                Curve::SCell pixel = KImage.sCell(*it, pixelModel);

                int iv = imageFlowData.getRefImage().operator()(*it);

                ListDigraph::Arc aSource;
                connectNodeToPixel(fg.sourceNode, pixel, aSource);

                fg.arcWeightMap[aSource] = inverseProp(iv, distrBkg) + 1.0 / directIncidentPixelWeight[pixel];
                fg.arcTypeMap[aSource] = FlowGraph::ArcType::SourceArc;

            }
        }

        //Target Arc
        {
            Curve::SCell pixelModel = KImage.sCell(KSpace::Point(1, 1), KSpace::POS);
            for (auto it = innerCurve.getInnerPointsRange().begin();
                 it != innerCurve.getInnerPointsRange().end();
                 ++it) {
                Curve::SCell pixel = KImage.sCell(*it, pixelModel);

                int iv = imageFlowData.getRefImage().operator()(*it);

                ListDigraph::Arc aTarget;
                connectPixelToNode(pixel,fg.targetNode, aTarget);

                fg.arcWeightMap[aTarget] = inverseProp(iv, distrFrg);
                fg.arcTypeMap[aTarget] = FlowGraph::ArcType::TargetArc;

            }
        }

    }

}