#include "ConnectorSeedRange.h"

template<typename CellularSpace, typename TIterator>
ConnectorSeedRange<CellularSpace,
                   TIterator>::ConnectorSeedRange(KSpace& KImage,
                                                  const SCellCirculatorType& internalCurveCirculator,
                                                  const SCellCirculatorType& externalCurveCirculator):
        KImage(KImage),
        internalCurveCirculator(internalCurveCirculator),
        externalCurveCirculator(externalCurveCirculator)
{
    setPointelGroup(internalCurveCirculator,POINTEL_GROUP_INTERNAL_CURVE);
    setPointelGroup(externalCurveCirculator,POINTEL_GROUP_EXTERNAL_CURVE);

    SCellCirculatorType ie = externalCurveCirculator;
    SCellCirculatorType ii = internalCurveCirculator;

    //Connection edges creation
    if(DGtal::isNotEmpty(externalCurveCirculator,externalCurveCirculator)){

        do
        {
            SCell externalLinel = *ie;

            SCell sourceExtLinel = KImage.sIndirectIncident(externalLinel,*KImage.sDirs(externalLinel));    //Source
            SCell targetExtLinel = KImage.sDirectIncident(externalLinel,*KImage.sDirs(externalLinel));    //Target

            SCell mainSpel = KImage.sDirectIncident(*ie,KImage.sOrthDir(externalLinel));
            SCells neighborhood = KImage.sNeighborhood(mainSpel);


            for(auto itp=neighborhood.begin();itp!=neighborhood.end();++itp) {
                SCells incidentEdges = KImage.sLowerIncident(*itp);

                for (auto ite = incidentEdges.begin(); ite != incidentEdges.end(); ite++) {
                    SCell potentialConnectionLinel = *ite;

                    if(gridLinels.find(potentialConnectionLinel)!=gridLinels.end()) continue;
                    SCell otherDirection(potentialConnectionLinel);
                    KImage.sSetSign(otherDirection,!KImage.sSign(otherDirection));
                    if(gridLinels.find(otherDirection)!=gridLinels.end()) continue;

                    SCell p1 = KImage.sIndirectIncident(potentialConnectionLinel,*KImage.sDirs(potentialConnectionLinel));    //Source
                    SCell p2 = KImage.sDirectIncident(potentialConnectionLinel,*KImage.sDirs(potentialConnectionLinel));      //Target

                    if ((pointelGroup.find(p1.preCell().coordinates) == pointelGroup.end()) ||
                        (pointelGroup.find(p2.preCell().coordinates) == pointelGroup.end())) {
                        continue;   //There is a pointel which does not belong to any PointelGroup
                    }

                    if ((pointelGroup[p1.preCell().coordinates] == pointelGroup[p2.preCell().coordinates] ) )
                    {
                        continue;   //Same PointelGroup. It is not a connection linel
                    }

                    bool from_intern_to_extern = pointelGroup[p1.preCell().coordinates ]==POINTEL_GROUP_INTERNAL_CURVE;

                    alignIterator(ii,potentialConnectionLinel);

                    if( from_intern_to_extern ){  //Internal Curve is Source

                        if( p2.preCell().coordinates != sourceExtLinel.preCell().coordinates ){
                            continue;   //Conectivity Error
                        }

                        connectorSeedList.push_back( ConnectorSeedType(potentialConnectionLinel,
                                                                       ii,
                                                                       ie,
                                                                       ConnectorType::internToExtern) );
                    }else{

                        if( p1.preCell().coordinates != targetExtLinel.preCell().coordinates ){
                            continue;   //Conectivity Error
                        }

                        connectorSeedList.push_back( ConnectorSeedType(potentialConnectionLinel,
                                                                       ie,
                                                                       ii,
                                                                       ConnectorType::externToIntern)
                        );
                    }

                }
            }


            ++ie;
        }while(ie!=externalCurveCirculator);
    }

    extensionConnectors();
}

template<typename CellularSpace, typename TIterator>
void ConnectorSeedRange<CellularSpace, TIterator>::createExtensionConnectorsSeeds(std::vector<MatchPair>& vMP,
                                                                                  int radius)
{
    typedef typename std::vector<MatchPair>::const_iterator MatchIteratorType;


    for(MatchIteratorType x = vMP.begin();x!=vMP.end();++x){
        SCell start = *x->itb;
        SCell end = *x->ite;

        SCell commonPointel = KImage.sIndirectIncident(end,*KImage.sDirs(end));
        Dimension d = *KImage.sDirs(start);

        SCell first = start;
        int v = KImage.sSign(first)?-1:1;

        std::vector<SCell> curveSCells;
        int i =0;
        do {
            first = KImage.sGetAdd(first,d,v);
            curveSCells.push_back(first);

            ++i;
            std::cout << KImage.sDirectIncident(first,d) << std::endl;

        }while(KImage.sDirectIncident(first,d).preCell().coordinates!=commonPointel.preCell().coordinates
               && i <= radius);

        //I can also check for the connectivity of the curve here
        if(KImage.sDirectIncident(first,d).preCell().coordinates==commonPointel.preCell().coordinates){

            connectorSeedList.push_back( ConnectorSeedType(curveSCells.begin(),
                                                           curveSCells.end(),
                                                           x->itb,
                                                           x->ite)
            );

        }
    }


};

template<typename CellularSpace, typename TIterator>
void ConnectorSeedRange<CellularSpace, TIterator>::extensionConnectors()
{
    int radius = 8;
    std::list<SCellCirculatorType> L;
    std::vector<MatchPair> vMP;

    SCellCirculatorType it = internalCurveCirculator;

    int i =0;
    SCellCirculatorType b=it;
    do{
        L.push_back(b);
        --b;
        ++i;
    }while(b!=internalCurveCirculator && i<radius);

    do{
        validateStack(L,vMP,radius);
        ++it;
        L.push_front(it);
    }while(it!=internalCurveCirculator);

    createExtensionConnectorsSeeds(vMP,radius);
};


template<typename CellularSpace, typename TIterator>
void ConnectorSeedRange<CellularSpace, TIterator>::validateStack(std::list<SCellCirculatorType>& L,
                                                                 std::vector<MatchPair>& vMP,
                                                                 int radius)
{
    SCellCirculatorType t = L.front();

    if(L.size()>radius){
        L.pop_back();
    }

    auto it =L.begin();
    do{
        MatchPair mp;
        if( isItMatch(t,*it,mp,radius) ){
            it=L.erase(it);
            vMP.push_back(mp);
            std::cout << *(vMP[vMP.size()-1].itb) << "::" << *(vMP[vMP.size()-1].ite) << std::endl;
        }
        ++it;
    }while(it!=L.end());

}

template<typename CellularSpace, typename TIterator>
LinelType ConnectorSeedRange<CellularSpace, TIterator>::getLinelType(const SCell& linel)
{
    unsigned int x = KImage.sKCoords(linel)[0];
    unsigned int y = KImage.sKCoords(linel)[1];

    bool xEven = x%2==0;
    bool yEven = y%2==0;
    bool positive = KImage.sSign(linel);

    if( !xEven && yEven){
        if(positive) return LinelType::Left;
        else return LinelType::Right;
    }

    if( xEven && !yEven ){
        if(positive) return LinelType::Down;
        else return LinelType::Up;
    }

};

template<typename CellularSpace, typename TIterator>
bool ConnectorSeedRange<CellularSpace, TIterator>::isItMatch(SCellCirculatorType& closure,
                                                             SCellCirculatorType& el,
                                                             MatchPair& mp,
                                                             int radius)
{
    LinelType ltClosure = getLinelType(*closure);
    LinelType ltEl = getLinelType(*el);


    if(ltClosure==( (ltEl+3)%4 ) )
    {
        Dimension changeDir = *KImage.sDirs(*el);
        Dimension zeroDir = (changeDir-1)%2;

        LinelType ltNext = getLinelType( *(el+1) );
        if( ltNext==ltEl || ltNext==ltClosure ) return false;

        Point v = KImage.sKCoords(*closure) - KImage.sKCoords(*el);
        if( abs(v[zeroDir])!=1) return false;

        mp.itb = el;
        mp.ite = closure;
        mp.ite++;

        return true;

    }

    return false;
}

template<typename CellularSpace, typename TIterator>
void ConnectorSeedRange<CellularSpace,
                        TIterator>::alignIterator(SCellCirculatorType& internalCirculator,
                                                  SCell& connectorLinel)
{
    SCell p1 = KImage.sIndirectIncident(connectorLinel,*KImage.sDirs(connectorLinel));    //Source
    SCell p2 = KImage.sDirectIncident(connectorLinel,*KImage.sDirs(connectorLinel));      //Target

    bool from_intern_to_extern = pointelGroup[p1.preCell().coordinates]==POINTEL_GROUP_INTERNAL_CURVE;

    --internalCirculator;   //To not incur the risk to traverse all the scells.

    SCellCirculatorType endInternalCirculator = internalCirculator;
    if(from_intern_to_extern) {
        do {
            SCell pInt = KImage.sDirectIncident(*internalCirculator, *KImage.sDirs(*internalCirculator));
            if (pInt.preCell().coordinates == p1.preCell().coordinates)//Internal target equals Connector source.
            {
                break;
            }
            ++internalCirculator;
        } while (internalCirculator!=endInternalCirculator);
    }else{
        do {
            SCell pInt = KImage.sIndirectIncident(*internalCirculator, *KImage.sDirs(*internalCirculator));
            if (pInt.preCell().coordinates == p2.preCell().coordinates)//Internal source equals Connector target.
            {
                break;
            }
            ++internalCirculator;
        } while (internalCirculator!=endInternalCirculator);
    }


}


template<typename CellularSpace, typename TIterator>
void ConnectorSeedRange<CellularSpace,
                        TIterator>::setPointelGroup(const SCellCirculatorType& curveCirculator,
                                                    int pointelGroupId)
{
    if( DGtal::isNotEmpty(curveCirculator,curveCirculator) ) {
        SCellCirculatorType ic = curveCirculator;
        do {
            gridLinels.insert(*ic);
            SCells incidentPointels = KImage.sLowerIncident(*ic);

            SCell p = KImage.sDirectIncident( *ic, *KImage.sDirs(*ic) );   //Target
            pointelGroup[p.preCell().coordinates] = pointelGroupId;

            ++ic;
        } while (ic != curveCirculator);
    }

}


template<typename CellularSpace, typename TIterator>
typename ConnectorSeedRange<CellularSpace,TIterator>::ConnectorSeedIteratorType ConnectorSeedRange<CellularSpace,TIterator>::begin() const
{
    return connectorSeedList.begin();
}

template<typename CellularSpace, typename TIterator>
typename ConnectorSeedRange<CellularSpace,TIterator>::ConnectorSeedIteratorType ConnectorSeedRange<CellularSpace,TIterator>::end() const
{
    return connectorSeedList.end();
}

