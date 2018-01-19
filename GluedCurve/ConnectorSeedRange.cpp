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

