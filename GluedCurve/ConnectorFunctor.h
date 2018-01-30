#ifndef CONNECTOR_FUNCTION
#define CONNECTOR_FUNCTION

#include "GluedCurveIterator.h"
#include "ConnectorSeed.h"

template<typename MyCirculator>
MyCirculator walkCirculator(const MyCirculator c, int w){
    MyCirculator cc = c;
    while(w!=0){
        w>0?cc++:cc--;
        w>0?w--:w++;
    }

    return cc;
}


template<typename TPConnectorSeedType>
class ConnectorSeedToGluedCurveRange
{
public:
    BOOST_CONCEPT_ASSERT( (ConnectorSeedConcept<TPConnectorSeedType>) );


    typedef TPConnectorSeedType ConnectorSeedType;
    typedef typename ConnectorSeedType::SCellCirculatorType SCellCirculatorType;
    typedef typename ConnectorSeedType::SCellIteratorType SCellIteratorType;

    typedef GluedCurveIterator<SCellCirculatorType,SCellIteratorType > GCIterator;

    typedef ConnectorSeedType Input;
    typedef std::pair<GCIterator,GCIterator> Output;

public:
    ConnectorSeedToGluedCurveRange(unsigned int distance):distance(distance){};

    Output operator()(const ConnectorSeedType& it) const;

private:
    unsigned int distance;
};

template<typename ConnectorSeedIteratorType>
typename ConnectorSeedToGluedCurveRange<ConnectorSeedIteratorType>::Output
ConnectorSeedToGluedCurveRange<ConnectorSeedIteratorType>::operator() (const Input& cc) const
{
    SCellCirculatorType it1b,it1e, it2b,it2e;


    it1b = walkCirculator(cc.firstCirculator, -distance+1);
    it1e = cc.firstCirculator;
    it2b = cc.secondCirculator;
    it2e = walkCirculator(it2b, distance);



    return std::pair<GCIterator,
                     GCIterator>(
            GCIterator(cc.connectors.begin(), --cc.connectors.end(), it1b, it1e, it2b, it2e,cc.cType),
            GCIterator(cc.connectors.begin(), --cc.connectors.end(), it1b, it1e, it2b, it2e,cc.cType,true)
    );
}

#endif