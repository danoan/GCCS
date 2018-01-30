#ifndef SEGBYCUT_CONNECTORSEED_H
#define SEGBYCUT_CONNECTORSEED_H

enum ConnectorType{
    internToExtern,
    externToIntern
};

template<typename T1,
         typename T2>
struct ConnectorSeed{
    typedef T1 ConnectorElementType;
    typedef T2 SCellCirculatorType;

    std::vector<ConnectorElementType> connectors;
    ConnectorType cType;

    typedef typename std::vector<ConnectorElementType>::const_iterator SCellIteratorType;


    SCellCirculatorType firstCirculator;
    SCellCirculatorType secondCirculator;

    ConnectorSeed(ConnectorElementType cet,
                  SCellCirculatorType& fit,
                  SCellCirculatorType& sit,
                  ConnectorType ct): firstCirculator(fit),
                                     secondCirculator(sit),
                                     cType(ct)
    {
        connectors.push_back(cet);
    };


    ConnectorSeed(SCellIteratorType& cItB,
                  SCellIteratorType& cItE,
                  SCellCirculatorType& fit,
                  SCellCirculatorType& sit,
                  ConnectorType ct): firstCirculator(fit),
                                     secondCirculator(sit),
                                     cType(ct)
    {
        connectors.insert(cItB,cItE);
    };
};

template< typename X >
struct ConnectorSeedConcept {
    typedef typename X::ConnectorElementType ConnectorElementType;
    typedef typename X::SCellCirculatorType SCellCirculatorType;
};

#endif //SEGBYCUT_CONNECTORSEED_H
