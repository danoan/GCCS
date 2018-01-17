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

    ConnectorElementType connector;
    ConnectorType cType;


    SCellCirculatorType firstCirculator;
    SCellCirculatorType secondCirculator;

    ConnectorSeed(ConnectorElementType& cet,
                  SCellCirculatorType& fit,
                  SCellCirculatorType& sit,
                  ConnectorType ct): connector(cet),
                                     firstCirculator(fit),
                                     secondCirculator(sit),
                                     cType(ct)
    {};
};

template< typename X >
struct ConnectorSeedConcept {
    typedef typename X::ConnectorElementType ConnectorElementType;
    typedef typename X::SCellCirculatorType SCellCirculatorType;
};

#endif //SEGBYCUT_CONNECTORSEED_H
