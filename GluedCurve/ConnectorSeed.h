#ifndef SEGBYCUT_CONNECTORSEED_H
#define SEGBYCUT_CONNECTORSEED_H

enum ConnectorType{
    internToExtern,
    externToIntern,
    makeConvex
};

template<typename T1,
         typename T2,
        typename T3>
struct ConnectorSeed{
    typedef T1 ConnectorElementType;
    typedef T2 SCellCirculatorType;
    typedef T3 KSpace;


    std::vector<ConnectorElementType> connectors;

    typedef typename SCellCirculatorType::value_type SCell;
    typedef DGtal::Dimension Dimension;


    bool isValid;

    KSpace KImage;


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

    ConnectorSeed(SCellCirculatorType& fit,
                  SCellCirculatorType& sit,
                  int radius): firstCirculator(fit),
                               secondCirculator(sit),
                               cType(makeConvex)
    {
        SCell start = *fit;
        SCell end = *sit;

        SCell commonPointel = KImage.sIndirectIncident(end, *KImage.sDirs(end));

        Dimension d = *KImage.sDirs(start);
        SCell first = start;
        int v = KImage.sSign(first) ? -1 : 1;

        std::vector<SCell> curveSCells;
        int i = 0;
        do {
            curveSCells.push_back(first);
            first = KImage.sGetAdd(first, d, v);

            ++i;
        } while (KImage.sDirectIncident(first, d).preCell().coordinates != commonPointel.preCell().coordinates &&
                 i <= radius);

        if (KImage.sDirectIncident(first, d).preCell().coordinates != commonPointel.preCell().coordinates) {
            isValid = false;
        } else {
            isValid = true;
        }


//        if (isValid) {
//            Curve C;
//            try {
//                C.initFromSCellsRange(curveSCells.begin(), curveSCells.end());
//                return true;
//            } catch (const std::exception &e) {
//                return false;
//            }
//        }
    }
};

template< typename X >
struct ConnectorSeedConcept {
    typedef typename X::ConnectorElementType ConnectorElementType;
    typedef typename X::SCellCirculatorType SCellCirculatorType;
};

#endif //SEGBYCUT_CONNECTORSEED_H
