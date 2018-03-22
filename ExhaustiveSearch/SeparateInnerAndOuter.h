#ifndef SEGBYCUT_SEPARATEINNERANDOUTER_H
#define SEGBYCUT_SEPARATEINNERANDOUTER_H

#include <DGtal/helpers/StdDefs.h>


#include <ConnectorSeedRange.h>

class SeparateInnerAndOuter
{
public:
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    typedef DGtal::Z2i::Curve Curve;
    typedef DGtal::Z2i::KSpace KSpace;

    typedef ConnectorSeedRange<KSpace,SCellCirculator> MyConnectorSeedRange;
    typedef MyConnectorSeedRange::ConnectorSeedType ConnectorSeed;
    typedef MyConnectorSeedRange::ConnectorSeedIteratorType ConnectorSeedIterator;

    typedef std::vector<ConnectorSeed> SeedVector;

public:
    SeparateInnerAndOuter(KSpace& KImage,
                          Curve& innerCurve,
                          Curve& outerCurve):KImage(KImage),
                                             innerCurve(innerCurve),
                                             outerCurve(outerCurve){}

    void operator()(SeedVector& fromInnerSeeds,
                    SeedVector& fromOuterSeeds)
    {
        SCellCirculator innerCirc(innerCurve.begin(),
                             innerCurve.begin(),
                             innerCurve.end());

        SCellCirculator outerCirc(outerCurve.begin(),
                             outerCurve.begin(),
                             outerCurve.end());

        MyConnectorSeedRange csr(KImage,innerCirc,outerCirc);

        innerAndOuterList(fromInnerSeeds,fromOuterSeeds,csr.begin(),csr.end());
    }

private:
    void innerAndOuterList(SeedVector& fromInnerSCells,
                           SeedVector& fromOuterSCells,
                           ConnectorSeedIterator begin,
                           ConnectorSeedIterator end)
    {
        for(ConnectorSeedIterator it = begin;it!=end;++it)
        {
            switch( it->cType )
            {
                case ConnectorType::internToExtern:
                    fromInnerSCells.push_back(*it);
                    break;
                case ConnectorType::externToIntern:
                    fromOuterSCells.push_back(*it);
                    break;
            }
        }
    }

public:
    KSpace& KImage;
    DGtal::Z2i::Curve& innerCurve;
    DGtal::Z2i::Curve& outerCurve;
};

#endif //SEGBYCUT_SEPARATEINNERANDOUTER_H
