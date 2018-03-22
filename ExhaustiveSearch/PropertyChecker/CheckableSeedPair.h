#ifndef SEGBYCUT_CHECKABLESEEDPAIR_H
#define SEGBYCUT_CHECKABLESEEDPAIR_H

#include <DGtal/helpers/StdDefs.h>
#include <ConnectorSeedRange.h>
#include "CheckableElementConcept.h"

class CheckableSeedPair
{
private:
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    typedef DGtal::Z2i::KSpace KSpace;
    typedef KSpace::SCell SCell;

    typedef ConnectorSeedRange<KSpace,SCellCirculator> MyConnectorSeedRange;
    typedef MyConnectorSeedRange::ConnectorSeedType ConnectorSeed;

public:
    typedef CheckableSeedPair Self;

    typedef std::pair<ConnectorSeed,ConnectorSeed> DataType;
    typedef std::pair<ConnectorSeed,ConnectorSeed> CheckedType;
    typedef SCell MarkedType;

    typedef
    class UnsignedSCellComparison{
    public:
        bool operator()(const MarkedType& s1, const MarkedType& s2) const{
            return s1.preCell().coordinates == s2.preCell().coordinates;
        }
    } ComparisonClass;

    BOOST_CONCEPT_ASSERT( (CheckableElementConcept<Self>) );

public:
    CheckableSeedPair(){}
    CheckableSeedPair(CheckedType data):_data(data){}

    const CheckedType& data() const  { return _data; }

private:
    CheckedType _data;
};

#endif //SEGBYCUT_CHECKABLESEEDPAIR_H
