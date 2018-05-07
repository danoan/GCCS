#ifndef SEGBYCUT_MINIMUMDISTANCECHECKER_H
#define SEGBYCUT_MINIMUMDISTANCECHECKER_H

#include "../CheckableSeedPair.h"
#include "MarkedMapCheckerInterface.h"

#include "DGtal/helpers/StdDefs.h"

#ifndef CHECKABLESEEDPAIR_HASH
#define CHECKABLESEEDPAIR_HASH
namespace std
{
    template<>
    struct hash<CheckableSeedPair::MarkedType>
    {
        std::size_t operator()(const CheckableSeedPair::MarkedType& k) const
        {
            return ( ( hash<int>()( k.preCell().coordinates[0] )
                       ^( hash<int>()( k.preCell().coordinates[1] ) << 1 ) >> 1 )
                     ^( hash<bool>()( k.preCell().positive) )
            );
        }
    };
}
#endif

class MinimumDistanceChecker: public MarkedMapCheckerInterface<CheckableSeedPair>
{
public:
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    typedef DGtal::Z2i::KSpace KSpace;

private:
    bool validConnector(const KSpace::SCell& c1) const
    {
        KSpace::SCell c1s = KImage.sIndirectIncident(c1,*KImage.sDirs(c1));
        KSpace::SCell c1t = KImage.sDirectIncident(c1,*KImage.sDirs(c1));

        if( this->_markMap.at( c1s ) ) return false;
        if( this->_markMap.at( c1t ) ) return false;

        return true;
    }

    void markConnector(const KSpace::SCell& c1,bool mark)
    {
        KSpace::SCell c1s = KImage.sIndirectIncident(c1,*KImage.sDirs(c1));
        KSpace::SCell c1t = KImage.sDirectIncident(c1,*KImage.sDirs(c1));

        this->_markMap[ c1s ] = mark;
        this->_markMap[ c1t ] = mark;
    }


public:
    MinimumDistanceChecker(const KSpace& KImage):KImage(KImage){}
    bool operator()(const CheckableSeedPair& sp) const
    {
        if( !validConnector( sp.data().first.connectors.at(0) ) ) return false;
        if( !validConnector( sp.data().second.connectors.at(0) ) ) return false;

        return true;
    }

    void mark(const CheckableSeedPair& sp)
    {
        markConnector( sp.data().first.connectors.at(0),true);
        markConnector( sp.data().second.connectors.at(0),true);
    }

    void unmark(const CheckableSeedPair& sp)
    {
        markConnector( sp.data().first.connectors.at(0),false);
        markConnector( sp.data().second.connectors.at(0),false);
    }

private:
    std::unordered_map<
            CheckableSeedPair::MarkedType,
            bool,
            std::hash<CheckableSeedPair::MarkedType>,
            CheckableSeedPair::ComparisonClass>
            _markMap;

    const KSpace& KImage;

};

#endif //SEGBYCUT_MINIMUMDISTANCECHECKER_H
