#ifndef SEGBYCUT_MINIMUMDISTANCECHECKER_H
#define SEGBYCUT_MINIMUMDISTANCECHECKER_H

#include "../CheckableSeedPair.h"
#include "MarkedMapCheckerInterface.h"

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

public:
    bool operator()(const CheckableSeedPair& sp) const
    {
        if( this->_markMap.at( sp.data().first.connectors.at(0) ) ) return false;
        if( this->_markMap.at( sp.data().second.connectors.at(0) ) ) return false;

        return true;
    }

    void mark(const CheckableSeedPair& sp)
    {
        this->_markMap[ sp.data().first.connectors.at(0) ] = true;
        this->_markMap[ sp.data().second.connectors.at(0) ] = true;
    }

    void unmark(const CheckableSeedPair& sp)
    {
        this->_markMap[ sp.data().first.connectors.at(0) ] = false;
        this->_markMap[ sp.data().second.connectors.at(0) ] = false;
    }

private:
    std::unordered_map<
            CheckableSeedPair::MarkedType,
            bool,
            std::hash<CheckableSeedPair::MarkedType>,
            CheckableSeedPair::ComparisonClass>
            _markMap;

};

#endif //SEGBYCUT_MINIMUMDISTANCECHECKER_H
