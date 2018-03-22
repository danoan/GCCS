#ifndef SEGBYCUT_GLUEDINTERSECTIONCHECKER_H
#define SEGBYCUT_GLUEDINTERSECTIONCHECKER_H

#include "MarkedMapCheckerInterface.h"
#include "../CheckableSeedPair.h"

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
            );
        }
    };
}
#endif


class GluedIntersectionChecker: public MarkedMapCheckerInterface<CheckableSeedPair>
{
public:
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

public:
    bool operator()(const CheckableSeedPair& sp) const
    {
        SCellCirculator extCircBegin = sp.data().first.secondCirculator; //First -> intToExt
        SCellCirculator extCircEnd = sp.data().second.firstCirculator; //Second -> extToInt

        SCellCirculator it = extCircBegin;
        do
        {
            if( this->_markMap.at(*it) ) return false;
            ++it;
        }while(it!=extCircEnd);

        if( this->_markMap.at(*it) ) return false;
        else return true;
    }

    void mark(const CheckableSeedPair& sp)
    {
        SCellCirculator extCircBegin = sp.data().first.secondCirculator; //First -> intToExt
        SCellCirculator extCircEnd = sp.data().second.firstCirculator; //Second -> extToInt

        SCellCirculator it = extCircBegin;
        do
        {
            this->_markMap[*it] = true;
            ++it;
        }while(it!=extCircEnd);

        this->_markMap[*it] = true;
    }

    void unmark(const CheckableSeedPair& sp)
    {
        SCellCirculator extCircBegin = sp.data().first.secondCirculator; //First -> intToExt
        SCellCirculator extCircEnd = sp.data().second.firstCirculator; //Second -> extToInt

        SCellCirculator it = extCircBegin;
        do
        {
            this->_markMap[*it] = false;
            ++it;
        }while(it!=extCircEnd);

        this->_markMap[*it] = false;
    }

private:
    std::unordered_map<
            CheckableSeedPair::MarkedType,
            bool,
            std::hash<CheckableSeedPair::MarkedType>,
            CheckableSeedPair::ComparisonClass>
            _markMap;

};

#endif //SEGBYCUT_GLUEDINTERSECTIONCHECKER_H
