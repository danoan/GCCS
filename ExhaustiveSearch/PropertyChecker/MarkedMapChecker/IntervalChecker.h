#ifndef SEGBYCUT_INTERVALCHECKER_H
#define SEGBYCUT_INTERVALCHECKER_H

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


class IntervalChecker: public MarkedMapCheckerInterface<CheckableSeedPair>
{
public:
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    typedef DGtal::Z2i::Curve::SCell SCell;

public:
    IntervalChecker(SCellCirculator itb, SCellCirculator ite, int gluedCurveLength)
    {
        this->gluedCurveLength = gluedCurveLength;
        auto it = itb;
        do
        {
            _allowedSCells.insert(*it);
            ++it;
        }while(it!=ite);
        _allowedSCells.insert(*ite);
    }
    
    bool operator()(const CheckableSeedPair& sp) const
    {
        SCellCirculator extCircBegin = sp.data().first.firstCirculator; //First -> intToExt
        SCellCirculator extCircEnd = sp.data().second.secondCirculator; //Second -> extToInt
        
        for(int i=0;i<gluedCurveLength;++i) ++extCircBegin;

        if( _allowedSCells.find(*extCircBegin) == _allowedSCells.end() ||
            _allowedSCells.find(*extCircEnd) == _allowedSCells.end() )
        {
            return false;
        }
        else
        {
            auto it = extCircBegin;
            do
            {
                if(_allowedSCells.find(*it)==_allowedSCells.end()) return false;
                ++it;
            }while(it!=extCircEnd);
        }

        return true;
    }

    void mark(const CheckableSeedPair& sp)
    {

    }

    void unmark(const CheckableSeedPair& sp)
    {

    }

private:
    std::unordered_map<
            CheckableSeedPair::MarkedType,
            bool,
            std::hash<CheckableSeedPair::MarkedType>,
            CheckableSeedPair::ComparisonClass>
            _markMap;
    
    std::set<SCell> _allowedSCells;
    int gluedCurveLength;
};


#endif //SEGBYCUT_INTERVALCHECKER_H
