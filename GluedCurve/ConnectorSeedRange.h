#ifndef SEGBYCUT_GLUEDCURVERANGE_H
#define SEGBYCUT_GLUEDCURVERANGE_H

#include "boost/iterator/iterator_facade.hpp"
#include "boost/iterator/iterator_concepts.hpp"

#include "DGtal/base/IteratorCirculatorTraits.h"
#include "DGtal/base/IteratorFunctions.h"

#include "DGtal/topology/CCellularGridSpaceND.h"

#include <map>
#include <list>

#include "ConnectorSeed.h"

namespace Development{
    extern bool makeConvexArcs;
}

enum LinelType{
    Up=0,Left=1,Down=2,Right=3
};

template<typename CellularSpace, typename TIterator>
class ConnectorSeedRange{
    BOOST_STATIC_ASSERT(( boost::is_same<
                             typename DGtal::IteratorCirculatorTraits<TIterator>::Type,
                             DGtal::CirculatorType>::value ));
public:

    const static int POINTEL_GROUP_INTERNAL_CURVE = 0;
    const static int POINTEL_GROUP_EXTERNAL_CURVE = 1;

    typedef CellularSpace KSpace;

    BOOST_CONCEPT_ASSERT(( DGtal::concepts::CCellularGridSpaceND< KSpace > ));

    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::SCells SCells;

    typedef typename DGtal::Dimension Dimension;
    typedef typename KSpace::Point Point;
    typedef DGtal::Z2i::Curve Curve;


    typedef typename KSpace::Point SCellPointelKey;

    typedef SCell ConnectorElementType;
    typedef TIterator SCellCirculatorType;
    typedef ConnectorSeed<ConnectorElementType,SCellCirculatorType,KSpace> ConnectorSeedType;
    typedef typename std::vector< ConnectorSeedType >::const_iterator ConnectorSeedIteratorType;


    BOOST_STATIC_ASSERT(( boost::is_same<
            typename DGtal::IteratorCirculatorTraits<TIterator>::Value,
            SCell>::value ));


    struct MatchPair{
        SCellCirculatorType itb,ite;
    };


    ConnectorSeedRange(KSpace& KImage,
                    const SCellCirculatorType& internCurveCirculator,
                    const SCellCirculatorType& externCurveCirculator);


    inline ConnectorSeedIteratorType begin() const;
    inline ConnectorSeedIteratorType end() const;


private:

    void identifyAlmostEnclosingLinels(std::set<SCell>& enclosingLinels,
                                       SCellCirculatorType externalCurveCirculator);

    bool enclosingPixelStack(std::queue<SCell> pixels,
                             const std::queue<SCell>& linels);


    void alignIterator(SCellCirculatorType& internalCirculator,
                       SCell& connectorLinel);


    void setPointelGroup(const SCellCirculatorType& curveCirculator,
                         int pointelGroupId);

    void extensionConnectors(SCellCirculatorType itBegin,
                             SCellCirculatorType counterClockCirc,
                             bool counterclock,
                             std::vector<MatchPair>& visitedPairs);

    void validateStack(std::list<SCellCirculatorType>& L,
                       std::vector<MatchPair>& vMP,
                       int radius,
                       bool counterclock
    );

    bool isItMatch(SCellCirculatorType& n, SCellCirculatorType& o, MatchPair& mp, int radius, bool counterclock);

    LinelType getLinelType(const SCell& linel);

    void createExtensionConnectorsSeeds(std::vector<MatchPair>& vMP,
                                        int radius);

    void createExtensionConnectorsSeeds(std::vector<MatchPair>& vMP,
                                        SCellCirculatorType& counterClCirc,
                                        int radius);


private:

    KSpace KImage;
    std::vector<ConnectorSeedType> connectorSeedList;

    std::map<SCellPointelKey,int> pointelGroup;

    SCellCirculatorType internalCurveCirculator,
            externalCurveCirculator;


    std::set<SCell> gridLinels;
};

#include "ConnectorSeedRange.cpp"

#endif