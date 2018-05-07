
#include <iterator>

#include "DGtal/helpers/StdDefs.h"

#include <ConnectorSeedRange.h>
#include "../FlowGraph/ImageFlowData.h"

#include "../ExhaustiveSearch/SeparateInnerAndOuter.h"
#include "../ExhaustiveSearch/ContainerCombinator.h"
#include "../ExhaustiveSearch/PropertyChecker/CheckableSeedPair.h"
#include "../FlowGraph/weightSettings.h"
#include "../ExhaustiveSearch/CombinationsEvaluator.h"


void checkAllPossibilities(KSpace& KImage,
                           Curve& innerCurve,
                           Curve& outerCurve)
{

    typedef SeparateInnerAndOuter::ConnectorSeed ConnectorSeed;
    typedef SeparateInnerAndOuter::ConnectorSeedIterator ConnectorSeedIterator;

    SeparateInnerAndOuter::SeedVector fromInnerSeeds;
    SeparateInnerAndOuter::SeedVector fromOuterSeeds;

    {
        SeparateInnerAndOuter _(KImage,innerCurve,outerCurve);
        _(fromInnerSeeds,fromOuterSeeds);
    }


    typedef std::vector< CheckableSeedPair::DataType > DataType;

    DataType pairSeedList;
    ContainerCombinator<ConnectorSeedIterator,ConnectorSeedIterator> combinator;

    combinator.operator()(fromInnerSeeds.begin(),
                          fromInnerSeeds.end(),
                          fromOuterSeeds.begin(),
                          fromOuterSeeds.end(),
                          std::back_inserter(pairSeedList));


    typedef std::vector< CheckableSeedPair > CheckableSeedList;
    CheckableSeedList filteredPairList;

    for(auto it= pairSeedList.begin();it!=pairSeedList.end();++it)
    {
        DGtal::PointVector<2,int> coordInt = it->first.connectors[0].preCell().coordinates;
        DGtal::PointVector<2,int> coordExt = it->second.connectors[0].preCell().coordinates;
        DGtal::PointVector<2,int> diff = coordExt - coordInt;
        int connectorsDistance = (abs(diff[0])+abs(diff[1]) )/2;

        //Internal Connector must be different from External Connector
        if(coordInt==coordExt) continue;

        //Internal and External Connector must cover between 10 and at most 20 external linels
        if( connectorsDistance >= 20 || connectorsDistance <= 6) continue;

        filteredPairList.push_back(*it);

    }

    int maxSimultaneousPairs = fromInnerSeeds.size()/4;
    CombinationsEvaluator<1> CE;
    CE.addChecker( new GluedIntersectionChecker() );
    CE.addChecker( new MinimumDistanceChecker(KImage) );
    CE(filteredPairList,KImage,maxSimultaneousPairs,"../output/combinations/L");

}

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};


int main()
{
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;

    std::string filepath = "../images/flow-evolution/single_square.pgm";
    Image2D image = DGtal::GenericReader<Image2D>::import(filepath);

    ImageFlowData imf(image);
    imf.init(ImageFlowData::DilationOnly,5);

    checkAllPossibilities(imf.getKSpace(),
                          imf.getMostInnerCurve(),
                          imf.getMostOuterCurve() );


    return 0;
}