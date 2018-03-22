
#include <iterator>

#include "DGtal/helpers/StdDefs.h"

#include <ConnectorSeedRange.h>
#include "../FlowGraph/ImageFlowData.h"

#include "../ExhaustiveSearch/SeparateInnerAndOuter.h"
#include "../ExhaustiveSearch/ContainerCombinator.h"
#include "../ExhaustiveSearch/LazyCombinations.h"
#include "../ExhaustiveSearch/PropertyChecker/CheckableSeedPair.h"
#include "../ExhaustiveSearch/PropertyChecker/MarkedMapChecker/GluedIntersectionChecker.h"
#include "../ExhaustiveSearch/PropertyChecker/MarkedMapChecker/MinimumDistanceChecker.h"
#include "../FlowGraph/weightSettings.h"

void addIntervalSCells(std::vector<KSpace::SCell>& vectorOfSCells,
                       SeparateInnerAndOuter::SCellCirculator begin,
                       SeparateInnerAndOuter::SCellCirculator end)
{
    SeparateInnerAndOuter::SCellCirculator it = begin;
    while (it != end)
    {
        vectorOfSCells.push_back(*it);
        ++it;
    }
    vectorOfSCells.push_back(*it);
}

void addSeedPairSCells(std::vector<KSpace::SCell>& vectorOfSCells,
                       CheckableSeedPair& currentPair,
                       CheckableSeedPair& nextPair)
{
    SeparateInnerAndOuter::ConnectorSeed inToExtSeed = currentPair.data().first;
    SeparateInnerAndOuter::ConnectorSeed extToIntSeed = currentPair.data().second;

    if( currentPair.data().first.cType != ConnectorType::internToExtern)
    {
        std::cout << "ERROR" << std::endl;
    }

    if( currentPair.data().second.cType != ConnectorType::externToIntern)
    {
        std::cout << "ERROR" << std::endl;
    }


    vectorOfSCells.push_back(inToExtSeed.connectors[0]);
    addIntervalSCells(vectorOfSCells,inToExtSeed.secondCirculator,extToIntSeed.firstCirculator);


    SeparateInnerAndOuter::ConnectorSeed nextIntToExtSeed = nextPair.data().first;
    vectorOfSCells.push_back(extToIntSeed.connectors[0]);
    addIntervalSCells(vectorOfSCells,extToIntSeed.secondCirculator,nextIntToExtSeed.firstCirculator);
}

void createCurve(SeparateInnerAndOuter::Curve& curve,
                 CheckableSeedPair* seedPairs,
                 int totalPairs)
{
    typedef KSpace::SCell SCell;
    std::vector<SCell> scells;

    CheckableSeedPair currentPair;
    CheckableSeedPair nextPair;
    for(int i=0;i<totalPairs;++i)
    {
        if(i==(totalPairs-1))
        {
            currentPair = seedPairs[totalPairs-1];
            nextPair = seedPairs[0];
        }else
        {
            currentPair = seedPairs[i];
            nextPair = seedPairs[i+1];
        }

        addSeedPairSCells(scells,currentPair,nextPair);
    }

    curve.initFromSCellsVector(scells);
}

double energyValue(SeparateInnerAndOuter::Curve& curve, std::map<KSpace::SCell,double>& weightMap)
{
    auto it = curve.begin();
    double v=0;
    do
    {
        v+=weightMap[*it];
        ++it;
    }while(it!=curve.end());

    return v;
}


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
    CheckableSeedList pairList;

    for(auto it= pairSeedList.begin();it!=pairSeedList.end();++it)
    {
        if(it->first.connectors[0].preCell().coordinates!=it->second.connectors[0].preCell().coordinates)
        {
            pairList.push_back(*it);
        }
    }


    GluedIntersectionChecker intersectionChecker;
    MinimumDistanceChecker minDistChecker;
    for(auto it=pairList.begin();it!=pairList.end();++it)
    {
        intersectionChecker.unmark(*it);
        minDistChecker.unmark(*it);
    }


    int maxSimultaneousPairs = fromInnerSeeds.size()/4;
    double minEnergyValue = 100;
    double currentEnergyValue;
    SeparateInnerAndOuter::Curve minCurve;
    for(int i=1;i<3;++i)
    {
        LazyCombinations<CheckableSeedList,2> myCombinations(pairList);
        myCombinations.addConsistencyChecker(&intersectionChecker);
        myCombinations.addConsistencyChecker(&minDistChecker);


        CheckableSeedPair seedCombination[2];
        int n=0;
        while(myCombinations.next(seedCombination))
        {
            SeparateInnerAndOuter::Curve curve;
            std::map<KSpace::SCell,double> weightMap;

            createCurve(curve,seedCombination,i);
            setGridCurveWeight(curve,KImage,weightMap);

            currentEnergyValue = energyValue(curve,weightMap);
            if(currentEnergyValue<minEnergyValue)
            {
                std::cout << "Updated min energy value: " << minEnergyValue << " -> " << currentEnergyValue << std::endl;
                minEnergyValue = currentEnergyValue;
                minCurve = curve;
            }

            ++n;
        }
        std::cout << i << "-Combination: " << n << std::endl;
    }

    std::cout << "Min Energy Value: " << minEnergyValue << std::endl;
    Board2D board;
    board << minCurve;
    board.saveEPS("minenergycurve.eps");


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