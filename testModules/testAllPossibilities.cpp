
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
                           Curve& outerCurve, Curve& minCurve,int itnum)
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
        if( connectorsDistance < 4 || connectorsDistance > 10) continue;

        filteredPairList.push_back(*it);

    }

    int maxSimultaneousPairs = fromInnerSeeds.size()/4;
    CombinationsEvaluator<1> CE;
    CE.addChecker( new GluedIntersectionChecker() );
    CE.addChecker( new MinimumDistanceChecker(KImage) );
    //CE(filteredPairList,KImage,maxSimultaneousPairs,"../output/combinations/L" + std::to_string(itnum));

    double currLength = computeLength(innerCurve,KImage);
    CheckableSeedPair seedCombination[10];
    CE(minCurve,
       seedCombination,
       filteredPairList,
       KImage,
       maxSimultaneousPairs,
       "../output/combinations/LFLIM" + std::to_string(itnum),
       currLength);
}

typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
void constructImageFromCurve(Image2D& image,
                             const Domain& domain,
                             const Curve& curve) {
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Curve Curve;

    Curve::InnerPointsRange innerPixelsRange = curve.getInnerPointsRange();
    Curve::OuterPointsRange outerPixelsRange = curve.getOuterPointsRange();

    std::stack<Point> toExplorePoints;
    Point filter[4] = {Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1)};
    DigitalSet dsOut(domain);
    DigitalSet dsForbidden(domain);
    {
        Curve::InnerPointsRange::ConstIterator it = innerPixelsRange.begin();

        do {
            dsForbidden.insert(*it);
            dsOut.insert(*it);
            ++it;
        } while (it != innerPixelsRange.end());
    }

    {
        Curve::OuterPointsRange::ConstIterator it = outerPixelsRange.begin();

        do {
            dsForbidden.insert(*it);
            ++it;
        } while (it != outerPixelsRange.end());
    }


    for (auto it = dsOut.begin(); it != dsOut.end(); ++it)
    {
        toExplorePoints.push(*it);
        for(int i=0;i<4;++i) toExplorePoints.push(*it + filter[i]);
    }

    while(!toExplorePoints.empty())
    {
        Point p = toExplorePoints.top(); toExplorePoints.pop();

        if(p[0]<domain.lowerBound()[0] || p[1]<domain.lowerBound()[1]) continue;
        if(p[0]>domain.upperBound()[0] || p[1]>domain.upperBound()[1]) continue;
        if(dsOut(p)) continue;
        if(dsForbidden(p)) continue;

        dsOut.insert(p);
        for(int i=0;i<4;++i) toExplorePoints.push(p+filter[i]);

    }

    for(auto it=image.domain().begin();it!=image.domain().end();++it)
    {
        image.setValue(*it,0);
    }

    for(auto it=dsOut.begin();it!=dsOut.end();++it)
    {
        image.setValue(*it,255);
    }
}

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};

double computeLength(ImageFlowData& imf)
{

    return 0;
}


int main()
{
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;

    std::string filepath = "../images/flow-evolution/single_square.pgm";
    Image2D image = DGtal::GenericReader<Image2D>::import(filepath);


    double l1,l2;
    for(int i=0;i<30;++i)
    {
        Curve minCurve;

        ImageFlowData imf(image,image);
        imf.init(ImageFlowData::DilationOnly,5);

        checkAllPossibilities(imf.getKSpace(),
                              imf.getMostInnerCurve(),
                              imf.getMostOuterCurve(),minCurve,i );

        constructImageFromCurve(image,image.domain(),minCurve);
    }



    return 0;
}