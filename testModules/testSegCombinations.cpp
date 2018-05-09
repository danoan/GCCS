
#include <iterator>

#include "DGtal/helpers/StdDefs.h"

#include <ConnectorSeedRange.h>
#include "../FlowGraph/ImageFlowData.h"
#include "../FlowGraph/weightSettings.h"

#include "../ExhaustiveSearch/SeparateInnerAndOuter.h"
#include "../ExhaustiveSearch/ContainerCombinator.h"
#include "../ExhaustiveSearch/PropertyChecker/CheckableSeedPair.h"
#include "../ExhaustiveSearch/CombinationsEvaluator.h"
#include "../utils/Artist.h"

typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;

void findOneExpansionMinimumEnergy(Curve& minCurve,
                                   KSpace& KImage,
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

        //Internal and External Connector must cover between 2 and at most 10 external linels
        if( connectorsDistance >= 20 || connectorsDistance <= 1) continue;

        filteredPairList.push_back(*it);

    }

    int maxSimultaneousPairs = fromInnerSeeds.size()/4;
    CombinationsEvaluator<1> CE;
    CE.addChecker( new GluedIntersectionChecker() );
    CE.addChecker( new MinimumDistanceChecker(KImage) );
    CE(minCurve,filteredPairList,KImage,maxSimultaneousPairs);

}

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

void computeOneExpansions()
{
    std::string filepath = "../images/flow-evolution/single_square.pgm";
    std::string outputFolder = "../output/segComb/fromCurve/";

    boost::filesystem::create_directories(outputFolder);

    Image2D image = DGtal::GenericReader<Image2D>::import(filepath);

    Curve minCurve;

    int i=0;
    while(i<100)
    {
        ImageFlowData imf(image);
        imf.init(ImageFlowData::DilationOnly, 5);


        findOneExpansionMinimumEnergy(minCurve,
                                      imf.getKSpace(),
                                      imf.getMostInnerCurve(),
                                      imf.getMostOuterCurve());

        constructImageFromCurve(image,image.domain(),minCurve);
        DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i+1) + ".pgm",image);
        ++i;
    }
}

void computeStabbingCircles()
{
    std::string inputFolder = "../output/segComb/fromCurve/";
    std::string outputFolder = "../output/segComb/fromCurve/";

    boost::filesystem::path p = inputFolder;
    boost::filesystem::directory_iterator it(p);

    Board2D board;
    while(it!=boost::filesystem::directory_iterator())
    {
        if(strcmp( it->path().extension().c_str(),".eps" )!=0) {

            Image2D image = DGtal::GenericReader<Image2D>::import(it->path().c_str());
            ImageFlowData imf(image);
            imf.init(ImageFlowData::FlowMode::DilationOnly, 5);


            Artist EA(imf.getKSpace(), board);
            EA.drawMaximalStabbingCircles(imf.getMostInnerCurve());
            std::string currentFilepath = outputFolder + it->path().stem().c_str() + "-stabbing-circles.eps";
            EA.board.saveEPS(currentFilepath.c_str());
        }
        ++it;
    }
}


int main()
{
    computeOneExpansions();
//    computeStabbingCircles();

    return 0;
}