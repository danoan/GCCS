#include <boost/filesystem/operations.hpp>
#include "DGtal/helpers/StdDefs.h"

#include <ConnectorSeedRange.h>
#include "../FlowGraph/ImageFlowData.h"
#include "../FlowGraph/weightSettings.h"
#include "../ExhaustiveSearch/CombinationsEvaluator.h"
#include "../ExhaustiveSearch/PropertyChecker/MarkedMapChecker/IntervalChecker.h"
#include "../ExhaustiveSearch/SeparateInnerAndOuter.h"
#include "../ExhaustiveSearch/ContainerCombinator.h"
#include "../utils/Artist.h"

typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
typedef DGtal::Z2i::Curve Curve;

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;
};

template<typename CIT>
int sizeCirc(CIT cbegin, CIT cend)
{
    int n=0;
    do
    {
        ++cbegin;
        ++n;
    }while(cbegin!=cend);

    return n;
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


KSpace::SCell scellFromIncidentPoints(const KSpace& KImage,
                                      const KSpace::Point& innerP1, //Not Khalimsky Coordinates
                                      const KSpace::Point& outerP2) //Not Khalimsky Coordinates
{
    KSpace::SCell pixelModelS1 = KImage.sCell(KSpace::Point(1,1),true);
    KSpace::SCell pixelModelS2 = KImage.sCell(KSpace::Point(1,1),false);

    KSpace::SCell s1 = KImage.sCell(innerP1,pixelModelS1);
    KSpace::SCell s2 = KImage.sCell(outerP2,pixelModelS2);

    SCells scellCollS1 = KImage.sLowerIncident(s1);
    SCells scellCollS2 = KImage.sLowerIncident(s2);
    for(auto it1=scellCollS1.begin();it1!=scellCollS1.end();++it1)
    {
        for(auto it2=scellCollS2.begin();it2!=scellCollS2.end();++it2)
        {
            if(it1->preCell().coordinates==it2->preCell().coordinates)
            {
                return *it1;
            }
        }
    }
}

void findBestSCellPair(KSpace::SCell& s1,
                       KSpace::SCell& s2,
                       const KSpace& KImage,
                       const Curve& curve)
{
    typedef typename Curve::IncidentPointsRange::ConstCirculator AdapterCirculator;
    typedef StabbingCircleComputer<AdapterCirculator> SegmentComputer;
    typedef CurvatureFromDCAEstimator<SegmentComputer, false> SCFunctor;

    Curve::IncidentPointsRange intRange = curve.getIncidentPointsRange();

    typedef DGtal::SaturatedSegmentation<SegmentComputer> Segmentation;
    Segmentation seg( intRange.c(), intRange.c(), SegmentComputer() );


    Segmentation::SegmentComputerIterator it = seg.begin();


    int windowSize = 1+3;
    std::list<SegmentComputer> q;

    for(int i=0;i<windowSize;++i,++it)
    {
        q.push_back(*it);
    }
    Segmentation::SegmentComputerIterator itEnd = it;

    int minimumDet=100;
    SegmentComputer chosen[2];
    do
    {
        SegmentComputer scc = q.front();
        q.pop_front();

        int currentSize = sizeCirc(scc.begin(),scc.end());

        for(std::list<SegmentComputer>::const_iterator itn=q.begin();itn!=q.end();++itn)
        {
            int nsize = sizeCirc(itn->begin(),itn->end());
            int diffAbs = std::abs(nsize-currentSize);

            if(std::min(currentSize,nsize)+diffAbs<minimumDet)
            {
                minimumDet =std::min(currentSize,nsize)+diffAbs;
                chosen[0] = scc;
                chosen[1] = *itn;
            }
        }

        q.push_back(*it);
        ++it;
        if(it==seg.end()) it = seg.begin();
    }while(it!=itEnd);

    std::cout << sizeCirc(chosen[0].begin(),chosen[0].end()) << std::endl;
    std::cout << sizeCirc(chosen[0].begin(),chosen[0].end()) << std::endl;

    s1 = scellFromIncidentPoints(KImage,chosen[0].begin()->first,chosen[0].begin()->second);

    AdapterCirculator beforeEndSecond = chosen[1].begin();
    AdapterCirculator temp = beforeEndSecond;
    ++temp;
    while(temp!=chosen[1].end())
    {
        ++beforeEndSecond;
        ++temp;
    }

    s2 = scellFromIncidentPoints(KImage,beforeEndSecond->first,beforeEndSecond->second);
}


typedef Circulator<Curve::ConstIterator> CurveCirculator;
IntervalChecker* selectMinimizationRegion(CurveCirculator& start,
                                          CurveCirculator& end,
                                          const KSpace& KImage,
                                          const Curve& curve)
{
    KSpace::SCell s1,s2;
    findBestSCellPair(s1,s2,KImage,curve);
    CurveCirculator cBegin(curve.begin(),curve.begin(),curve.end());

    auto it = cBegin;
    do
    {
        if(*it==s1)
        {
            start = it;
        }
        if(*it==s2)
        {
            end = it;
            ++end;
        }
        ++it;
    }while(it!=cBegin);


    Board2D board;
    auto itj=start;
    do
    {
        board << *itj;
        ++itj;
    }while(itj!=end);


    return new IntervalChecker(start,end,5);
}

void findOneExpansionMinimumEnergy(Curve& minCurve,
                                   IntervalChecker* intervalChecker,
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
    CE.addChecker( intervalChecker );
    CE(minCurve,filteredPairList,KImage,maxSimultaneousPairs);

}

int main()
{
    std::string filepath = "../images/flow-evolution/single_square.pgm";
    std::string outputFolder = "../output/smallestStabPair/fromCurve/";

    boost::filesystem::create_directories(outputFolder);

    Image2D image = DGtal::GenericReader<Image2D>::import(filepath);

    int i =0;
    while(i<100)
    {
        ImageFlowData imf(image);
        imf.init(ImageFlowData::DilationOnly, 5);


        CurveCirculator start, end;
        IntervalChecker *intervalChecker = selectMinimizationRegion(start,
                                                                    end,
                                                                    imf.getKSpace(),
                                                                    imf.getMostInnerCurve());

        Curve minCurve;
        findOneExpansionMinimumEnergy(minCurve,
                                      intervalChecker,
                                      imf.getKSpace(),
                                      imf.getMostInnerCurve(),
                                      imf.getMostOuterCurve());

        constructImageFromCurve(image,image.domain(),minCurve);
        DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i+1) + ".pgm",image);
        ++i;
    }



    return 0;
}