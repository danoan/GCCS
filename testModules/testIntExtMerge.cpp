#include <boost/filesystem/operations.hpp>
#include "DGtal/helpers/StdDefs.h"

#include <ConnectorSeedRange.h>
#include "../FlowGraph/ImageFlowData.h"
#include "../FlowGraph/weightSettings.h"
#include "../ExhaustiveSearch/CombinationsEvaluator.h"
#include "../ExhaustiveSearch/PropertyChecker/MarkedMapChecker/IntervalChecker.h"
#include "../ExhaustiveSearch/SeparateInnerAndOuter.h"
#include "../ExhaustiveSearch/ContainerCombinator.h"
#include "../Artist.h"


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


struct RelevantCurve
{
    Curve curve;
    std::set<Curve::SCell> not_allowed;
    std::vector<CheckableSeedPair> selectedPairs;
};

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


typedef StabbingCircleComputer<AdapterCirculator> SegmentComputer;
bool is_not_allowed(const KSpace& KImage, const SegmentComputer& scc1, const SegmentComputer& scc2, const RelevantCurve& rv)
{
    KSpace::SCell s1 = scellFromIncidentPoints(KImage,scc1.begin()->first,scc1.begin()->second);

    AdapterCirculator beforeEndSecond = scc2.begin();
    AdapterCirculator temp = beforeEndSecond;
    ++temp;
    while (temp != scc2.end()) {
        ++beforeEndSecond;
        ++temp;
    }

    KSpace::SCell s2 = scellFromIncidentPoints(KImage, beforeEndSecond->first, beforeEndSecond->second);

    CurveCirculator cBegin(rv.curve.begin(),rv.curve.begin(),rv.curve.end());
    CurveCirculator start,end;
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

    auto itk = start;
    ++itk;
    --end;
    do
    {
        if( rv.not_allowed.find( *itk )!=rv.not_allowed.end() ) return true;
        ++itk;
    }while(itk!=end);


    return false;
}

class GoodSCellPair
{
public:
    typedef typename Curve::IncidentPointsRange::ConstCirculator AdapterCirculator;
    typedef StabbingCircleComputer<AdapterCirculator> SegmentComputer;
    typedef CurvatureFromDCAEstimator<SegmentComputer, false> SCFunctor;
    typedef DGtal::SaturatedSegmentation<SegmentComputer> Segmentation;

    struct SolutionCandidate
    {
        SolutionCandidate(){}
        SolutionCandidate(SegmentComputer scc1, SegmentComputer scc2, int value):scc1(scc1),
                                                                                 scc2(scc2),
                                                                                 value(value){}
        SegmentComputer scc1;
        SegmentComputer scc2;

        int value;
    };

    struct SolutionCompare
    {
        bool operator()(const SolutionCandidate& s1, const SolutionCandidate& s2)
        {
            return s1.value < s2.value;
        }
    };

public:
    GoodSCellPair(const KSpace& KImage, RelevantCurve& rv):KImage(KImage),
                                                     curve(rv.curve)
    {
        Curve::IncidentPointsRange intRange = curve.getIncidentPointsRange();

        Segmentation seg( intRange.c(), intRange.c(), SegmentComputer() );
        Segmentation::SegmentComputerIterator it = seg.begin();

        maxIterations = curve.size();
        windowSize = 1+2;
        for(auto it=seg.begin();it!=seg.end();++it)
        {
            q.push_back(*it);
        }

        int i=0;
        for(auto it=seg.begin();i<windowSize-1;++i)
        {
            q.push_back(*it);
        }

        findPotentialSolutions();
        solutionCandidates.sort(SolutionCompare());
    }

    bool next(SolutionCandidate& sc)
    {
        if(solutionCandidates.size()==0) return false;
        sc = solutionCandidates.front();
        solutionCandidates.pop_front();

        return true;
    }

private:
    void findPotentialSolutions()
    {
        int currentIteration=0;
        int minimumDet=100;
        do {
            SegmentComputer scc = q.front();
            q.pop_front();

            int currentSize = sizeCirc(scc.begin(), scc.end());

            for (std::list<SegmentComputer>::const_iterator itn = q.begin(); itn != q.end(); ++itn) {
                int nsize = sizeCirc(itn->begin(), itn->end());
                int diffAbs = std::abs(nsize - currentSize);

                solutionCandidates.push_back(SolutionCandidate(scc,
                                                               *itn,
                                                               std::min(currentSize, nsize) + diffAbs) );

            }

            ++currentIteration;
        }while(currentIteration < maxIterations);
    }

private:
    const KSpace& KImage;
    Curve& curve;

    int windowSize;
    int maxIterations;
    std::list<SegmentComputer> q;
    std::list<SolutionCandidate> solutionCandidates;
};


class PrepareSolutionCandidate
{
public:
    typedef GoodSCellPair::SolutionCandidate SolutionCandidate;

    typedef typename Curve::IncidentPointsRange::ConstCirculator AdapterCirculator;
    typedef Circulator<Curve::ConstIterator> CurveCirculator;
    typedef KSpace::SCell SCell;

public:
    PrepareSolutionCandidate(const KSpace& KImage,
                             const Curve& curve,
                             int gluedCurveLength,
                             SolutionCandidate& sc)
    {
        AdapterCirculator incidentPointsSeg1Begin = sc.scc1.begin();
        AdapterCirculator incidentPointsSeg1End = sc.scc1.end();

        AdapterCirculator incidentPointsSeg2Begin = sc.scc2.begin();
        AdapterCirculator incidentPointsSeg2End = sc.scc2.begin();

        std::cout << sizeCirc(incidentPointsSeg1Begin, incidentPointsSeg1End) << std::endl;
        std::cout << sizeCirc(incidentPointsSeg2Begin, incidentPointsSeg2End) << std::endl;

        SCell s1 = scellFromIncidentPoints(KImage,
                                           incidentPointsSeg1Begin->first,
                                           incidentPointsSeg1Begin->second);

        AdapterCirculator beforeEndSecond = incidentPointsSeg2Begin;
        AdapterCirculator temp = beforeEndSecond;
        ++temp;
        while (temp != incidentPointsSeg2End) {
            ++beforeEndSecond;
            ++temp;
        }

        SCell s2 = scellFromIncidentPoints(KImage, beforeEndSecond->first, beforeEndSecond->second);
        createIntervalChecker(curve,s1,s2);
    }

private:
    void createIntervalChecker(const Curve& curve, const SCell& s1, const SCell& s2)
    {
        CurveCirculator cBegin(curve.begin(),curve.begin(),curve.end());
        CurveCirculator start;
        CurveCirculator end;

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

        intervalChecker = new IntervalChecker(start,end,5);
    }

public:
    IntervalChecker* intervalChecker;

};


IntervalChecker* selectMinimizationRegion(CurveCirculator& start,
                                          CurveCirculator& end,
                                          bool& successFlag,
                                          const KSpace& KImage,
                                          RelevantCurve& rv)
{
    GoodSCellPair goodPair(KImage,rv);
    GoodSCellPair::SolutionCandidate sc;
    while( goodPair.next(sc) )
    {
        PrepareSolutionCandidate psc(KImage,rv.curve,5,sc);
    }
}

void findOneExpansionMinimumEnergy(Curve& minCurve,
                                   IntervalChecker* intervalChecker,
                                   RelevantCurve& rv,
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
        if( connectorsDistance >= 20 || connectorsDistance <= 3) continue;

        filteredPairList.push_back(*it);

    }

    int maxSimultaneousPairs = fromInnerSeeds.size()/4;
    CombinationsEvaluator<1> CE;
    CE.addChecker( new GluedIntersectionChecker() );
    CE.addChecker( new MinimumDistanceChecker(KImage) );
    CE.addChecker( intervalChecker );

    CheckableSeedPair seedCombination[10];
    CE(minCurve,seedCombination,filteredPairList,KImage,maxSimultaneousPairs);

    rv.selectedPairs.push_back(seedCombination[0]);
}


int main() {
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    std::string filepath = "../images/flow-evolution/39-join.pgm";
    std::string outputFolder = "../output/intExtMerge/fromCurve/";

    boost::filesystem::create_directories(outputFolder);

    Image2D image = DGtal::GenericReader<Image2D>::import(filepath);

    int i =0;
    while(i<100)
    {
        ImageFlowData imf(image);

        imf.init(ImageFlowData::DilationOnly, 5);

        Curve &outerCurve = imf.getMostOuterCurve();
        Curve &innerCurve = imf.getMostInnerCurve();

        RelevantCurve rv;
        rv.curve = innerCurve;

        bool successFlag;
        do {
            CurveCirculator start, end;

            IntervalChecker *intervalChecker = selectMinimizationRegion(start,
                                                                        end,
                                                                        successFlag,
                                                                        imf.getKSpace(),
                                                                        rv);

            if (successFlag) {
                Curve minCurve;
                findOneExpansionMinimumEnergy(minCurve,
                                              intervalChecker,
                                              rv,
                                              imf.getKSpace(),
                                              rv.curve,
                                              imf.getMostOuterCurve());


                CheckableSeedPair& csp = rv.selectedPairs.back();
                SCellCirculator itStart = csp.data().first.firstCirculator;


                ++itStart;
                do
                {
                    rv.not_allowed.insert(*itStart);
                    ++itStart;
                }while(itStart!=csp.data().second.secondCirculator);


                Board2D board;
                for(auto it=rv.not_allowed.begin();it!=rv.not_allowed.end();++it)
                {
                    board << *it;
                    board.saveEPS( (std::to_string(i) + ".eps").c_str() );
                }


                constructImageFromCurve(image, image.domain(), minCurve);
                DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i + 1) + ".pgm", image);
            } else {
                std::cout << "THE END" << std::endl;
            }

            ++i;

        } while (successFlag);

        Curve newCurve;
        for(auto it=rv.curve.begin();it!=rv.curve.end();++it)
        {
            if( rv.not_allowed.find(*it)!=rv.not_allowed.end() ) continue;
            newCurve.push_back(*it);
        }

        for(auto it=rv.selectedPairs.begin();it!=rv.selectedPairs.end();++it)
        {
            newCurve.push_back(it->data().first.connectors[0]);
            SCellCirculator itStart = it->data().first.secondCirculator;

            do
            {
                newCurve.push_back(*itStart);
                ++itStart;
            }while(itStart!=it->data().second.firstCirculator);
            newCurve.push_back(*itStart);

            newCurve.push_back(it->data().second.connectors[0]);
        }

//        assert(newCurve.isValid());
//        Board2D board;
//        board << newCurve;
//        board.saveEPS(("teste.eps"));
//        exit(1);

        constructImageFromCurve(image, image.domain(), newCurve);
        DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i) + "-join.pgm", image);

    }


}