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
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;
    typedef ConnectorSeed<SCell,SCellCirculator,KSpace> SCellConnectorSeed;

    typedef std::pair<SCellConnectorSeed,SCellConnectorSeed> ConnectorSeedPair;

    Curve curve;
    std::set<Curve::SCell> not_allowed;

    std::list<ConnectorSeedPair> selectedPairs;
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


class GoodSCellPair
{
public:
    typedef typename Curve::IncidentPointsRange::ConstCirculator AdapterCirculator;
    typedef StabbingCircleComputer<AdapterCirculator> SegmentComputer;
    typedef CurvatureFromDCAEstimator<SegmentComputer, false> SCFunctor;
    typedef DGtal::SaturatedSegmentation<SegmentComputer> Segmentation;

    typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;

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

    struct PreCandidate
    {
        PreCandidate(SegmentComputer scc):scc(scc)
        {
            SCEstimator estimator;
            estimator.init(0.01,scc.begin(),scc.end());
            estimator.attach(scc);
            curvatureValue = estimator.eval(scc.begin());
        }

        SegmentComputer scc;
        double curvatureValue;
    };

    struct SolutionCompare
    {
        bool operator()(const SolutionCandidate& s1, const SolutionCandidate& s2)
        {
            return s1.value > s2.value;
        }
    };

public:
    GoodSCellPair(const KSpace& KImage, RelevantCurve& rv):KImage(KImage),
                                                     curve(rv.curve)
    {
        Board2D board;
        KSpace K2 = KImage;
        Artist EA(K2,board);
        EA.drawMaximalStabbingCircles(rv.curve);
        EA.board.saveEPS("teste.eps");
//        exit(1);

        Curve::IncidentPointsRange intRange = curve.getIncidentPointsRange();

        Segmentation seg( intRange.c(), intRange.c(), SegmentComputer() );
        Segmentation::SegmentComputerIterator it = seg.begin();


        maxIterations = sizeCirc(seg.begin(),seg.end());
        windowSize = 1+3;
        int i=0;
        for(;i<windowSize;++i)
        {
            q.push_back(PreCandidate(*it));
            ++it;
        }


        findPotentialSolutions(it);
        solutionCandidates.sort(SolutionCompare());

        while(solutionCandidates.size()>10)
        {
            solutionCandidates.pop_back();
        }
    }

    bool next(SolutionCandidate& sc)
    {

        if(solutionCandidates.size()==0) return false;
        sc = solutionCandidates.front();
        solutionCandidates.pop_front();

//        Board2D board;
//        board << sc.scc1;
//        board << sc.scc2;
//        board.save("teste.eps");
//        exit(1);

        return true;
    }

private:
    void findPotentialSolutions(Segmentation::SegmentComputerIterator it)
    {   Board2D board;
        board << DGtal::SetMode(SegmentComputer().className(), "Annulus");
        int currentIteration=0;
        int minimumDet=100;
        do {
            PreCandidate pc = q.front();
            q.pop_front();

            int distance=0;
            for (std::list<PreCandidate>::const_iterator itn = q.begin(); itn != q.end(); ++itn)
            {
                board.clear();
                board << pc.scc;
                board << itn->scc;
                board.saveEPS( (std::to_string(distance) + ".eps").c_str());

                solutionCandidates.push_back(SolutionCandidate(pc.scc,
                                                               itn->scc,
                                                               2*(pc.curvatureValue + itn->curvatureValue) + distance ) );
                distance++;
            }

            q.push_back(PreCandidate(*it));
            ++it;
            ++currentIteration;
        }while(currentIteration < maxIterations);
    }

private:
    const KSpace& KImage;
    Curve& curve;

    int windowSize;
    int maxIterations;
    std::list<PreCandidate> q;
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
        AdapterCirculator incidentPointsSeg2End = sc.scc2.end();

        std::cout << sizeCirc(incidentPointsSeg1Begin, incidentPointsSeg1End) << "::" << sc.value << std::endl;
        std::cout << sizeCirc(incidentPointsSeg2Begin, incidentPointsSeg2End) << std::endl;


        typedef StabbingCircleComputer<AdapterCirculator> SegmentComputer;

        Board2D board;
        board << DGtal::SetMode(SegmentComputer().className(), "Annulus");
        board << sc.scc1;
        board.saveEPS("tempStab.eps");


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

        board.save("interval.eps");

        intervalChecker = new IntervalChecker(start,end,5);
    }

public:
    IntervalChecker* intervalChecker;

};

class OneExpansionMinimumEnergy
{
public:
    typedef SeparateInnerAndOuter::ConnectorSeed ConnectorSeed;
    typedef SeparateInnerAndOuter::ConnectorSeedIterator ConnectorSeedIterator;

    typedef std::vector< CheckableSeedPair::DataType > DataType;
    typedef std::vector< CheckableSeedPair > CheckableSeedList;
    typedef std::vector< CheckableSeedPair > SelectedPairs;

    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    typedef std::pair<ConnectorSeed,ConnectorSeed> ConnectorSeedPair;
    typedef std::set<Curve::SCell> NotAllowedSet;

public:
    OneExpansionMinimumEnergy(KSpace& KImage,
                              Curve& innerCurve,
                              Curve& outerCurve):KImage(KImage),
                                                       innerCurve(innerCurve),
                                                       outerCurve(outerCurve)
    {

    }

    Curve operator()(IntervalChecker* intervalChecker)
    {
        selectedPairs.clear();

        SeparateInnerAndOuter::SeedVector fromInnerSeeds;
        SeparateInnerAndOuter::SeedVector fromOuterSeeds;

        {
            SeparateInnerAndOuter _(KImage,innerCurve,outerCurve);
            _(fromInnerSeeds,fromOuterSeeds);
        }



        DataType pairSeedList;
        ContainerCombinator<ConnectorSeedIterator,ConnectorSeedIterator> combinator;

        combinator.operator()(fromInnerSeeds.begin(),
                              fromInnerSeeds.end(),
                              fromOuterSeeds.begin(),
                              fromOuterSeeds.end(),
                              std::back_inserter(pairSeedList));



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

        Curve minCurve;
        CheckableSeedPair seedCombination[10];
        CE(minCurve,seedCombination,filteredPairList,KImage,maxSimultaneousPairs);

        selectedPairs.push_back(seedCombination[0]);
        return minCurve;
    }

    bool allowedSolution(NotAllowedSet& nas, std::list<ConnectorSeedPair>& updateList)
    {
        CheckableSeedPair& csp = selectedPairs.back();
        SCellCirculator itStart = csp.data().first.firstCirculator;

        ++itStart;
        do
        {
            if(nas.find(*itStart)!=nas.end()) return false;
            ++itStart;
        }while(itStart!=csp.data().second.secondCirculator);

        itStart = csp.data().first.firstCirculator;
        ++itStart;
        do
        {
            nas.insert(*itStart);
            ++itStart;
        }while(itStart!=csp.data().second.secondCirculator);

        updateList.push_back(csp.data());
        return true;
    }

private:
    KSpace& KImage;
    Curve& innerCurve;
    Curve& outerCurve;

    SelectedPairs selectedPairs;

};


bool nextBestCurve(GoodSCellPair& goodPair,
                   ImageFlowData& imf,
                   RelevantCurve& rv,
                   Curve& minCurve)
{
    GoodSCellPair::SolutionCandidate sc;
    OneExpansionMinimumEnergy oeme(imf.getKSpace(),imf.getMostInnerCurve(),imf.getMostOuterCurve());
    while( goodPair.next(sc) )
    {
        PrepareSolutionCandidate psc(imf.getKSpace(),rv.curve,5,sc);
        minCurve = oeme(psc.intervalChecker);

        if( oeme.allowedSolution( rv.not_allowed, rv.selectedPairs  ) ) return true;

    }

    return false;
}

class CurveUpdater
{
public:
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    typedef RelevantCurve::ConnectorSeedPair ConnectorSeedPair;
    typedef std::list< ConnectorSeedPair > SelectedPairs;

    CurveUpdater(Curve& newCurve, SelectedPairs& selectedPairs,int gluedCurveLength)
    {
        selectedPairs.push_back(selectedPairs.front());

        ConnectorSeedPair current;
        while(selectedPairs.size()>1)
        {
            current = selectedPairs.front();
            selectedPairs.pop_front();

            createSCellsFromConnectionSeedPair(newCurve,current,selectedPairs.front(),gluedCurveLength);
        }
        selectedPairs.pop_front();

    }

private:
    void addSCellsInBetween(Curve& newCurve, SCellCirculator itb, SCellCirculator ite)
    {
        do{
            newCurve.push_back(*itb);
            ++itb;
        }while(itb!=ite);
    }

    void createSCellsFromConnectionSeedPair(Curve& newCurve,
                                            RelevantCurve::ConnectorSeedPair& currentCSP,
                                            RelevantCurve::ConnectorSeedPair& nextCSP,
                                            int gluedCurveLength)
    {
        newCurve.push_back(currentCSP.first.connectors[0]);


        SCellCirculator itExternB = currentCSP.first.secondCirculator;
        SCellCirculator itExternE = currentCSP.second.firstCirculator;
        ++itExternE;

        addSCellsInBetween(newCurve,itExternB,itExternE);
        newCurve.push_back(currentCSP.second.connectors[0]);


        SCellCirculator itBridge = currentCSP.second.secondCirculator;
        SCellCirculator itNextB = nextCSP.first.firstCirculator;
        ++itNextB;

        addSCellsInBetween(newCurve,itBridge,itNextB);


    }

};




int main() {
    typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
    typedef DGtal::Circulator<SCellIterator> SCellCirculator;

    std::string filepath = "../images/flow-evolution/10-join.pgm";
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
        GoodSCellPair goodPair(imf.getKSpace(),rv);
        do {
            Curve minCurve;
            successFlag = nextBestCurve(goodPair,imf,rv,minCurve);


            if(successFlag)
            {
//                Board2D board;
//                for(auto it=rv.not_allowed.begin();it!=rv.not_allowed.end();++it)
//                {
//                    board << *it;
//                    board.saveEPS( (std::to_string(i) + ".eps").c_str() );
//                }


                constructImageFromCurve(image, image.domain(), minCurve);
                DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i + 1) + ".pgm", image);
            }

            ++i;
        } while (successFlag);

        Curve newCurve;
        CurveUpdater(newCurve,rv.selectedPairs,5);

//        Board2D board;
//        board << newCurve;
//        board.saveEPS(("teste.eps"));

        assert(newCurve.isValid());


        constructImageFromCurve(image, image.domain(), newCurve);
        DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i) + "-join.pgm", image);

    }


}