
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
#include "../FreezeOneExp/CombinationsEvaluator.h"

typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;

typedef std::vector<FreezeOneExp::JointSolution> JointSolutionVector;

void findOneExpansionMinimumEnergy(JointSolutionVector& jsv,
                                   KSpace& KImage,
                                   Curve& innerCurve,
                                   Curve& outerCurve,
                                   std::map<SCell,int>& jointOrder,
                                   int maxJoints
)
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
        int connectorsDistance = (jointOrder[it->second.connectors[0]] - jointOrder[it->first.connectors[0]] + maxJoints)%maxJoints;

        //Internal Connector must be different from External Connector
//        if(coordInt==coordExt) continue;

        //Internal and External Connector must cover between 2 and at most 10 external linels
        if( connectorsDistance<=1 || connectorsDistance >= 6) continue;

        filteredPairList.push_back(*it);

    }

    int maxSimultaneousPairs = fromInnerSeeds.size()/4;
    FreezeOneExp::CombinationsEvaluator<1> CE;
    CE.addChecker( new GluedIntersectionChecker() );
    CE.addChecker( new MinimumDistanceChecker(KImage) );
    CE(jsv,filteredPairList,KImage,maxSimultaneousPairs);

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

typedef DGtal::Z2i::Curve::ConstIterator SCellIterator;
typedef DGtal::Circulator<SCellIterator> SCellCirculator;

void addIntervalSCells(std::vector<KSpace::SCell> &vectorOfSCells,
                       SCellCirculator begin,
                       SCellCirculator end) {
    SCellCirculator it = begin;
    while (it != end) {
        vectorOfSCells.push_back(*it);
        ++it;
    }
    vectorOfSCells.push_back(*it);
}

void addSeedPairSCells(std::vector<KSpace::SCell> &vectorOfSCells,
                       CheckableSeedPair &currentPair,
                       CheckableSeedPair &nextPair) {

    typedef ConnectorSeedRange<KSpace, SCellCirculator> MyConnectorSeedRange;
    typedef MyConnectorSeedRange::ConnectorSeedType ConnectorSeed;

    ConnectorSeed inToExtSeed = currentPair.data().first;
    ConnectorSeed extToIntSeed = currentPair.data().second;

    if (currentPair.data().first.cType != ConnectorType::internToExtern) {
        std::cout << "ERROR" << std::endl;
    }

    if (currentPair.data().second.cType != ConnectorType::externToIntern) {
        std::cout << "ERROR" << std::endl;
    }


    vectorOfSCells.push_back(inToExtSeed.connectors[0]);
    addIntervalSCells(vectorOfSCells, inToExtSeed.secondCirculator, extToIntSeed.firstCirculator);


    ConnectorSeed nextIntToExtSeed = nextPair.data().first;
    vectorOfSCells.push_back(extToIntSeed.connectors[0]);
    addIntervalSCells(vectorOfSCells, extToIntSeed.secondCirculator, nextIntToExtSeed.firstCirculator);
}

void createCurve(Curve &curve,
                 CheckableSeedPair *seedPairs,
                 int totalPairs) {
    typedef KSpace::SCell SCell;
    std::vector<SCell> scells;

    CheckableSeedPair currentPair;
    CheckableSeedPair nextPair;
    for (int i = 0; i < totalPairs; ++i) {
        if (i == (totalPairs - 1)) {
            currentPair = seedPairs[totalPairs - 1];
            nextPair = seedPairs[0];
        } else {
            currentPair = seedPairs[i];
            nextPair = seedPairs[i + 1];
        }

        addSeedPairSCells(scells, currentPair, nextPair);
    }

    curve.initFromSCellsVector(scells);
}

void createJointOrder(std::map<SCell,int>& jointOrder, std::vector<SCell>& jointsVector,ImageFlowData& ifd)
{
    std::map<SCell,int> indexOrder;
    int index=0;
    for(auto it=ifd.getMostOuterCurve().begin();it!=ifd.getMostOuterCurve().end();++it,++index)
    {
        indexOrder[*it]=index;
    }

    auto it = ifd.curvePairBegin();

    ImageFlowData::GluedCurveSetRange gcsRange = it->getGCSRange();
    for(auto iu=gcsRange.begin();iu!=gcsRange.end();++iu)
    {


        if(iu->first.connectorType()==ConnectorType::internToExtern)
        {
            jointOrder[*iu->first.connectorsBegin()] = indexOrder[ *iu->first.curveSegment2Begin() ];
        }
        else
        {
            jointOrder[*iu->first.connectorsBegin()] = indexOrder[ *iu->first.curveSegment1End() ];
        }

        jointsVector.push_back(*iu->first.connectorsBegin());


    }
}

template<typename IteratorType>
bool nextSolution(IteratorType& start, IteratorType end, int maxJoints, std::map<SCell,int>& jointOrder ,std::map<int,bool>& validInnerJoints,std::map<int,bool>& validOuterJoints,int& tries)
{
    //it: Iterator to CheckableSeedPair
    auto &it = start;

    int innerIndex;
    int outerIndex;

    if(tries>=3) return false;

    while(true)
    {
        if(it==end) break;

        innerIndex = jointOrder[ it->csp.data().first.connectors[0] ];
        outerIndex = jointOrder[ it->csp.data().second.connectors[0] ];

        bool notGood=true;
        if(validInnerJoints[innerIndex]==true && validOuterJoints[outerIndex]==true)
        {
            notGood=false;
            for(int i = innerIndex;i!=outerIndex;)
            {
                if(validInnerJoints[i]==false || validOuterJoints[i]==false){
                    notGood=true;
                    break;
                }
                ++i;
                i=i%maxJoints;
            }
        }
        if(notGood){++it;++tries;}
        else break;
    }

    if(it==end) return false;



    while(innerIndex!=outerIndex)
    {
        validInnerJoints[innerIndex] = false;
        validOuterJoints[innerIndex] = false;
        innerIndex = (innerIndex+1)%maxJoints;
    }

    validInnerJoints[outerIndex] = false;
    validOuterJoints[outerIndex] = false;

    return true;

}

void computeOneExpansions()
{
    std::string filepath = "../images/flow-evolution/single_square.pgm";
    std::string outputFolder = "../output/freeze-one-exp/fromCurve/";

    boost::filesystem::create_directories(outputFolder);

    Image2D image = DGtal::GenericReader<Image2D>::import(filepath);

    JointSolutionVector jsv;

    int i=0;
    while(i<100)
    {
        jsv.clear();

        ImageFlowData imf(image);
        imf.init(ImageFlowData::DilationOnly, 5);


        std::map<SCell,int> jointOrder;
        std::vector<SCell> jointsVector;
        createJointOrder(jointOrder,jointsVector,imf);
        int maxJoints = imf.getMostOuterCurve().size();

        findOneExpansionMinimumEnergy(jsv,
                                      imf.getKSpace(),
                                      imf.getMostInnerCurve(),
                                      imf.getMostOuterCurve(),
                                      jointOrder,
                                      maxJoints);

        auto it = jsv.begin();
        auto end = jsv.end();
        std::vector<CheckableSeedPair> selectedPairs;


        std::map<int,bool> validInnerJoints;
        std::map<int,bool> validOuterJoints;
        for(int i=0;i<maxJoints;++i)
        {
            validInnerJoints[i] = true;
            validOuterJoints[i] = true;
        }

        {
            std::cout << "Joints Order" << std::endl;
            int i =0;
            for (auto it =jsv.begin();it!=jsv.end();++it,++i)
            {
                std::cout << i << ":" << jointOrder[it->csp.data().first.connectors[0]]
                          << " and "
                          << jointOrder[it->csp.data().second.connectors[0]]
                          << std::endl;
            }
        }

        int tries=0;
        while(nextSolution(it,end,maxJoints,jointOrder,validInnerJoints,validOuterJoints,tries))
        {
            selectedPairs.push_back(it->csp);
        }

        {

            Board2D board;
            DGtal::Color colorList[4] = {DGtal::Color::Blue, DGtal::Color::Yellow, DGtal::Color::Green, DGtal::Color::Red};
            int j=0;
            for(auto it=jointsVector.begin();it!=jointsVector.end();++it)
            {
                if(validInnerJoints[jointOrder[*it]] && validOuterJoints[jointOrder[*it]])
                {
                    board << DGtal::CustomStyle(it->className(),
                                                new DGtal::CustomColors(DGtal::Color::Gray,DGtal::Color::Gray));
                    board << *it;
                }
                else
                {
                    board << DGtal::CustomStyle(it->className(),
                                                new DGtal::CustomColors(DGtal::Color::Black, DGtal::Color::Black));
                    board << *it;
                }

            }
            for(auto it=selectedPairs.begin();it!=selectedPairs.end();++it,++j)
            {
                board << DGtal::CustomStyle(it->data().first.connectors[0].className(),
                                            new DGtal::CustomColors(colorList[j % 4], colorList[j % 4]));
                board << it->data().first.connectors[0];
                board << it->data().second.connectors[0];
                std::cout << it->data().first.connectors[0].preCell().coordinates << " and "
                          << it->data().second.connectors[0].preCell().coordinates << std::endl;
            }
            board.saveEPS( (outputFolder + "selectedPairs" + std::to_string(i+1) + ".eps").c_str() );
        }

        Curve curve;
        createCurve(curve,selectedPairs.data(),selectedPairs.size());

        Image2D tempImage(image.domain());
        constructImageFromCurve(tempImage,tempImage.domain(),curve);
        image = tempImage;
        DGtal::GenericWriter<Image2D>::exportFile(outputFolder + std::to_string(i+1) + ".pgm",image);
        ++i;
    }
}

void computeStabbingCircles()
{
    std::string outputFolder = "../output/freeze-one-exp/fromCurve/";

    boost::filesystem::path p = outputFolder;
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