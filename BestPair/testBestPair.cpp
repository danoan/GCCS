#include "../FlowGraph/ImageFlowData.h"
#include "../Flows/PreprocessImage.h"
#include "../BestPair/API.h"
#include "../FlowGraph/FlowGraphUtils.h"
#include "JointSet.h"

namespace Development{
    bool solveShift = false;
    bool crossElement = false;

    bool lambdaEstimator = false;
    bool pessimistEstimator = false;

    bool makeConvexArcs = false;
    bool invertGluedArcs = false;

    bool iteractive = false;
    std::string windowName = "simpleFlow";
};

void mapGDets(BestLinel::WeightMap& mapOfGdets, ImageFlowData& ifd, int g)
{
    BestLinel::WeightMap weightMap;
    setArcsWeight(ifd,weightMap);

    ImageFlowData::CurvePairIterator it = ifd.curvePairBegin();
    for(;it!=ifd.curvePairEnd();++it)
    {
        ImageFlowData::GluedCurveSetRange gcsRange = it->getGCSRange();

        for(auto it=gcsRange.begin();it!=gcsRange.end();++it)
        {
            mapOfGdets[*it->first.connectorsBegin()] = BestLinel::gdet( it,weightMap,ifd.getKSpace(),g);
        }
    }
}

void testGDets(ImageData& ID)
{
    int g = 5;
    ImageFlowData ifd(ID.originalImage,ID.originalImage);
    ifd.init(ImageFlowData::FlowMode::DilationOnly,g);

    BestLinel::WeightMap mapOfGdets;

    mapGDets(mapOfGdets,ifd,g);
}

void testJointSet(ImageData& ID)
{
    int g = 5;
    ImageFlowData ifd(ID.originalImage,ID.originalImage);
    ifd.init(ImageFlowData::FlowMode::DilationOnly,g);

    BestLinel::JointSet js(ifd);

    std::function<void(BestLinel::JointSet::JointIterator,
                       BestLinel::JointSet::JointIterator,
                       std::string)> drawIntercalate( [&](BestLinel::JointSet::JointIterator begin,
                                                         BestLinel::JointSet::JointIterator end,
                                                         std::string outputName){
        DGtal::Board2D board;
        auto it = begin;
        DGtal::Color colorList[3] = {DGtal::Color::Blue, DGtal::Color::Yellow, DGtal::Color::Green};
        for (int i = 0; it != end; ++it, ++i) {
            board << DGtal::CustomStyle(it->scell.className(),
                                        new DGtal::CustomColors(colorList[i % 3], colorList[i % 3]));
            board << it->scell;
        }

        board.saveEPS(outputName.c_str());
    });

    drawIntercalate(js.inJointsVector.begin(),js.inJointsVector.end(),"inJoints.eps");
    drawIntercalate(js.outJointsVector.begin(),js.outJointsVector.end(),"outJoints.eps");
}

void addUntil(Curve& c, BestLinel::JointSet::CurveCirculator ccirc, BestLinel::JointSet::SCell until)
{
    while(*ccirc!=until)
    {
        c.push_back(*ccirc);
        ++ccirc;
    }
    c.push_back(*ccirc);
}

void constructCurveFromJoints(DGtal::Z2i::Curve& curve,
                              ImageFlowData& imf,
                              BestLinel::JointSet::Joint inJoint,
                              BestLinel::JointSet::Joint outJoint)
{
    curve.push_back(inJoint.scell);
    addUntil(curve,inJoint.sucCurve,*outJoint.antCurve);

    curve.push_back(outJoint.scell);
    addUntil(curve,outJoint.sucCurve,*inJoint.antCurve);



}

struct Determinant
{
    Determinant(BestLinel::JointSet::Joint j1, BestLinel::JointSet::Joint j2, double sigma):j1(j1),j2(j2),sigma(sigma){}
    BestLinel::JointSet::Joint j1;
    BestLinel::JointSet::Joint j2;
    double sigma;
};

typedef BestLinel::JointSet::JointIterator JointIterator;
typedef DGtal::Circulator<JointIterator> JointCirculator;

void chooseBestPair(int g,
                    std::vector<Determinant>& detList,
                    JointCirculator firstCirc,
                    JointCirculator scndCirc,
                    BestLinel::WeightMap& mapOfGdets)
{
    JointCirculator cinit = firstCirc;
    JointCirculator coutit = walkCirculator(scndCirc,g);//scndCirc;//walkCirculator(joutcirc,2);
    do
    {
        JointCirculator tempOut = coutit;
        double best = mapOfGdets[tempOut->scell];
        BestLinel::JointSet::Joint bestJoint = *tempOut;
        for(int i=0;i<g;++i)
        {
            if(mapOfGdets[tempOut->scell]>best)
            {
                best = mapOfGdets[tempOut->scell];
                bestJoint = *tempOut;
            }
            ++tempOut;
        }
        detList.push_back( Determinant(*cinit,bestJoint,best+mapOfGdets[cinit->scell]) );
        ++cinit;
        ++coutit;
    }while(cinit!=firstCirc);


}


void testBestPair(ImageData& ID)
{
    int g = 5;
    ImageFlowData ifd(ID.originalImage,ID.originalImage);
    ifd.init(ImageFlowData::FlowMode::DilationOnly,g);

    BestLinel::JointSet js(ifd);
    BestLinel::WeightMap mapOfGdets;

    mapGDets(mapOfGdets,ifd,g);

    JointCirculator jincirc(js.inJointsVector.begin(),
                            js.inJointsVector.begin(),
                            js.inJointsVector.end());

    JointCirculator joutcirc(js.outJointsVector.begin(),
                             js.outJointsVector.begin(),
                             js.outJointsVector.end());



    std::vector< Determinant > detListFromInner;
    std::vector< Determinant > detListFromOuter;

    chooseBestPair(g,detListFromInner,jincirc,joutcirc,mapOfGdets);
    chooseBestPair(g,detListFromOuter,joutcirc,jincirc,mapOfGdets);

    std::sort(detListFromInner.begin(),detListFromInner.end(),[&](const struct Determinant& d1, const struct Determinant& d2){return d1.sigma > d2.sigma;});
    std::sort(detListFromOuter.begin(),detListFromOuter.end(),[&](const struct Determinant& d1, const struct Determinant& d2){return d1.sigma > d2.sigma;});

    int i=0;
    for(auto it=detListFromInner.begin();it!=detListFromInner.end();++it,++i)
    {
        std::cout << it->sigma << std::endl;

        Curve ncurve;
        constructCurveFromJoints(ncurve,ifd,detListFromInner[i].j1,detListFromInner[i].j2);
        DGtal::Board2D board;
        board << ncurve;
        board.saveEPS( ("bestCurve" + std::to_string(i) +".eps").c_str() );
    }


}

int main()
{
    std::string imagePath = "../images/segSet/single_square.pgm";
    ImageData ID(imagePath);

    testGDets(ID);
    testJointSet(ID);
    testBestPair(ID);

    return 0;
}