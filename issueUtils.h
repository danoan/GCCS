#ifndef SEGBYCUT_ISSUEUTILS_H
#define SEGBYCUT_ISSUEUTILS_H

#include "boost/filesystem.hpp"

#include <DGtal/io/boards/Board2D.h>
#include "DGtal/geometry/curves/FP.h"
#include <boost/filesystem/path.hpp>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "utils.h"
#include "typeDefs.h"

namespace Issue {

    using namespace DGtal;
    using namespace DGtal::Z2i;

    using namespace SegCut;

    template<typename Iterator, typename Board>
    void drawCCP(const Iterator &itb, const Iterator &ite, Board &aBoard) {

        //choose the drawing mode
        aBoard << SetMode("ArithmeticalDSS", "BoundingBox");
        //prepare the drawing style and the pen color
        std::string aStyleName = "ArithmeticalDSS/BoundingBox";
        CustomPenColor *aPenColor;

        //for each maximal segment
        for (Iterator i(itb); i != ite; ++i) {

            //get the current maximal segment
            typedef typename Iterator::SegmentComputer::Primitive DSS;
            DSS maximalDSS = i->primitive();

            //if located at the end of a connected part
            if (!(i.intersectNext() && i.intersectPrevious())) {

                aPenColor = new CustomPenColor(DGtal::Color::Black);

                //otherwise
            } else {

                //get the points located before and after the maximal segment
                typedef typename DSS::Point Point;
                Point beforeFirst = *(--(i->begin()));
                Point afterLast = *(i->end());

                //remainders and bounds
                typedef typename DSS::Integer Integer;
                Integer r1 = maximalDSS.remainder(beforeFirst);
                Integer r2 = maximalDSS.remainder(afterLast);
                Integer mu = maximalDSS.mu();
                Integer omega = maximalDSS.omega();

                //configurations
                if ((r1 <= mu - 1) && (r2 <= mu - 1)) {                    //concave
                    aPenColor = new CustomPenColor(DGtal::Color::Green);
                } else if ((r1 >= mu + omega) && (r2 >= mu + omega)) {     //convex
                    aPenColor = new CustomPenColor(DGtal::Color::Blue);
                } else if ((r1 >= mu + omega) && (r2 <= mu - 1)) {         //convex to concave
                    aPenColor = new CustomPenColor(DGtal::Color::Yellow);
                } else if ((r1 <= mu - 1) && (r2 >= mu + omega)) {         //concave to convex
                    aPenColor = new CustomPenColor(DGtal::Color::Yellow);
                } else {                                           //pb
                    aPenColor = new CustomPenColor(DGtal::Color::Red);
                }

            }

            // draw the maximal segment on the board
            aBoard << CustomStyle(aStyleName, aPenColor)
                   << maximalDSS;

        }

    }


    template<typename Iterator, typename Board>
    void segmentationIntoMaximalDSSs(const Iterator &itb, const Iterator &ite,
                                     Board &aBoard) {
        typedef typename Issue::IteratorCirculatorTraits<Iterator>::Value::Coordinate Coordinate;

        //choose the primitive computer and the segmentation
        typedef Issue::ArithmeticalDSSComputer<Iterator, Coordinate, 4> RecognitionAlgorithm;
        typedef Issue::SaturatedSegmentation<RecognitionAlgorithm> Segmentation;

        //create the segmentation
        RecognitionAlgorithm algo;
        Segmentation s(itb, ite, algo);

        //draw the result
        drawCCP(s.begin(), s.end(), aBoard);
    }

    void drawFaithfulPolygon(Issue::Image2D &image,
                             std::string outputFolder) {
        Board2D board;
        Curve intCurve, extCurve;
        KSpace KImage;

        computeBoundaryCurve(intCurve, KImage, image, 100);

        GridCurve<KSpace>::PointsRange range = intCurve.getPointsRange();

        FP<GridCurve<KSpace>::PointsRange::ConstCirculator, Integer, 4> fp(range.c(), range.c());
        board << fp;

        segmentationIntoMaximalDSSs(range.c(), range.c(), board);


        std::string fpOutputFilePath = outputFolder + "/faithful-polygon";
        board.saveEPS(fpOutputFilePath.c_str());
    }

    template<typename RangeIterator>
    void howManyInflections(const RangeIterator &pointRangeBegin, const RangeIterator &pointRangeEnd, int &ninflections,
                            int &nconcavities) {
        typedef typename IteratorCirculatorTraits<RangeIterator>::Value::Coordinate Coordinate;

        //choose the primitive computer and the segmentation
        typedef ArithmeticalDSSComputer<RangeIterator, Coordinate, 4> RecognitionAlgorithm;
        typedef SaturatedSegmentation<RecognitionAlgorithm> Segmentation;

        typedef typename SaturatedSegmentation<RecognitionAlgorithm>::SegmentComputerIterator Iterator;

        //create the segmentation
        RecognitionAlgorithm algo;
        Segmentation s(pointRangeBegin, pointRangeEnd, algo);

        Iterator itb = s.begin();
        Iterator ite = s.end();


        ninflections = 0;
        nconcavities = 0;

        for (Iterator i(itb); i != ite; ++i) {

            //get the current maximal segment
            typedef typename Iterator::SegmentComputer::Primitive DSS;
            DSS maximalDSS = i->primitive();

            //if located at the end of a connected part
            if (!(i.intersectNext() && i.intersectPrevious())) {

                //otherwise
            } else {

                //get the points located before and after the maximal segment
                typedef typename DSS::Point Point;
                Point beforeFirst = *(--(i->begin()));
                Point afterLast = *(i->end());

                //remainders and bounds
                typedef typename DSS::Integer Integer;
                Integer r1 = maximalDSS.remainder(beforeFirst);
                Integer r2 = maximalDSS.remainder(afterLast);
                Integer mu = maximalDSS.mu();
                Integer omega = maximalDSS.omega();

                //configurations
                if ((r1 <= mu - 1) && (r2 <= mu - 1)) {                    //concave
                    ++nconcavities;
                } else if ((r1 >= mu + omega) && (r2 >= mu + omega)) {     //convex

                } else if ((r1 >= mu + omega) && (r2 <= mu - 1)) {         //convex to concave
                    ++ninflections;
                } else if ((r1 <= mu - 1) && (r2 >= mu + omega)) {         //concave to convex
                    ++ninflections;
                } else {                                           //pb

                }

            }


        }

    }

    void updateGluedWeightUsingFP(ConnectorSeedRangeType &seedRange,
                                  KSpace &KImage,
                                  std::map<Z2i::SCell, double> &weightMap,
                                  double cmin,
                                  double cmax)
    {
        typedef ConstRangeAdapter<SCellGluedCurveIterator,
                DGtal::functors::SCellToPoint<KSpace>,
                DGtal::Z2i::KSpace::Point> GluedCurvePointsRange;

        DGtal::functors::SCellToPoint<KSpace> myFunc(KImage);

        SeedToGluedCurveRangeFunctor stgcF(10);
        GluedCurveSetRange gcsRange(seedRange.begin(),
                                    seedRange.end(),
                                    stgcF);

        int ninflections, nconcavities;
        double f = 0.5*cmax;
        for (GluedCurveIteratorPair it = gcsRange.begin(); it != gcsRange.end(); ++it) {
            SCellGluedCurveIterator gcBegin = it->first;
            SCellGluedCurveIterator gcEnd = it->second;


            GluedCurvePointsRange range(gcBegin, gcEnd, myFunc);

            howManyInflections(range.begin(), range.end(), ninflections, nconcavities);
            std::cout << ninflections << "::" << nconcavities << std::endl;

            weightMap[gcBegin.linkSurfel()] += f*nconcavities;


        };
    };

}

#endif //SEGBYCUT_ISSUEUTILS_H
