#ifndef SEGBYCUT_IMAGEPROC2_H
#define SEGBYCUT_IMAGEPROC2_H

#include <DGtal/shapes/implicit/ImplicitBall.h>
#include <DGtal/shapes/Shapes.h>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"

#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include "FourNeighborhoodPredicate.h"
#include "EightNeighborhoodPredicate.h"


namespace ImageProc2 {
    struct ImageAsDigitalSet {
        typedef DGtal::Z2i::KSpace KImage;
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;

        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

        ImageAsDigitalSet(DigitalSet &digitalSet, std::string imagePath) {
            int threshValue = 100;
            Image2D image = DGtal::GenericReader<Image2D>::import(imagePath);

            DigitalSetInserter inserter(digitalSet);
            DGtal::setFromImage(image, inserter, threshValue, 255);
        }

        ImageAsDigitalSet(DigitalSet &digitalSet, const Image2D &image) {
            int threshValue = 100;

            DigitalSetInserter inserter(digitalSet);
            DGtal::setFromImage(image, inserter, threshValue, 255);
        }

    };

    struct ImageToCVMat {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

        ImageToCVMat(cv::Mat &cvImg, const Image2D &dgtalImg) {
            int ubY = dgtalImg.domain().upperBound()[1];

            for (auto it = dgtalImg.domain().begin(); it != dgtalImg.domain().end(); ++it) {
                Point p = *it;
                unsigned char v(dgtalImg(*it));
                cvImg.at<unsigned char>((ubY - p[1]), p[0]) = v;
            }
        }
    };

    struct CVMatToImage {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

        CVMatToImage(Image2D &dgtalImg, const cv::Mat &cvImg) {
            int ubY = cvImg.rows - 1;
            for (int i = 0; i < cvImg.rows; i++) {
                for (int j = 0; j < cvImg.cols; j++) {
                    unsigned char v(cvImg.at<unsigned char>(i, j));
                    dgtalImg.setValue(Point(j, ubY - i), v);
                }
            }
        }
    };

    struct DigitalSetToCVMat {
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::Z2i::DigitalSet DigitalSet;

        DigitalSetToCVMat(cv::Mat &cvImg, const DigitalSet &dgtalSet) {
            int ubY = dgtalSet.domain().upperBound()[1];

            for (auto it = dgtalSet.begin(); it != dgtalSet.end(); ++it) {
                Point p = *it;
                unsigned char v = (unsigned char) (dgtalSet(*it)) ? 255 : 0;
                cvImg.at<unsigned char>((ubY - p[1]), p[0]) = v;
            }
        }
    };


    struct CVMatToDigitalSet {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

        CVMatToDigitalSet(DigitalSet &dgtalSet, const cv::Mat &cvImg) {
            Domain domain(Point(0, 0), Point(cvImg.cols - 1, cvImg.rows - 1));

            Image2D temp(domain);
            int ubY = cvImg.rows - 1;
            for (int i = 0; i < cvImg.rows; i++) {
                for (int j = 0; j < cvImg.cols; j++) {
                    unsigned char v(cvImg.at<unsigned char>(i, j));
                    temp.setValue(Point(j, ubY - i), v);
                }
            }
            ImageAsDigitalSet(dgtalSet, temp);
        }
    };


    struct Dilate {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

        typedef enum {
            RECT = cv::MORPH_RECT, CROSS = cv::MORPH_CROSS
        } StructuringElement;

        Dilate(Image2D &imgOut, const DigitalSet &dsIn, StructuringElement se, int size = 1) {
            int r = dsIn.domain().upperBound()[1] + 1;
            int c = dsIn.domain().upperBound()[0] + 1;

            cv::Mat cvSrc(r, c, CV_8UC1);
            cv::Mat dilation_dst(r, c, CV_8UC1);
            dilation_dst = 0;
            cvSrc = 0;

            DigitalSetToCVMat(cvSrc, dsIn);


            cv::Mat element = cv::getStructuringElement(se,
                                                        cv::Size(2 * size + 1, 2 * size + 1),
                                                        cv::Point(size, size));

            cv::dilate(cvSrc, dilation_dst, element);
            CVMatToImage(imgOut, dilation_dst);
        }

        Dilate(DigitalSet &dgtalSetOut, const DigitalSet &dsIn, StructuringElement se, int size = 1) {
            int r = dsIn.domain().upperBound()[1] + 1;
            int c = dsIn.domain().upperBound()[0] + 1;

            cv::Mat cvSrc(r, c, CV_8UC1);
            cv::Mat dilation_dst(r, c, CV_8UC1);
            dilation_dst = 0;
            cvSrc = 0;

            DigitalSetToCVMat(cvSrc, dsIn);


            cv::Mat element = cv::getStructuringElement(se,
                                                        cv::Size(2 * size + 1, 2 * size + 1),
                                                        cv::Point(size, size));

            cv::dilate(cvSrc, dilation_dst, element);
            CVMatToDigitalSet(dgtalSetOut, dilation_dst);
        }

    };

    struct Erode {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

        typedef enum {
            RECT = cv::MORPH_RECT, CROSS = cv::MORPH_CROSS
        } StructuringElement;

        Erode(Image2D &imgOut, const DigitalSet &dsIn, StructuringElement se, int size = 1) {
            int r = dsIn.domain().upperBound()[1] + 1;
            int c = dsIn.domain().upperBound()[0] + 1;

            cv::Mat cvSrc(r, c, CV_8UC1);
            cv::Mat dilation_dst(r, c, CV_8UC1);
            dilation_dst = 0;
            cvSrc = 0;

            DigitalSetToCVMat(cvSrc, dsIn);


            cv::Mat element = cv::getStructuringElement(se,
                                                        cv::Size(2 * size + 1, 2 * size + 1),
                                                        cv::Point(size, size));

            cv::erode(cvSrc, dilation_dst, element);
            CVMatToImage(imgOut, dilation_dst);
        }

        Erode(DigitalSet &dgtalSetOut, const DigitalSet &dsIn, StructuringElement se, int size = 1) {
            int r = dsIn.domain().upperBound()[1] + 1;
            int c = dsIn.domain().upperBound()[0] + 1;

            cv::Mat cvSrc(r, c, CV_8UC1);
            cv::Mat dilation_dst(r, c, CV_8UC1);
            dilation_dst = 0;
            cvSrc = 0;

            DigitalSetToCVMat(cvSrc, dsIn);


            cv::Mat element = cv::getStructuringElement(se,
                                                        cv::Size(2 * size + 1, 2 * size + 1),
                                                        cv::Point(size, size));

            cv::erode(cvSrc, dilation_dst, element);
            CVMatToDigitalSet(dgtalSetOut, dilation_dst);
        }

    };

    template<typename TNeighborhood>
    struct DigitalBoundary {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::Z2i::DigitalSet DigitalSet;

        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
        typedef TNeighborhood NeighborhoodPredicate;


        typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

        DigitalBoundary(DigitalSet &boundaryDS, DigitalSet &originalDS) {
            NeighborhoodPredicate NP(originalDS);

            DigitalSetInserter inserter(boundaryDS);

            std::remove_copy_if(originalDS.begin(), originalDS.end(), inserter, NP);
        }
    };

    struct DigitalSetExplorer {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::Z2i::DigitalSet DigitalSet;

        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
        typedef FourNeighborhoodPredicate<DigitalSet> NeighborhoodPredicate;


        typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

        DigitalSetExplorer(DigitalSet &digitalSetOut, const DigitalSet &digitalSetIn, Point startingPoint,
                           bool stopWhenIn = true) {
            Domain d = digitalSetIn.domain();
            Point filter[4] = {Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1)};

            Point lowerBound = digitalSetIn.domain().lowerBound();
            Point upperBound = digitalSetIn.domain().upperBound();

            std::stack<Point> stack;
            std::set<Point> markedSet;
            stack.push(startingPoint);

            Point p;
            while (!stack.empty()) {
                p = stack.top();
                stack.pop();
                if (markedSet.find(p) != markedSet.end()) continue;
                markedSet.insert(p);

                digitalSetOut.insert(p);

                for (int i = 0; i < 4; ++i) {
                    Point np = p + filter[i];
                    if (np[0] < lowerBound[0] || np[1] < lowerBound[1]) continue;
                    if (np[0] > upperBound[0] || np[1] > upperBound[1]) continue;

                    if (stopWhenIn) {
                        if (digitalSetIn(np)) continue;
                    } else {
                        if (!digitalSetIn(np)) continue;
                    }


                    stack.push(np);
                }
            }


        }
    };

    struct NoHoles {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;

        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

        NoHoles(DigitalSet &noHoles, const DigitalSet &dsIn) {
            Domain domain = dsIn.domain();

            DigitalSet exploredSet(domain);
            DigitalSetExplorer(exploredSet, dsIn, Point(0, 0), true);

            DigitalSetInserter diNoHoles(noHoles);
            exploredSet.computeComplement(diNoHoles);
        }
    };

    struct SetDifference {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;

        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

        SetDifference(DigitalSet &differenceSet, const DigitalSet &A, const DigitalSet &B) {
            Domain domain = A.domain();

            DigitalSet temp(domain);
            DigitalSetInserter diTemp(temp);
            A.computeComplement(diTemp);
            temp += B;

            DigitalSetInserter diDiff(differenceSet);
            temp.computeComplement(diDiff);
        }
    };

    struct ThickBorder {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;

        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

        ThickBorder(DigitalSet &thickBorder, const DigitalSet &dsIn, int thickness) {
            Domain domain = dsIn.domain();

            DigitalSet eroded(domain);
            Erode(eroded, dsIn, Erode::RECT, thickness);

            SetDifference(thickBorder, dsIn, eroded);
        }
    };

    struct DigitalIntersection {
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::DigitalSet DigitalSet;

        DigitalIntersection(DigitalSet &digitalIntersection, const DigitalSet &A, const DigitalSet &B) {
            Domain domain = A.domain();

            for (auto it = A.begin(); it != A.end(); ++it) {
                if (B(*it)) digitalIntersection.insert(*it);
            }
        }
    };

    class DigitalBallIntersection {
    public:
        typedef unsigned int Radius;

        typedef DGtal::Z2i::Space Space;
        typedef DGtal::Z2i::Point Point;
        typedef DGtal::Z2i::Domain Domain;


        typedef DGtal::Z2i::DigitalSet DigitalSet;
        typedef DGtal::ImplicitBall<Space> EuclideanBall;

    public:
        DigitalBallIntersection(Radius r, const DigitalSet &intersectWith) : _r(r), _ds(intersectWith) {}

        void operator()(DigitalSet &intersectionSet, Point center) {
            DigitalSet db(intersectionSet.domain());
            digitalBall(db, center, _r);
            DigitalIntersection(intersectionSet, db, _ds);
        }

        void digitalBall(DigitalSet &db, Point center, int radius) {
            EuclideanBall eb(center, _r);
            DGtal::Shapes<Domain>::euclideanShaper(db, eb, 1);
        }

    private:
        Radius _r;
        const DigitalSet &_ds;
    };

}

#endif //SEGBYCUT_IMAGEPROC2_H
