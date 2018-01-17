#ifndef SEGBYCUT_MOSTCENTEREDMAXIMALSEGMENTESTIMATORWRAPPER_H
#define SEGBYCUT_MOSTCENTEREDMAXIMALSEGMENTESTIMATORWRAPPER_H

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"


#include "MostCenteredMaximalSegmentEstimatorWrapper.h"

namespace Patch{

    template<typename IteratorType>
    IteratorType getMiddleIterator(const IteratorType& itb, const IteratorType& ite, bool& flag)
    {
        IteratorType b( itb );
        IteratorType f( ite );
        flag = true;
        while (b != f) {
            if (flag) {
                --f;
                flag = false;
            } else {
                ++b;
                flag = true;
            }
        }
        return b;
    }


    template<typename SegmentComputer, typename SCEstimator>
    class MostCenteredMaximalSegmentEstimator{

    public:
        typedef DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MyMCMSE;

        typedef typename MyMCMSE::ConstIterator ConstIterator;
        typedef typename MyMCMSE::Quantity Quantity;

        typedef DGtal::SaturatedSegmentation<SegmentComputer> Segmentation;
        typedef typename Segmentation::SegmentComputerIterator SegmentIterator;

        typedef typename std::back_insert_iterator< std::vector< Quantity > > OutputIterator;

        SegmentComputer mySC;
        SCEstimator mySCEstimator;
        double myH;
        ConstIterator myBegin,myEnd;


        MostCenteredMaximalSegmentEstimator(const MostCenteredMaximalSegmentEstimator& other)
        {
            myH = other.myH;
            mySC = other.mySC;
            mySCEstimator = other.mySC;
        }

        MostCenteredMaximalSegmentEstimator(const SegmentComputer& aSegmentComputer,
                                            const SCEstimator& aSCEstimator)
        {
            myH = 0;
            mySC = aSegmentComputer;
            mySCEstimator = aSCEstimator;
        }


        inline
        bool isValid()
        {
            return ( (myH > 0)&&(isNotEmpty(myBegin, myEnd)) );
        }


        inline void init(const double h, const ConstIterator& itb, const ConstIterator& ite)
        {

            myH = h;
            myBegin = itb;
            myEnd = ite;

            if (isValid())
                mySCEstimator.init( myH, myBegin, myEnd );

        }


        inline OutputIterator endEval(const ConstIterator& itb,
                                      const ConstIterator& ite,
                                      ConstIterator& itCurrent,
                                      SegmentIterator& first,
                                      SegmentIterator& last,
                                      OutputIterator result,
                                      DGtal::CirculatorType )
        {
            if ( (itb == ite) && (first.intersectPrevious() && last.intersectNext() ) )
            {//if first and last segment intersect (whole range)
                //last segment
                bool full_cycle;
                auto sEnd(last->end());
                ConstIterator itEnd = Patch::getMiddleIterator( first->begin(), sEnd, full_cycle );//(floor)


                mySCEstimator.attach( *last );
                result = mySCEstimator.eval( itCurrent, itEnd, result );

                itCurrent = itEnd;

                if (itCurrent != ite)
                {
                    //first segment
                    mySCEstimator.attach( *first );
                    result = mySCEstimator.eval( itCurrent, ite, result );
                }
            }
            else
            { //(sub range)
                mySCEstimator.attach( *last );
                result = mySCEstimator.eval( itCurrent, ite, result );
            }
            return result;
        }

        inline
        OutputIterator endEval(const ConstIterator& /*itb*/,
                               const ConstIterator& ite,
                               ConstIterator& itCurrent,
                               SegmentIterator& /*first*/,
                               SegmentIterator& last,
                               OutputIterator result,
                               DGtal::IteratorType )
        {
            mySCEstimator.attach( *last );
            result = mySCEstimator.eval( itCurrent, ite, result );
            return result;
        }

        inline OutputIterator endEval(const ConstIterator& itb,
                                      const ConstIterator& ite,
                                      ConstIterator& itCurrent,
                                      SegmentIterator& first,
                                      SegmentIterator& last,
                                      OutputIterator result)
        {
            typedef typename DGtal::IteratorCirculatorTraits<ConstIterator>::Type Type;
            return endEval (itb, ite, itCurrent, first, last, result, Type() );
        }

        inline OutputIterator
        eval(const ConstIterator& itb,
             const ConstIterator& ite,
             OutputIterator result)
        {


            Segmentation seg(myBegin, myEnd, mySC);
            seg.setSubRange(itb, ite);
            if ((myBegin != itb) || (myEnd != ite))
            { //if subrange
                seg.setMode("MostCentered++");
            }
            else
            {//whole range
                seg.setMode("MostCentered");
            }

            if (isValid()) {

                SegmentIterator segItBegin = seg.begin();
                SegmentIterator segItEnd = seg.end();
                SegmentIterator segIt = segItBegin;
                SegmentIterator nextSegIt = segIt;

                if (nextSegIt != segItEnd )
                {  //at least one maximal segment
                    ++nextSegIt;

                    if (nextSegIt == segItEnd )
                    {    //only one maximal segment
                        mySCEstimator.attach( *segIt );
                        result = mySCEstimator.eval( itb, ite, result );
                    }
                    else
                    {   //strictly more than one maximal segment
                        ConstIterator itCurrent = itb;

                        //main loop
                        while (nextSegIt != segItEnd)
                        {
                            bool full_cycle;
                            auto sEnd(segIt->end());
                            ConstIterator itEnd = Patch::getMiddleIterator( nextSegIt->begin(), sEnd, full_cycle );//(floor)


                            mySCEstimator.attach( *segIt );
                            result = mySCEstimator.eval( itCurrent, itEnd, result );

                            itCurrent = itEnd;

                            segIt = nextSegIt;
                            ++nextSegIt;
                        }

                        //end
                        result = endEval(itb, ite, itCurrent, segItBegin, segIt, result);

                    }//end one or more maximal segments test
                }//end zero or one maximal segment test
                return result;

            }
            else
            {//nothing is done without correct initialization
                std::cerr << "[DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& itb, const ConstIterator& ite,OutputIterator result)]"
                          << " ERROR. Object is not initialized." << std::endl;
                throw DGtal::InputException();
                return result;
            }
        }

        inline Quantity eval(const ConstIterator& it) {

            if ( isValid() )
            {

                if (isNotEmpty(it,myEnd))
                {
                    mostCenteredMaximalSegment( mySC, it, myBegin, myEnd );
                    mySCEstimator.attach( mySC );
                    return mySCEstimator.eval( it );
                }
                else
                {
                    std::cerr << "[DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& it)]"
                              << " ERROR. Iterator is invalid (==myEnd)." << std::endl;
                    throw DGtal::InputException();
                    return Quantity();
                }

            }
            else
            {
                std::cerr << "[DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& it)]"
                          << " ERROR. Object is not initialized." << std::endl;
                throw DGtal::InputException();
                return Quantity();
            }
        }

    };

}

#endif //SEGBYCUT_MOSTCENTEREDMAXIMALSEGMENTESTIMATORWRAPPER_H
