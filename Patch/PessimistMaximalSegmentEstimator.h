#ifndef SEGBYCUT_PESSIMISTMAXIMALSEGMENTESTIMATOR_H
#define SEGBYCUT_PESSIMISTMAXIMALSEGMENTESTIMATOR_H

#include <boost/concept/assert.hpp>
#include <boost/static_assert.hpp>
#include <DGtal/base/CUnaryFunctor.h>

#include "DGtal/base/Common.h"
#include "DGtal/base/Exceptions.h"
#include "DGtal/base/Circulator.h"
#include "DGtal/base/IteratorFunctions.h"

#include "DGtal/geometry/curves/estimation/CSegmentComputerEstimator.h"
#include "DGtal/geometry/curves/CForwardSegmentComputer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"

//SegmentComputerUtils -- mostCenteredMaximalSegment

namespace Patch {
    template<typename SegmentComputer, typename SCEstimator>
    class PessimistMaximalSegmentEstimator {

        BOOST_CONCEPT_ASSERT((DGtal::concepts::CForwardSegmentComputer<SegmentComputer>));
        BOOST_CONCEPT_ASSERT((DGtal::concepts::CSegmentComputerEstimator<SCEstimator>));

        BOOST_STATIC_ASSERT((boost::is_same<SegmentComputer,
                typename SCEstimator::SegmentComputer>::value));

        // ----------------------- Types ------------------------------
    public:

        typedef typename SegmentComputer::ConstIterator ConstIterator;
        typedef typename SCEstimator::Quantity Quantity;

        typedef DGtal::SaturatedSegmentation<SegmentComputer> Segmentation;
        typedef typename Segmentation::SegmentComputerIterator SegmentIterator;

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor. Not valid.
         */
        PessimistMaximalSegmentEstimator();

        /**
         * Constructor.
         * @param aSegmentComputer a segment computer
         * @param aSCEstimator an estimator
         */
        PessimistMaximalSegmentEstimator(const SegmentComputer &aSegmentComputer,
                                      const SCEstimator &aSCEstimator);

        /**
         * Destructor.
         */
        ~PessimistMaximalSegmentEstimator() {};

        // ----------------------- Interface --------------------------------------
    public:

        /**
         * Initialisation.
         * @param h grid size (must be >0).
         * @param itb begin iterator
         * @param ite end iterator
         */
        void init(const double h, const ConstIterator &itb, const ConstIterator &ite);

        /**
         * Unique estimation
         * @param it any valid iterator
         * @return the estimated quantity at *it
         *
         * NB: the whole range [@e myBegin , @e myEnd)|
         * is scanned in the worst case
         */
        Quantity eval(const ConstIterator &it);

        /**
         * Estimation for a subrange [@e itb , @e ite )
         *
         * @param itb subrange begin iterator
         * @param ite subrange end iterator
         * @param result output iterator on the estimated quantity
         *
         * @return the estimated quantity
         * from itb till ite (excluded)
         *
         * NB: the whole range [@e myBegin , @e myEnd)|
         * is scanned in the worst case
         */
        template<typename OutputIterator>
        OutputIterator eval(const ConstIterator &itb, const ConstIterator &ite,
                            OutputIterator result);


        /**
         * Checks the validity/consistency of the object.
         * @return 'true' if the object is valid, 'false' otherwise.
         */
        bool isValid() const;

        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:

        /** grid step */
        double myH;

        /** begin and end iterators */
        ConstIterator myBegin, myEnd;

        /** segmentComputer used to segment */
        SegmentComputer mySC;

        /** object estimating the quantity from segmentComputer */
        SCEstimator mySCEstimator;

        // ------------------------- Internal services ------------------------------


    private:

        /**
         * Copy constructor.
         * @param other the object to clone.
         * Forbidden by default.
         */
        PessimistMaximalSegmentEstimator(const PessimistMaximalSegmentEstimator &other);

        /**
         * Assignment.
         * @param other the object to copy.
         * @return a reference on 'this'.
         * Forbidden by default.
         */
        PessimistMaximalSegmentEstimator &operator=(const PessimistMaximalSegmentEstimator &other);


    };
}

#include "PessimistMaximalSegmentEstimator.ih"


#endif //SEGBYCUT_PESSIMISTMAXIMALSEGMENTESTIMATOR_H
