#ifndef SEGBYCUT_LAMBDAMAXIMALSEGMENTESTIMATOR_H
#define SEGBYCUT_LAMBDAMAXIMALSEGMENTESTIMATOR_H


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
#include "DGtal/geometry/curves/estimation/FunctorsLambdaMST.h"

//SegmentComputerUtils -- mostCenteredMaximalSegment

namespace Patch {
    /**
     *  Polynomial lambda functor \cite LachaudIVC2007.
     * \f$ 16 ( x^4 - 2 x^3 + x^2 ) \f$
     *
     */
    struct Lambda16Function : std::unary_function < double, double >
    {
        double operator() (double x) const
        {
            double e2 = x * x;
            double e3 = e2 * x;
            return 16.0 * ( e3 *x - 2.0 * e3 + e2 );
        }
    };


    template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction=DGtal::functors::Lambda64Function>
    class LambdaMaximalSegmentEstimator {

        BOOST_CONCEPT_ASSERT((DGtal::concepts::CForwardSegmentComputer<SegmentComputer>));
        BOOST_CONCEPT_ASSERT((DGtal::concepts::CSegmentComputerEstimator<SCEstimator>));
        BOOST_CONCEPT_ASSERT((DGtal::concepts::CUnaryFunctor<LambdaFunction, double, double>));

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
        LambdaMaximalSegmentEstimator();

        /**
         * Constructor.
         * @param aSegmentComputer a segment computer
         * @param aSCEstimator an estimator
         */
        LambdaMaximalSegmentEstimator(const SegmentComputer &aSegmentComputer,
                                      const SCEstimator &aSCEstimator);

        /**
         * Destructor.
         */
        ~LambdaMaximalSegmentEstimator() {};

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

        LambdaFunction lambdaFunction;

        // ------------------------- Internal services ------------------------------


    private:

        /**
         * Copy constructor.
         * @param other the object to clone.
         * Forbidden by default.
         */
        LambdaMaximalSegmentEstimator(const LambdaMaximalSegmentEstimator &other);

        /**
         * Assignment.
         * @param other the object to copy.
         * @return a reference on 'this'.
         * Forbidden by default.
         */
        LambdaMaximalSegmentEstimator &operator=(const LambdaMaximalSegmentEstimator &other);


    };
}

#include "LambdaMaximalSegmentEstimator.ih"


#endif //SEGBYCUT_LAMBDAMAXIMALSEGMENTESTIMATOR_H
