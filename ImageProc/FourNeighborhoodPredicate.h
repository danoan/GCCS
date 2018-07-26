#ifndef GLOBALPAIRWISESEGMENTATION_FOURNEIGHBORHOODPREDICATE_H
#define GLOBALPAIRWISESEGMENTATION_FOURNEIGHBORHOODPREDICATE_H


#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/images/CConstImage.h"
#include "DGtal/base/ConstAlias.h"

#include "boost/concept/assert.hpp"
#include "boost/assert.hpp"

namespace ImageProc
{
    template<typename DigitalSet>
    class FourNeighborhoodPredicate {
    public:
        typedef DGtal::Z2i::Domain Domain;
        typedef DGtal::Z2i::Point Point;

        FourNeighborhoodPredicate(DigitalSet &DS) :
                myDigitalSet(DS) {
            lowerBound = DS.domain().lowerBound();
            upperBound = DS.domain().upperBound();
        };

        bool operator()(const Point &aPoint) const {
            Point np;
            int s = 0;
            for (int i = 0; i < 4; ++i) {
                np = aPoint + filter[i];
                if (np[0] < lowerBound[0] || np[1] < lowerBound[1]) continue;
                if (np[0] > upperBound[0] || np[1] > upperBound[1]) continue;


                s += myDigitalSet(np) ? 1 : 0;
            }

            return !(s > 0 && s < 4);
        }

        bool operator()(const Domain::ConstIterator &it) const {
            return (*this)(*it);
        }


    private:
        const DigitalSet &myDigitalSet;
        Point filter[4] = {Point(0, 1), Point(1, 0),
                           Point(-1, 0), Point(0, -1)};

        Point lowerBound, upperBound;


    protected:
        FourNeighborhoodPredicate();

    };
}

#endif //GLOBALPAIRWISESEGMENTATION_FOURNEIGHBORHOODPREDICATE_H
