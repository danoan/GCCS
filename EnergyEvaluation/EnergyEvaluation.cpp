#include "EnergyEvaluation.h"

namespace EnergyEvaluation
{
    double integratedSquaredCurvature(DigitalSet &ds)
    {
        Domain domain = ds.domain();
        KSpace KImage;
        KImage.init(domain.lowerBound(),domain.upperBound(),true);

        DigitalSet boundary(domain);
        MyBoundary(boundary,ds);

        SurfelAdjacency SAdj(true);
        SCell aBel = Surfaces::findABel(KImage, boundary, 10000);

        Curve curve;
        std::vector<SCell> scells;
        Surfaces::track2DBoundary(scells,KImage,SAdj,boundary,aBel);

        curve.initFromSCellsVector(scells);

        WeightMap weightMap;
        setGridCurveWeight(curve,
                           KImage,
                           weightMap);

        double energy = 0;
        for(auto it=weightMap.begin();it!=weightMap.end();++it)
        {
            energy += it->second;
        }

        DGtal::Board2D board;
        board << curve;
        board.saveEPS("curve.eps");

        return energy;
    }
}