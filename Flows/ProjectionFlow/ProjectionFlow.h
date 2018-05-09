#ifndef SEGBYCUT_PROJECTIONFLOW_H
#define SEGBYCUT_PROJECTIONFLOW_H

#include <opencv/highgui.h>

#include "../../FlowGraph/ImageFlowData.h"
#include "../../FlowGraph/FlowGraphBuilder.h"
#include "../../FlowGraph/weightSettings.h"
#include "../../FlowGraph/FlowGraphQuery.h"
#include "../../FlowGraph/ImageUpdater.h"

#include "CurveProjector.h"
#include "../../utils/io.h"

class ProjectionFlow
{
public:
    friend class FlowGraphDebug;

    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
    typedef FlowGraphBuilder::LinelWeightMap LinelWeightMap;

public:
    ProjectionFlow(ImageFlowData& imageFlowData):imageFlowData(imageFlowData)
    {

    }

    void operator()(Image2D& imageSolution)
    {
        LinelWeightMap weightMap;
        setArcsWeight(imageFlowData,weightMap);

        ImageUpdater imgUpdater(imageFlowData.getKSpace());

        double currentEnergyValue=0;
        double previousEnergyValue;
        do
        {
            previousEnergyValue = currentEnergyValue;

            FlowGraph fg;
            FlowGraphBuilder fgb(fg,imageFlowData,weightMap);

            currentEnergyValue = FlowGraphQuery::cutValue(fg);
            imgUpdater(imageSolution, fg, imageFlowData.getOriginalImage());

            std::cout << "Current: " << currentEnergyValue << std::endl;
            std::cout << "Previous: " << previousEnergyValue << std::endl;

            ImageFlowData partialIMF(imageSolution);
            partialIMF.init(imageFlowData.getFlowMode(), imageFlowData.getGluedCurveLength());

            LinelWeightMap partialWeightMap;
            setArcsWeight(partialIMF,partialWeightMap);

            //Update Glued Arcs weights
            std::cout << "Weight Update:" << std::endl;
            for(auto itc =partialIMF.curvePairBegin();itc!=partialIMF.curvePairEnd();++itc)
            {
                for(auto it = itc->gcsRangeBegin();it!=itc->gcsRangeEnd();++it)
                {
                    SCell connector = it->first.linkSurfel();
                    if(weightMap.find(connector)!=weightMap.end())
                    {
                        if(weightMap[connector]!=partialWeightMap[connector])
                            std::cout << weightMap[connector] << " -> " << partialWeightMap[connector] << std::endl;
                        weightMap[connector] = partialWeightMap[connector];
                    }

                }

            }

            switch(imageFlowData.getFlowMode())
            {
                case ImageFlowData::DilationOnly:
                    dilationProjectUpdate(partialIMF,imageFlowData,partialWeightMap,weightMap);
                    break;
                case ImageFlowData::ErosionOnly:
                    erosionProjectUpdate(partialIMF,imageFlowData,partialWeightMap,weightMap);
                    break;
                default:
                    dilationProjectUpdate(partialIMF,imageFlowData,partialWeightMap,weightMap);
                    erosionProjectUpdate(partialIMF,imageFlowData,partialWeightMap,weightMap);
            }


        }while(currentEnergyValue!=previousEnergyValue);

        _cutValue = currentEnergyValue;
    }
    
    double cutValue(){return _cutValue;}

private:

    void projectAndUpdateWeights(KSpace& KImage,
                                 Curve& toProject,
                                 Curve& projectOn,
                                 LinelWeightMap& toProjectWeightMap,
                                 LinelWeightMap& projectOnWeightMap)
    {
        CurveProjector projector(KImage,
                                 toProject,
                                 projectOn);

        CurveProjector::Projection projectionMap;
        projector(projectionMap);

        std::cout << "Projection Update:" << std::endl;
        auto it = toProject.begin();
        do
        {
            SCell projectedSCell = projectionMap[*it];
            if(projector.validProjection(projectedSCell))
            {
                if(projectOnWeightMap[projectedSCell]!=toProjectWeightMap[*it])
                    std::cout << projectOnWeightMap[projectedSCell] << " -> " << toProjectWeightMap[*it] << std::endl;
                projectOnWeightMap[projectedSCell] = toProjectWeightMap[*it];
            }

            ++it;
        }while(it!=toProject.end());
    }

    void dilationProjectUpdate(ImageFlowData& toProjectIMF,
                               ImageFlowData& projectOnIMF,
                               LinelWeightMap& toProjectWeightMap,
                               LinelWeightMap& projectOnWeightMap)
    {
        Curve& toProject = toProjectIMF.getMostOuterCurve();
        Curve& projectOn = projectOnIMF.getMostOuterCurve();

        projectAndUpdateWeights(projectOnIMF.getKSpace(),
                                toProject,
                                projectOn,
                                toProjectWeightMap,
                                projectOnWeightMap);
    }

    void erosionProjectUpdate(ImageFlowData& toProjectIMF,
                              ImageFlowData& projectOnIMF,
                              LinelWeightMap& toProjectWeightMap,
                              LinelWeightMap& projectOnWeightMap)
    {
        Curve& toProject = toProjectIMF.getMostInnerCurve();
        Curve& projectOn = projectOnIMF.getMostInnerCurve();

        projectAndUpdateWeights(projectOnIMF.getKSpace(),
                                toProject,
                                projectOn,
                                toProjectWeightMap,
                                projectOnWeightMap);
    }

    void displayImage(std::string windowName,
                      std::string imagePath)
    {
        cv::imshow( windowName,
                    cv::imread(imagePath,CV_8U)
        );
        cv::waitKey(0);
    }



private:
    ImageFlowData& imageFlowData;
    double _cutValue;
};

#endif //SEGBYCUT_PROJECTIONFLOW_H
