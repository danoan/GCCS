#ifndef SEGBYCUT_IMAGEFLOWDATADEBUG_H
#define SEGBYCUT_IMAGEFLOWDATADEBUG_H

#include "ImageFlowData.h"
#include "FlowGraphBuilder.h"

class ImageFlowDataDebug{

public:
    ImageFlowDataDebug(ImageFlowData& imageFlowData);

    typedef DGtal::Z2i::SCell SCell;

    typedef std::pair<ImageFlowData::CurveType, ImageFlowData::CurveType> GluedConnectorsMapKey;

    typedef struct _GluedConnectors{
        std::vector<SCell> intConnections;
        std::vector<SCell> extConnections;
        std::vector<SCell> makeConvexConnections;

        ImageFlowData::CurveData intCurveData;
        ImageFlowData::CurveData extCurveData;
    } GluedConnectors;

    void drawNoConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                       std::string outputFolder,
                                       std::string suffix);

    void drawInteriorConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                             std::string outputFolder,
                                             std::string suffix);

    void drawExteriorConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                             std::string outputFolder,
                                             std::string suffix);

    void drawMakeConvexConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                               std::string outputFolder,
                                               std::string suffix);


private:
    void categorizeGluedConnectors();
    GluedConnectorsMapKey getGluedConnectorMapKey(ImageFlowData::CurvePair& cp);
    std::string getCurveNameFromType(ImageFlowData::CurveType ct);

private:
    ImageFlowData& imageFlowData;

    std::map< GluedConnectorsMapKey, GluedConnectors> gluedConnectorsMap;
};

#endif //SEGBYCUT_IMAGEFLOWDATADEBUG_H
