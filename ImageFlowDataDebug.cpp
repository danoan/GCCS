#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "ImageFlowDataDebug.h"
#include "utils.h"

ImageFlowDataDebug::ImageFlowDataDebug(ImageFlowData& imageFlowData):imageFlowData(imageFlowData)
{
    categorizeGluedConnectors();
}

void ImageFlowDataDebug::categorizeGluedConnectors()
{
    for(auto itCP =imageFlowData.curvesPairVector.begin();itCP!=imageFlowData.curvesPairVector.end();++itCP)
    {
        GluedConnectors gluedConnectors;

        ImageFlowData::GluedCurveIteratorPair begin = itCP->gcsRangeBegin();
        ImageFlowData::GluedCurveIteratorPair end = itCP->gcsRangeEnd();

        for (ImageFlowData::GluedCurveIteratorPair itGC = begin; itGC != end; ++itGC) {
            ConnectorType ct = itGC->first.connectorType();

            switch (ct) {

                case internToExtern:
                    gluedConnectors.intConnections.push_back(itGC->first.linkSurfel());
                    break;
                case externToIntern:
                    gluedConnectors.extConnections.push_back(itGC->first.linkSurfel());
                    break;
                case makeConvex:
                    auto itC = itGC->first.connectorsBegin();
                    do{
                        gluedConnectors.makeConvexConnections.push_back(*itC);
                        if(itC==itGC->first.connectorsEnd()) break;
                        ++itC;
                    }while(true);
            }

        }

        gluedConnectors.intCurveData = itCP->intCurveData;
        gluedConnectors.extCurveData = itCP->extCurveData;

        GluedConnectorsMapKey key = getGluedConnectorMapKey(*itCP);

        gluedConnectorsMap[key] = gluedConnectors;

    }

}

ImageFlowDataDebug::GluedConnectorsMapKey ImageFlowDataDebug::getGluedConnectorMapKey(ImageFlowData::CurvePair& cp)
{
    return  GluedConnectorsMapKey(cp.intCurveData.curveType,
                                  cp.extCurveData.curveType);
}

std::string ImageFlowDataDebug::getCurveNameFromType(ImageFlowData::CurveType ct)
{
    switch(ct){
        case ImageFlowData::CurveType::OriginalCurve:
            return "Original";
        case ImageFlowData::CurveType::DilatedCurve:
            return "Dilated";
        case ImageFlowData::CurveType::ErodedCurve:
            return "Eroded";
        default:
            return "Noname";
    }
}


void ImageFlowDataDebug::drawNoConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                                       std::string outputFolder,
                                                       std::string suffix)
{
    double cmin=100;
    double cmax=-100;
    Board2D board;
    for(int i=0;i<2;++i){
        auto itCP =imageFlowData.curvesPairVector.begin();
        for(;itCP!=imageFlowData.curvesPairVector.end();++itCP)
        {
            drawCurvatureMap(itCP->intCurveData.curve.begin(),
                             itCP->intCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);
        }
        --itCP;

        drawCurvatureMap(itCP->extCurveData.curve.begin(),
                         itCP->extCurveData.curve.end(),
                         cmin,
                         cmax,
                         board,
                         weightMap);
    }


    std::string outputFilepath = outputFolder + "/noConnections" + suffix + ".eps";
    boost::filesystem::path p1(outputFilepath.c_str());
    p1.remove_filename();
    boost::filesystem::create_directories(p1);

    board.saveEPS(outputFilepath.c_str());
}

void ImageFlowDataDebug::drawInteriorConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                                             std::string outputFolder,
                                                             std::string suffix)
{
    double cmin=100;
    double cmax=-100;
    Board2D board;
    for(auto itCP =imageFlowData.curvesPairVector.begin();itCP!=imageFlowData.curvesPairVector.end();++itCP)
    {

        for(int i=0;i<2;++i){
            drawCurvatureMap(itCP->intCurveData.curve.begin(),
                             itCP->intCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);

            drawCurvatureMap(itCP->extCurveData.curve.begin(),
                             itCP->extCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);

            GluedConnectorsMapKey key = getGluedConnectorMapKey(*itCP);
            GluedConnectors gluedConnectors = gluedConnectorsMap[key];

            drawCurvatureMap(gluedConnectors.intConnections.begin(),
                             gluedConnectors.intConnections.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);
        }

        std::string filename = "/intConnections-";
        filename += getCurveNameFromType( itCP->intCurveData.curveType );
        filename += "-to-";
        filename += getCurveNameFromType( itCP->extCurveData.curveType ) + suffix;
        filename+= ".eps";

        std::string outputFilepath = outputFolder + filename;
        boost::filesystem::path p1(outputFilepath.c_str());
        p1.remove_filename();
        boost::filesystem::create_directories(p1);

        board.saveEPS(outputFilepath.c_str());
    }

}


void ImageFlowDataDebug::drawExteriorConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                                             std::string outputFolder,
                                                             std::string suffix)
{
    double cmin=100;
    double cmax=-100;
    Board2D board;
    for(auto itCP =imageFlowData.curvesPairVector.begin();itCP!=imageFlowData.curvesPairVector.end();++itCP)
    {

        for(int i=0;i<2;++i){
            drawCurvatureMap(itCP->intCurveData.curve.begin(),
                             itCP->intCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);

            drawCurvatureMap(itCP->extCurveData.curve.begin(),
                             itCP->extCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);

            GluedConnectorsMapKey key = getGluedConnectorMapKey(*itCP);
            GluedConnectors gluedConnectors = gluedConnectorsMap[key];

            drawCurvatureMap(gluedConnectors.extConnections.begin(),
                             gluedConnectors.extConnections.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);
        }

        std::string filename = "/extConnections-";
        filename += getCurveNameFromType( itCP->extCurveData.curveType );
        filename += "-to-";
        filename += getCurveNameFromType( itCP->intCurveData.curveType ) + suffix;
        filename+= ".eps";

        std::string outputFilepath = outputFolder + filename;
        boost::filesystem::path p1(outputFilepath.c_str());
        p1.remove_filename();
        boost::filesystem::create_directories(p1);

        board.saveEPS(outputFilepath.c_str());
    }

}

void ImageFlowDataDebug::drawMakeConvexConnectionsCurvatureMap(std::map<SCell,double>& weightMap,
                                                               std::string outputFolder,
                                                               std::string suffix)
{
    double cmin=100;
    double cmax=-100;
    Board2D board;
    for(auto itCP =imageFlowData.curvesPairVector.begin();itCP!=imageFlowData.curvesPairVector.end();++itCP)
    {

        for(int i=0;i<2;++i){
            drawCurvatureMap(itCP->intCurveData.curve.begin(),
                             itCP->intCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);

            drawCurvatureMap(itCP->extCurveData.curve.begin(),
                             itCP->extCurveData.curve.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);

            GluedConnectorsMapKey key = getGluedConnectorMapKey(*itCP);
            GluedConnectors gluedConnectors = gluedConnectorsMap[key];

            drawCurvatureMap(gluedConnectors.makeConvexConnections.begin(),
                             gluedConnectors.makeConvexConnections.end(),
                             cmin,
                             cmax,
                             board,
                             weightMap);
        }

        std::string filename = "/makeConvexConnections-";
        filename += suffix;
        filename += ".eps";

        std::string outputFilepath = outputFolder + filename;
        boost::filesystem::path p1(outputFilepath.c_str());
        p1.remove_filename();
        boost::filesystem::create_directories(p1);

        board.saveEPS(outputFilepath.c_str());
    }
}
