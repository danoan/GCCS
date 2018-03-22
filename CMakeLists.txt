cmake_minimum_required(VERSION 2.6)
project(segByCut)

set(CMAKE_PREFIX_PATH /home/dantu/Documents/local_libs)
set(CMAKE_LIBRARY_PATH /home/dantu/Documents/local_libs/lib)

find_package(DGtal REQUIRED)
find_package(OpenCV REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(GluedCurve)
INCLUDE_DIRECTORIES(Patch)

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


set(SOURCE_FILES imageProc.cpp Artist.cpp utils.cpp io.cpp)
set(FLOW_FILES FlowGraph/FlowGraphBuilder.cpp FlowGraph/FlowGraphDebug.cpp FlowGraph/FlowGraphQuery.cpp
        FlowGraph/ImageFlowData.cpp FlowGraph/ImageFlowDataDebug.cpp FlowGraph/weightSettings.cpp)

add_executable(testEstimators testModules/testEstimators.cpp ${SOURCE_FILES})
add_executable(testFlowGraphQuery testModules/testFlowGraphQuery.cpp ${SOURCE_FILES} ${FLOW_FILES} )
add_executable(testCompImages testModules/testCompImages.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testEnergyDilation testModules/testEnergyDilation.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testAllPossibilities testModules/testAllPossibilities.cpp ${SOURCE_FILES} ${FLOW_FILES})

add_executable(issue issue.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(presentationFigures presentationFigures.cpp ${SOURCE_FILES} ${FLOW_FILES})

add_executable(simpleFlow simpleFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(refundFlow refundFlow.cpp RefundFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})

add_executable(binaryDataset binaryDataset.cpp RefundFlow.cpp ${SOURCE_FILES} ${FLOW_FILES} ExhaustiveSearch/CombinationsEvaluator.h)





target_link_libraries(testEstimators ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testFlowGraphQuery ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testCompImages ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testEnergyDilation ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testAllPossibilities ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(issue ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(presentationFigures ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(simpleFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(refundFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(binaryDataset ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
