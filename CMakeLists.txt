cmake_minimum_required(VERSION 2.6)
project(segByCut)

set(CMAKE_PREFIX_PATH /home-local/dantu1/Documents/local_libs)
set(CMAKE_LIBRARY_PATH /home-local/dantu1/Documents/local_libs/lib)

find_package(DGtal REQUIRED)
find_package(OpenCV REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(GluedCurve)
INCLUDE_DIRECTORIES(Patch)

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


set(SOURCE_FILES ImageProc/imageProc.cpp utils/Artist.cpp utils/utils.cpp utils/io.cpp)
set(FLOW_FILES FlowGraph/FlowGraphBuilder.cpp FlowGraph/FlowGraphDebug.cpp FlowGraph/FlowGraphQuery.cpp
        FlowGraph/ImageFlowData.cpp FlowGraph/ImageFlowDataDebug.cpp FlowGraph/weightSettings.cpp
        FlowGraph/Patch/SimpleFlowGraphBuilder.cpp FlowGraph/Patch/SimpleFlowGraphBuilder.h)

add_executable(testEstimators testModules/testEstimators.cpp ${SOURCE_FILES})
add_executable(testFlowGraphQuery testModules/testFlowGraphQuery.cpp ${SOURCE_FILES} ${FLOW_FILES} )
add_executable(testCompImages testModules/testCompImages.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testEnergyDilation testModules/testEnergyDilation.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testAllPossibilities testModules/testAllPossibilities.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testInvariantRotation testModules/testInvariantRotation.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testSegCombinations testModules/testSegCombinations.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testFindSmallestStabbingCircles testModules/testFindSmallestStabbingCircles.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testIntExtMerge testModules/testIntExtMerge.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testSegCombQuixote testModules/testSegCombQuixote.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testBestPair BestPair/testBestPair.cpp BestPair/API.cpp BestPair/API.h BestPair/JointSet.h BestPair/JointSet.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(testFreezeOneExp testModules/testFreezeOneExpansion.cpp FreezeOneExp/CombinationsEvaluator.h ${SOURCE_FILES} ${FLOW_FILES} )
add_executable(testCapture testModules/testCapture.cpp)


add_executable(testQHalf testModules/testQHalf.cpp ${SOURCE_FILES} ${FLOW_FILES})

add_executable(presentationFigures presentationFigures.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(reportFigures reportFigures.cpp ${SOURCE_FILES} ${FLOW_FILES})

add_executable(simpleFlow Flows/simpleFlow.cpp ${SOURCE_FILES} ${FLOW_FILES} Flows/PreprocessImage.h)
add_executable(refundFlow Flows/RefundFlow/refundFlow.cpp Flows/RefundFlow/RefundFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(projectionFlow Flows/ProjectionFlow/projectionFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(adaptativeScaleFlow Flows/AdaptativeScaleFlow/adaptativeScaleFlow.cpp ${SOURCE_FILES} ${FLOW_FILES} ExhaustiveSearch/PropertyChecker/MarkedMapChecker/IntervalChecker.h ImageProc/ImageProc.h FreezeOneExp/CombinationsEvaluator.h)
add_executable(boundedFlow Flows/BoundedFlow/boundedFlow.cpp Flows/RefundFlow/RefundFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})
add_executable(modifiedFlow Flows/modifiedFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})

add_executable(binaryDataset binaryDataset.cpp Flows/RefundFlow/RefundFlow.cpp ${SOURCE_FILES} ${FLOW_FILES})




target_link_libraries(testEstimators ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testFlowGraphQuery ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testCompImages ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testEnergyDilation ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testAllPossibilities ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testInvariantRotation ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testSegCombinations ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testFindSmallestStabbingCircles ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testIntExtMerge ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(testQHalf ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(testSegCombQuixote ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(testBestPair ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(testFreezeOneExp ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(testCapture ${OpenCV_LIBS} -lseed-capture -ldrawsics)


target_link_libraries(reportFigures ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(presentationFigures ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(simpleFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(refundFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(projectionFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(adaptativeScaleFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(boundedFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(modifiedFlow ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lseed-capture -ldrawsics -lemon -lm -lboost_system -lboost_filesystem)

target_link_libraries(binaryDataset ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
