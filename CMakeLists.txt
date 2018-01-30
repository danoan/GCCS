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


set(SOURCE_FILES imageProc.cpp Artist.cpp Patch/estimatorMethods.cpp utils.cpp)

add_executable(moduleTest test.cpp  ${SOURCE_FILES})
add_executable(flowGraph flowGraph.cpp ${SOURCE_FILES} FlowGraphBuilder.cpp GluedCurve/dumb.cpp)
add_executable(minimalCase minimalCase.cpp ${SOURCE_FILES})

add_executable(dumb GluedCurve/dumb.cpp GluedCurve/GluedCurveIterator.cpp ${SOURCE_FILES})
add_executable(issue issue.cpp ${SOURCE_FILES})

target_link_libraries(moduleTest ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(flowGraph ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)
target_link_libraries(minimalCase ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lm -lboost_system -lboost_filesystem)

target_link_libraries(dumb ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lm -lboost_system -lboost_filesystem)
target_link_libraries(issue ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lm -lboost_system -lboost_filesystem)