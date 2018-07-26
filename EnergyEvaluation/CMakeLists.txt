cmake_minimum_required(VERSION 2.6)
project(energy-evaluation)

set(LOCAL_INCLUDE /usr/local/lib CACHE PATH "Alternative header directory")
set(LOCAL_LIB /usr/local/include CACHE PATH "Alternative library directory")

INCLUDE_DIRECTORIES(${LOCAL_INCLUDE})
LINK_DIRECTORIES(${LOCAL_LIB})

find_package(DGtal REQUIRED)
find_package(OpenCV REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

INCLUDE_DIRECTORIES(../)

add_executable(testEnergyEvaluation testEnergyEvaluation.cpp EnergyEvaluation.cpp utils.cpp)
target_link_libraries(testEnergyEvaluation ${DGTAL_LIBRARIES} ${OpenCV_LIBS} -lemon -lm -lboost_system -lboost_filesystem)

include(GNUInstallDirs)
add_library(energy-eval STATIC EnergyEvaluation.cpp utils.cpp )

set_target_properties(energy-eval PROPERTIES COMPILE_FLAGS "-std=c++11")

configure_file(energy-eval.pc.in energy-eval.pc @ONLY)

install(TARGETS energy-eval ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES EnergyEvaluation.h utils.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EnergyEvaluation)
install(DIRECTORY ../Patch DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EnergyEvaluation)
