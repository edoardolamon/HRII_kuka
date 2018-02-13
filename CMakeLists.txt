cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(kuka_HRII)

find_package(XBotLogger REQUIRED)

set(CMAKE_CXX_FLAGS 	"${CMAKE_CXX_FLAGS} -std=c++11 -O3")

# we need Eigen for the algebra in the code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# for KUKA Robot
add_definitions(-DHAVE_SOCKLEN_T)

include_directories(include FRI_Remote)

# add the source files to be compiled
add_executable(kuka_firstapp src/fri_firstapp.cpp src/friudp.cpp src/friremote.cpp)
add_executable(kuka_example src/fri_simpleExample.cpp src/friudp.cpp src/friremote.cpp)
add_executable(kuka_control src/kuka_cartesianImpedenceControl.cpp src/friudp.cpp src/friremote.cpp)
add_executable(utils_test src/utils_test.cpp)
add_executable(kuka_control_test src/KukaLWRControl_test.cpp src/KukaLWRInterface.cpp src/friudp.cpp src/friremote.cpp)

target_link_libraries(kuka_control ${XBotLogger_LIBRARIES})
target_link_libraries(utils_test ${XBotLogger_LIBRARIES})
target_link_libraries(kuka_control_test ${XBotLogger_LIBRARIES})


