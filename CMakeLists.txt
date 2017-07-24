cmake_minimum_required(VERSION 3.7)
project(Nodesv5)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp MultiPathProcessing.cpp random.cpp functions.cpp)
add_executable(Nodesv5 ${SOURCE_FILES})
find_package(GSL REQUIRED)
target_link_libraries(Nodesv5 GSL::gsl GSL::gslcblas)

include_directories(${GSL_INCLUDE_DIR})