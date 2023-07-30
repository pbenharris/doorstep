# To build with debug
# cmake -DCMAKE_BUILD_TYPE=Debug .. -B .
cmake_minimum_required(VERSION 3.22)
project(doorstep)

find_package(Matplot++ REQUIRED)

set(CMAKE_BUILD_TYPE Debug)

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.74.0 REQUIRED)

add_executable(doorstep doorstep.cpp)
target_include_directories(doorstep PUBLIC
  ../dependencies/glob/glob-master/single_include/glob
  ../dependencies/json-develop/include
  ../dependencies/spdlog-1.x/include
  ${Boost_INCLUDE_DIRS}
  )

target_link_libraries(doorstep Matplot++::matplot ${Boost_LIBRARIES})

