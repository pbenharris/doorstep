# To build with debug
# cmake -DCMAKE_BUILD_TYPE=Debug .. -B .
cmake_minimum_required(VERSION 3.22)
project(Simple)

find_package(Matplot++ REQUIRED)

add_executable(doorstep doorstep.cpp)
target_include_directories(doorstep PUBLIC
  ../dependencies/glob/glob-master/single_include/glob
  ../dependencies/json-develop/include
  ../dependencies/spdlog-1.x/include
  )

target_link_libraries(doorstep Matplot++::matplot)

