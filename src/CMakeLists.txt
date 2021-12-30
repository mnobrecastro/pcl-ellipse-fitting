cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(pcl-ellipse-fitting)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PCL 1.11 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (ellipse_fitting ellipse_fitting.h ellipse_fitting.cpp main.cpp)
target_link_libraries (ellipse_fitting ${PCL_LIBRARIES})

#add_library(ellipse_fitting_lib STATIC ellipse_fitting.h ellipse_fitting.cpp#