cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME diff_geom)

project(${PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC_LIST main.cpp algebra.cpp color.cpp ellipsoid.cpp genpoints.cpp graphics.cpp manifold.cpp metric.cpp)

add_executable(${PROJECT_NAME} ${SRC_LIST})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()


add_compile_options(-O3 -fopenmp)




