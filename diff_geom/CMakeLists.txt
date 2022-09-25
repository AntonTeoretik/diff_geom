cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(diff_geom LANGUAGES CXX CUDA)
 
 
add_library(dg_lib STATIC
    algebra.cu
    color.cu
    ellipsoid.cu
    graphics.cu
    manifold.cu
    metric.cu
  )
 
target_compile_features(dg_lib PUBLIC cxx_std_17)
target_compile_options(dg_lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                       --use_fast_math
                    #    --relocatable-device-code=true
                       --extended-lambda
                       -Xptxas -O3,-v
                       >)

 
# We need to explicitly state that we need all CUDA files in the 
# particle library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties(dg_lib 
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
 
add_executable(diff_geom main.cu)
 
set_property(TARGET  diff_geom
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(diff_geom PRIVATE dg_lib)
 
