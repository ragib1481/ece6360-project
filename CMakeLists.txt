cmake_minimum_required(VERSION 3.23)
project(project CUDA)

set(CMAKE_CUDA_STANDARD 20)

set(CMAKE_CUDA_FLAGS_RELEASE "-O3")

find_package(OpenCV REQUIRED)

include_directories(${OPENCV_INCLUDE_DIRS})

add_executable(project main.cu Propagator.cuh helper.h Beam.cuh)

configure_file(parameter.txt parameter.txt COPYONLY)

set_target_properties(project PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(project fftw3 ${OpenCV_LIBS})