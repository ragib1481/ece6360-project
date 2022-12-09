cmake_minimum_required(VERSION 3.22)
project(project CUDA)

set(CMAKE_CUDA_STANDARD 20)

set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -fmad=false --expt-relaxed-constexpr --compiler-options='-ffp-contract=fast'")

find_package(OpenCV REQUIRED)
find_package(Boost REQUIRED)

include_directories(${OPENCV_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(project main.cu helper.h myMath.h Vec3.cuh Mesh.cuh Scatter.cuh)

configure_file(parameter.txt parameter.txt COPYONLY)

set_target_properties(project PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(project fftw3 ${OpenCV_LIBS} ${Boost_LIBRARIES})