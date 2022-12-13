cmake_minimum_required(VERSION 3.22)
project(project CUDA)

set(CMAKE_CUDA_STANDARD 20)

set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -fmad=false --expt-relaxed-constexpr --compiler-options='-ffp-contract=fast'")
set(CMAKE_CUDA_FLAGS_DEBUG "-fmad=false --expt-relaxed-constexpr --compiler-options='-ffp-contract=fast'")

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIR})

add_executable(project main.cu helper.cuh myMath.cuh Vec3.cuh Mesh.cuh Scatter.cuh Sphere.cuh Vec2.cuh Bessel.cuh)

configure_file(parameter.txt parameter.txt COPYONLY)

set_target_properties(project PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(project fftw3 ${OpenCV_LIBS})