cmake_minimum_required(VERSION 3.13)
project(final_lagrange)

find_package(CUDA QUIET REQUIRED)
message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")

## Pass options to NVCC
set(
        CUDA_NVCC_FLAGS
        ${CUDA_NVCC_FLAGS};
        -G
)

# Specify include directories
include_directories(
        /home/matt/CUDA_SAMPLES/NVIDIA_CUDA-10.0_Samples/common/inc
)

#cuda_add_library(
#        cuda_lagrange
#        lagrange_cuda.cu
#        lagrange_cuda.cuh
#)

cuda_add_executable(
        lagrange
        lagrange_cuda.cu
        main.cu
        common.cpp
        World.cpp
        globals.cpp
        cpu.cu
)
#
#
#target_link_libraries(lagrange cuda_lagrange)