cmake_minimum_required(VERSION 2.6)

project(PET CXX)

add_definitions(-std=c++11)

include_directories(src lib/cmdline)

file(GLOB COMMON
  std/geometry/*.h
  std/math/*.h
  std/util/*.h
  src/util/png_writer.cpp
  src/util/cmdline_hooks.cpp)

find_package(OpenMP)

# add libpng
find_package(PNG)
if(PNG)
  add_definitions(-DHAVE_LIBPNG=1)
endif()

# add CUDA specific flags
find_package(CUDA)
if(CUDA_FOUND)
  add_definitions(-DHAVE_CUDA=1)
  if(APPLE)
    set(CUDA_HOST_COMPILER /usr/bin/clang)
  endif()
  list(APPEND CUDA_NVCC_FLAGS         -arch=compute_30 -code=sm_30)
  list(APPEND CUDA_NVCC_FLAGS_RELEASE -Xptxas -O3 --use_fast_math)
endif()

# creates executable for specific PET command
macro(pet_target TYPE CMD)
  # handle CUDA wrappers
  if(CUDA_FOUND AND EXISTS src/${TYPE}/cuda/${CMD}_kernel.cu)
    # append CUDA specific headers
    file(GLOB ${TYPE}_CUDA_HEADERS src/${TYPE}/cuda/*.cuh)
    # compile kernel sources
    cuda_compile(${TYPE}_${CMD}_CUDA src/${TYPE}/cuda/${CMD}_kernel.cu)
    # add runner wrapper
    if(EXISTS src/${TYPE}/cuda/${CMD}_cuda.cpp)
      list(APPEND ${TYPE}_${CMD}_CUDA
        ${${TYPE}_CUDA_HEADERS}
        src/${TYPE}/cuda/${CMD}_cuda.cpp)
    endif()
  endif()

  # add target specific headers
  file(GLOB ${TYPE}_HEADERS src/${TYPE}/*.h)

  # add target executable
  add_executable(${TYPE}_${CMD}
    src/${TYPE}/${CMD}_cmd.cpp
    ${${TYPE}_HEADERS}
    ${COMMON}
    ${${TYPE}_${CMD}_CUDA})

  target_link_libraries(${TYPE}_${CMD} ${PNG_LIBRARIES} ${CUDA_LIBRARIES})
endmacro()

pet_target(2d_xy matrix)
# pet_target(2d_xy phantom)
pet_target(2d_xy reconstruction)
pet_target(2d_strip phantom_strip)
# pet_target(2d_strip reconstruction_strip)
