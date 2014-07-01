cmake_minimum_required(VERSION 2.6)

project(PET CXX)

add_definitions(-std=c++11 -Wall -Wunused-parameter)

include_directories(src lib/cmdline)

file(GLOB COMMON
  src/geometry/*.h
  src/math/*.h
  src/util/*.h
  src/util/png_writer.cpp
  src/util/cmdline_hooks.cpp)

# OpenMP support
find_package(OpenMP)
if(OPENMP_FOUND)
  list(APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
endif()

# add libpng
find_package(PNG)
if(PNG_FOUND)
  add_definitions(-DHAVE_LIBPNG=1)
  include_directories(${PNG_INCLUDE_DIRS})
  link_libraries(${PNG_LIBRARIES})
endif()

# add CUDA specific flags
find_package(CUDA)
if(CUDA_FOUND)
  add_definitions(-DHAVE_CUDA=1)
  include_directories(${CUDA_INCLUDE_DIRS})
  list(APPEND CUDA_NVCC_FLAGS
    -Xcompiler -Wall,-Wunused-parameter,-Wno-unused-function
    -arch=compute_30 -code=sm_30)
  list(APPEND CUDA_NVCC_FLAGS_RELEASE -Xptxas -O3 --use_fast_math)
  # don't add host flags to CUDA
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  # workaround FindCUDA bug
  if(NOT CMAKE_C_COMPILER)
    include(CMakeDetermineCCompiler)
    get_filename_component(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}" REALPATH)
  endif()
endif()

# creates executable for specific PET command
function(pet_target TYPE CMD)
  # if we pass NO_CUDA then CUDA is skipped
  set(OPTIONS "${ARGN}")
  list(FIND OPTIONS NO_CUDA NO_CUDA)
  # find all possible kernels
  file(GLOB CUDA src/${TYPE}/cuda/${CMD}*.cu)
  if(CUDA_FOUND AND CUDA AND NO_CUDA EQUAL -1)
    # compile kernel sources
    cuda_compile(CUDA ${CUDA})
    # append CUDA specific headers
    file(GLOB CUDA_HEADERS
      src/${TYPE}/cuda/*.h
      src/${TYPE}/cuda/*.cuh)
    list(APPEND CUDA ${CUDA_HEADERS})
    # add runner wrapper
    if(EXISTS src/${TYPE}/cuda/${CMD}_cuda.cpp)
      list(APPEND CUDA src/${TYPE}/cuda/${CMD}_cuda.cpp)
    endif()
  endif()
  # add target specific headers
  file(GLOB HEADERS src/${TYPE}/*.h)
  # add target executable
  add_executable(${TYPE}_${CMD}
    src/${TYPE}/${CMD}_cmd.cpp
    ${HEADERS}
    ${COMMON}
    ${CUDA})
  # add CUDA objects to the target
  target_link_libraries(${TYPE}_${CMD} ${CUDA_LIBRARIES})
  # fix Intel OpenMP runtime location
  if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_custom_command(TARGET ${TYPE}_${CMD}
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL}
        -add_rpath /opt/intel/lib
        -change    libiomp5.dylib @rpath/libiomp5.dylib
        ${TYPE}_${CMD})
  endif()
endfunction()

pet_target(2d_xy matrix)
pet_target(2d_xy phantom)
pet_target(2d_xy reconstruction)
pet_target(2d_strip phantom)
pet_target(2d_strip reconstruction)
