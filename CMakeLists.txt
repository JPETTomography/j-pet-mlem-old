# PET project main description file
#
# Author: Adam Strzelecki <adam.strzelecki@uj.edu.pl>
#
# Description:
#   Builds all PET project commands using CMake build generator.

cmake_minimum_required(VERSION 2.8.10)

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW) # do not dereference quoted variables
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

project(PET CXX) # using only C++

# configuration time options
set(PET_MAX_DETECTORS 512
    CACHE INTEGER
    "maximum number of detectors allowed in the scanner (impacts performance)"  )
set(PET_GRANULARITY_TYPE "WARP"
    CACHE STRING
    "warp, thread (or simple just for testing) granularity"  )
set(PET_USE_STATISTICS OFF
    CACHE BOOL
    "collect statistics and print them after program is finished"  )
set(PET_CACHE_ELLIPSE_PIXELS OFF
    CACHE BOOL
    "shared memory pixel buffer in error ellipse"  )
set(PET_USE_SENSITIVITY ON
    CACHE BOOL
    "use sensitivity when calculating reconstruction"  )
set(PET_USE_KERNEL ON
    CACHE BOOL
    "use kernel (when off it uses constant 1 for kernel value)"  )
set(PET_THROW_ON_ZERO_DENOMINATOR OFF
    CACHE BOOL
    "throw on zero denominator (event outside FOV)"  )
set(PET_MAX_PIXELS_PER_THREAD 12
    CACHE INTEGER
    "chosen arbitrarily, no more than 11 pixels per thread in current config"  )
set(PET_WARP_SIZE 32
    CACHE INTEGER
    "normally it should not be changed to anything else than 32"  )
set(PET_MAX_THREADS_PER_BLOCK 512
    CACHE INTEGER
    "not used in current implementation"  )
set(PET_USE_FAST_TEXT_PARSER ON
    CACHE BOOL
    "use fast text parser for loading events (around ~5x faster than C++)"  )

# extra verbosity
option(CUDA_PTXAS_VERBOSE "show verbose information about CUDA registers" OFF)

# use Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type Debug or Release" FORCE)
endif()

# check supported compiler versions
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU"   AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8.1) OR
   (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3) OR
   (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0) OR
   (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC"  AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.0))
   # we got supported compiler
else()
  message(FATAL_ERROR
    "Selected compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} is not supported.\n"
    "GNU/GCC >= 4.8.1, Clang >= 3.3, Intel/ICC >= 15.0 or MSVC >= 2015 is required for this project.")
endif()

if(MLEM_SILENT_BUILD)
  set(WARNING_FLAGS "")
  set(CUDA_WARNING_FLAGS "")
  message("#################################################")
  message("##         All warnings are supressed          ##")
  message("#################################################")
else()
  set(WARNING_FLAGS "-Wall -Wunused-parameter -Wno-comment")
  set(CUDA_WARNING_FLAGS "-Wall,-Wno-unused-function")
  message("#################################################")
  message("##Set MLEM_SILENT_BUILD to supress all warnings##")
  message("#################################################")
endif()
# enable C++11 and warnings
if(NOT MSVC)
  add_definitions(-std=c++11 ${WARNING_FLAGS})
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
  if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
  elseif(UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
  endif()
else()
  add_definitions(/D_USE_MATH_DEFINES /wd4244 /wd4800 /wd4267 /wd4996)
endif()

include_directories(src lib/cmdline lib/json)

# common headers and utility sources used by all commands
file(GLOB COMMON
  src/geometry/*.h
  src/2d/geometry/*.h
  src/3d/geometry/*.h
  src/math/*.h
  src/util/*.h
  src/util/png_writer.cpp
  src/util/cmdline_hooks.cpp
  src/util/backtrace.cpp
  src/util/string.cpp
  src/common/*.h
  src/common/options.cpp)

# OpenMP support
find_package(OpenMP)
if(OPENMP_FOUND)
  message(STATUS "OpenMP found")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# CUDA support
if(CMAKE_VERSION VERSION_GREATER 3.3)
  set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "use the static CUDA runtime")
endif()
find_package(CUDA)
if(CUDA_FOUND AND CUDA_VERSION VERSION_LESS 7.5)
  message(STATUS "Error: CUDA 7.5 or greater required for this project")
  unset(CUDA_FOUND)
endif()
if(CUDA_FOUND)
  add_definitions(-DHAVE_CUDA=1)
  include_directories(${CUDA_INCLUDE_DIRS})
  list(APPEND CUDA_NVCC_FLAGS
    -gencode arch=compute_62,code=compute_62 # GTX 1080
    -gencode arch=compute_52,code=compute_52 # GTX 980Ti
    -gencode arch=compute_30,code=compute_30 # GTX 770
    -std=c++11)
  list(APPEND CUDA_NVCC_FLAGS_RELEASE
    -Xptxas -O3,-lineinfo
    --use_fast_math)
  if(CUDA_PTXAS_VERBOSE)
    list(APPEND CUDA_NVCC_FLAGS_RELEASE -Xptxas -v)
  endif()
  list(APPEND CUDA_NVCC_FLAGS_DEBUG -g -G)
  if(NOT MSVC)
    if(NOT MLEM_SILENT_BUILD)
      list(APPEND CUDA_NVCC_FLAGS -Xcompiler ${CUDA_WARNING_FLAGS})
    else()
      list(APPEND CUDA_NVCC_FLAGS -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-deprecated-declarations)
     
    endif()
  else()
    list(APPEND CUDA_NVCC_FLAGS -Xcompiler /wd4244,/wd4800,/wd4267,/wd4996)
    list(APPEND CUDA_NVCC_FLAGS_RELEASE -Xcompiler /MD)
    list(APPEND CUDA_NVCC_FLAGS_DEBUG   -Xcompiler /MDd)
  endif()
  # don't add host flags to CUDA
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  # workaround FindCUDA bug when no C compiler is specified
  if(NOT CMAKE_C_COMPILER AND NOT CUDA_HOST_COMPILER)
    if(NOT MSVC)
      include(CMakeDetermineCCompiler)
      get_filename_component(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}" REALPATH)
    else()
      set(CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
    endif()
    set(CUDA_HOST_COMPILER "${CUDA_HOST_COMPILER}" CACHE PATH "CUDA host compiler path" FORCE)
  endif()
endif()

# libpng for output image generation
find_package(PNG)
if(NOT PNG_FOUND AND MSVC)
  include(ExtPNG) # builds libpng automatically
endif()
if(PNG_FOUND)
  add_definitions(-DHAVE_LIBPNG=1 ${PNG_DEFINITIONS})
  include_directories(${PNG_INCLUDE_DIRS})
  link_libraries(${PNG_LIBRARIES})
endif()

# optionally use Boost
find_package(Boost 1.58.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  add_definitions(-DHAVE_BOOST=1 ${Boost_DEFINITIONS})
endif()

# we need -ldl for backtrace on Linux
if(UNIX AND NOT APPLE)
  link_libraries(-ldl)
endif()

# add options above as compile flags removing PET_ prefix
get_cmake_property(vars VARIABLES)
foreach(var ${vars})
  if(var MATCHES "^PET_" AND NOT var MATCHES "_DIR$")
    string(REGEX REPLACE "^PET_" "" name ${var})
    if(name MATCHES "_TYPE$")
      string(REGEX REPLACE "_TYPE$" "" type_name ${name})
      add_definitions(-DUSE_${${var}}_${type_name}=1)
      add_definitions(-D${name}=${${var}})
      message(STATUS "Using ${${var}} ${type_name}")
    elseif(${var})
      if(${var} STREQUAL ON)
        add_definitions(-D${name}=1)
        message(STATUS "Enabled ${name}")
      else()
        message(STATUS "Using ${name}=${${var}}")
        add_definitions(-D${name}=${${var}})
      endif()
    else()
      add_definitions(-D${name}=0)
      message(STATUS "Disabled ${name}")
    endif()
  endif()
endforeach()



file(GLOB CUDA_COMMON_SRC src/common/cuda/*.cu)
if(CUDA_FOUND AND CUDA_COMMON_SRC )
  cuda_compile(_CUDA_COMMON ${CUDA_COMMON_SRC})
  set(CUDA_COMMON ${_CUDA_COMMON})
  message(${CUDA_COMMON})
endif()




# custom function that creates executable for specific PET command
# see below for command list generated by this function
function(pet_target DIM TYPE CMD)
  # if we pass NO_CUDA then CUDA is skipped
  list(FIND ARGN NO_CUDA NO_CUDA)
  # find all possible kernels
  file(GLOB CUDA src/${DIM}/${TYPE}/cuda/${CMD}*.cu)
  if(CUDA_FOUND AND CUDA AND NO_CUDA EQUAL -1)
    #file(GLOB CUDA_COMMON src/common/cuda/*.cu)
    # compile kernel sources
    cuda_compile(CUDA ${CUDA} )
    #message(${CUDA})
    # append CUDA specific headers & wrappers
    file(GLOB CUDA_HEADERS
      src/${DIM}/${TYPE}/cuda/*.h
      src/util/cuda/*.h
      src/${DIM}/${TYPE}/cuda/*.cuh
      src/${DIM}/${TYPE}/cuda/**/*.cuh
      src/${DIM}/${TYPE}/cuda/${CMD}*.cpp
      src/common/cuda/*.h
      src/common/cuda/*.cuh)
    list(APPEND CUDA ${CUDA_HEADERS})
  endif()
  # add other files
  file(GLOB others src/${DIM}/${TYPE}/*.cpp)
  foreach(other ${others})
    if(NOT other MATCHES "_(test|cmd).cpp$")
      list(APPEND OTHER ${other})
    endif()
  endforeach()
  # add extra files
  foreach(extra_glob ${ARGN})
    file(GLOB extras src/${extra_glob})
    list(APPEND OTHER ${extras})
  endforeach()
  # add target specific headers
  file(GLOB HEADERS src/${DIM}/${TYPE}/*.h)
  # add target executable
  add_executable(${DIM}_${TYPE}_${CMD}
    src/${DIM}/${TYPE}/${CMD}_cmd.cpp
    ${OTHER}
    ${HEADERS}
    ${COMMON}
    ${CUDA}
	${CUDA_COMMON})
  # add CUDA objects to the target
  if(CUDA_FOUND)
    target_link_libraries(${DIM}_${TYPE}_${CMD} ${CUDA_LIBRARIES})
  endif()
  # fix Intel OpenMP runtime location
  if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_custom_command(TARGET ${DIM}_${TYPE}_${CMD}
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL}
        -add_rpath /opt/intel/lib
        -change    libiomp5.dylib @rpath/libiomp5.dylib
        ${DIM}_${TYPE}_${CMD})
  endif()
endfunction()

# list all built PET commands using pet_target(subproj cmdname), where main
# source file is in src/subproj/cmdname_cmd.cpp

pet_target(1d toy gauss)

pet_target(2d barrel matrix)
pet_target(2d barrel phantom)
pet_target(2d barrel reconstruction)
pet_target(2d barrel lm_reconstruction)

pet_target(2d strip phantom)
pet_target(2d strip reconstruction)
pet_target(2d strip kernel_monte_carlo)

pet_target(2d toy phantom)

pet_target(2d barrel describe_scanner)

set(CL_2D_OPTIONS 2d/barrel/options.cpp 2d/strip/options.cpp)
set(CL_3D_OPTIONS 3d/hybrid/options.cpp)

pet_target(3d hybrid matrix ${CL_2D_OPTIONS})
pet_target(3d hybrid phantom ${CL_2D_OPTIONS})
pet_target(3d hybrid reconstruction ${CL_2D_OPTIONS})
pet_target(3d hybrid sensitivity ${CL_2D_OPTIONS})

pet_target(3d tool crop ${CL_2D_OPTIONS} ${CL_3D_OPTIONS})
pet_target(3d tool psf ${CL_2D_OPTIONS} ${CL_3D_OPTIONS})
pet_target(3d tool convert ${CL_2D_OPTIONS} ${CL_3D_OPTIONS})

if(Boost_FOUND)
  pet_target(2d barrel geometry)
endif()

# custom documentation and config files to be shown in project tree
file(GLOB EXTRA scripts/* cmake/* config/* phantoms/* src/3d/gate_input/* src/3d/full/* src/2d/gate/* testbed/*.py testbed/wscript)
add_custom_target(extra_mlem SOURCES
  README.md
  Doxyfile
  DoxygenLayout.xml
  .clang-format
  .gitignore
  .gitattributes
  ${EXTRA})

# tests
file(GLOB_RECURSE TESTS src/*_test.cpp)
include_directories(lib/catch/include)
if(POLICY CMP0037)
  cmake_policy(SET CMP0037 OLD) # allow "test" target
endif()
add_executable(utest src/util/test.cpp ${TESTS} ${COMMON})
set_target_properties(utest PROPERTIES EXCLUDE_FROM_ALL 1)

# clang doesn't print colored diagnostics when invoked from Ninja
if(UNIX AND
   CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
   CMAKE_GENERATOR STREQUAL "Ninja")
  add_definitions(-fcolor-diagnostics)
endif()

get_target_property(jpetmlem_include 3d_tool_convert INCLUDE_DIRECTORIES)

get_directory_property(has_parent PARENT_DIRECTORY)
if(has_parent)
  set(MLEM_INCLUDE_DIRS
    ${jpetmlem_include}
    PARENT_SCOPE
    )
endif()
