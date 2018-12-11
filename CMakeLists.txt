# CMake file for Resummino.
#
# Copyright 2008-2010 Jonathan Debove.
# Copyright 2011-2013 David R. Lamprea.
# Copyright 2011-2013 Marcel Rothering.
#
# Licensed under the terms of the EUPL version 1.1 or later.
# See the LICENCE file for more information.

cmake_minimum_required(VERSION 2.6)
project(resummino C CXX Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
include_directories(${CMAKE_SOURCE_DIR}/include)
link_directories(${CMAKE_SOURCE_DIR}/lib)

# Boost headers.
message(STATUS "Checking Boost headers")

set(BOOST_REQUIRED_HEADERS
  boost/algorithm/string/classification.hpp
  boost/algorithm/string/join.hpp
  boost/algorithm/string/predicate.hpp
  boost/algorithm/string/split.hpp
  boost/lexical_cast.hpp
)

foreach(headername ${BOOST_REQUIRED_HEADERS})
  string(TOUPPER "LIB_${headername}" var_headername)
  find_path(${var_headername} ${headername})
  if(NOT ${var_headername})
    message(FATAL_ERROR "Required Boost header ${headername} was not found.")
  else()
    include_directories(${${var_headername}})
  endif()
endforeach()

message(STATUS "Checking Boost headers - done")

# LHAPDF library and headers.
message(STATUS "Checking LHAPDF library and headers")

if(LHAPDF)
  file(TO_CMAKE_PATH "${LHAPDF}/include" LHAPDF_INCLUDE_DIR)
  file(TO_CMAKE_PATH "${LHAPDF}/lib" LHAPDF_LIBDIR)
endif()

find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h)
if(NOT LHAPDF_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find headers for library LHAPDF.")
else()
  include_directories(${LHAPDF_INCLUDE_DIR})
endif()

find_library(LHAPDF_LIB LHAPDF PATHS ${LHAPDF_LIBDIR})
if(NOT LHAPDF_LIB)
  message(FATAL_ERROR "Could not find library LHAPDF.")
else()
  message(STATUS "LHAPDF found: ${LHAPDF_LIB}")
endif()

message(STATUS "Checking LHAPDF library and headers - done")

# GSL library.
message(STATUS "Checking GNU Scientific Library")

set(GSL_REQUIRED_HEADERS
  gsl/gsl_blas.h
  gsl/gsl_math.h
  gsl/gsl_vector.h
  gsl/gsl_matrix.h
  gsl/gsl_integration.h
  gsl/gsl_monte_vegas.h
  gsl/gsl_multifit_nlin.h
)

foreach(headername ${GSL_REQUIRED_HEADERS})
  string(TOUPPER "LIB_${headername}" var_headername)
  find_path(${var_headername} ${headername})
  if(NOT ${var_headername})
    message(FATAL_ERROR "Required GNU Scientific Library header ${headername} was not found.")
  else()
    include_directories(${${var_headername}})
  endif()
endforeach()

find_library(GSL_LIB gsl)
find_library(CBLAS_LIB NAMES cblas gslcblas)
if(NOT (GSL_LIB AND CBLAS_LIB))
  message(FATAL_ERROR "Could not find GNU Scientific Library or CBLAS library.")
else()
  message(STATUS "GSL found: ${GSL_LIB}")
  message(STATUS "CBLAS found: ${CBLAS_LIB}")
endif()

message(STATUS "Checking GNU Scientific Library - done")

add_subdirectory(lib/ff)
add_subdirectory(lib/ql)
add_subdirectory(src)
