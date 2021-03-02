# CMake file for Resummino.
#
# Copyright 2008-2010 Jonathan Debove.
# Copyright 2011-2014 David R. Lamprea.
# Copyright 2011-2014 Marcel Rothering.
#
# Licensed under the terms of the EUPL version 1.1 or later.
# See the LICENCE file for more information.

cmake_minimum_required(VERSION 2.6)
project(resummino C CXX Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
include_directories(${CMAKE_SOURCE_DIR}/include)

include(ExternalProject)

## compiles looptools 213 (configure script changed so that the folder is always calles lib)
ExternalProject_Add(lt
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/LoopTools-2.13
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/lib/LoopTools-2.13/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE})

# Boost headers
include(cmake/RequireBoostHeaders.cmake)

# LHAPDF
find_package(LHAPDF REQUIRED)
include_directories(${LHAPDF_INCLUDE_DIRS})

# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

add_subdirectory(src)
