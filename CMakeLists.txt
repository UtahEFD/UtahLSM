#
# UtahLSM
# 
# Copyright (c) 2019 Jeremy A. Gibbs
# Copyright (c) 2019 Pete Willemsen
# Copyright (c) 2019 Rob Stoll
# Copyright (c) 2019 Eric Pardyjak
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#

# set cmake version and config file location
cmake_minimum_required (VERSION 3.1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)

# install to the run directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Start the project only after the system specific settings are loaded.
project(UtahLSM VERSION 1.0 DESCRIPTION "UtahLSM is a lightweight, fast-response land-surface model" LANGUAGES CXX)

#
# This code used C++ 11 and 14 features so enable these on CXX Flags
#
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#
# NETCDF is Required! Search for it.
#
SET(NETCDF_CXX "YES")
FIND_PACKAGE(NetCDF REQUIRED)
IF(NetCDF_FOUND)
  MESSAGE(STATUS "Found NetCDF: ${NETCDF_INCLUDES}: ${NETCDF_LIBRARIES}, ${NETCDF_INCLUDES_CXX}, ${NETCDF_LIBRARIES_CXX}")
  INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})
  INCLUDE_DIRECTORIES(${NETCDF_INCLUDES_CXX})
ELSE(NetCDF_FOUND)
  # Get system settings file
  set(SYST default)
  include(${SYST} OPTIONAL RESULT_VARIABLE SYSTINC)
ENDIF(NetCDF_FOUND)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(config/doxygen.conf.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
    COMMENT "Generating API documentation with Doxygen." VERBATIM
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc DESTINATION ${CMAKE_SOURCE_DIR})
else(DOXYGEN)
  message(STATUS "WARNING: Doxygen not found - Reference manual will not be created.")
endif()

# Print the C++ and CUDA compiler flags to the screen.
if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_RELEASE})
else()
  message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_DEBUG})
endif()

file(GLOB_RECURSE HEADERS "include/*")
file(GLOB_RECURSE SOURCES "src/*")

include_directories("include" SYSTEM ${INCLUDE_DIRS})
link_libraries(${LIBS})

add_subdirectory(src)
add_subdirectory(utahlsmOffline)