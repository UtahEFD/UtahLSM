#
#  utah_lsm cmake file
#  
#  These are settings for cmake to build utahlsm
#
#  Created by Jeremy Gibbs on 10/30/17

# set cmake version and config file location
cmake_minimum_required (VERSION 2.8.1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)

# install to the run directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# Set the default build type to RELEASE.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
endif()

# set compilers
#set(CMAKE_C_COMPILER "cc")
#set(CMAKE_CXX_COMPILER "c++")

# Start the project only after the system specific settings are loaded.
project(utahlsm)

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

# Only set the compiler flags when the cache is created
# to enable editing of the flags in the CMakeCache.txt file.
if(NOT HASCACHE)
  set(CMAKE_CXX_FLAGS ${USER_CXX_FLAGS} CACHE STRING
    "Flags used by the CXX-compiler during all build types." FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG ${USER_CXX_FLAGS_DEBUG} CACHE STRING
    "Flags used by the CXX-compiler during debug builds." FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE ${USER_CXX_FLAGS_RELEASE} CACHE STRING
    "Flags used by the CXX-compiler during release builds." FORCE)

  message(STATUS "Build Type: " ${CMAKE_BUILD_TYPE})
  set(HASCACHE TRUE CACHE BOOL "CMakeCache.txt created." FORCE)
  # Make sure that ccmake only contains build type.
  mark_as_advanced(HASCACHE)
  mark_as_advanced(CMAKE_INSTALL_PREFIX)
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
add_subdirectory(utahlsmMain)

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
