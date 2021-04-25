#
# UtahLSM
# 
# Copyright (c) 2021 Jeremy A. Gibbs
# Copyright (c) 2021 Rob Stoll
# Copyright (c) 2021 Eric Pardyjak
# Copyright (c) 2021 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(lsm
  input.cpp
  matrix.cpp
  most.cpp
  output.cpp
  radiation.cpp
  soil.cpp
  soil_brookscorey.cpp
  soil_campbell.cpp
  soil_vangenuchten.cpp
  soil_properties.cpp
  utah_lsm.cpp
)
