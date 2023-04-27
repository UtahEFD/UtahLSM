#
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(lsm
  input.cpp
  logger.cpp
  matrix.cpp
  output.cpp
  radiation.cpp
  settings.cpp
  sfc.cpp
  sfc_most.cpp
  soil.cpp
  soil_brookscorey.cpp
  soil_campbell.cpp
  soil_vangenuchten.cpp
  soil_properties.cpp
  utah_lsm.cpp
)
