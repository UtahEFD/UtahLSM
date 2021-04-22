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

LINK_DIRECTORIES(/usr/local/lib)

add_executable(utahlsm_offline 
  utahlsm_offline.cpp
)

target_link_libraries(utahlsm_offline lsm)
target_link_libraries(utahlsm_offline ${LIBS})
target_link_libraries(utahlsm_offline ${NETCDF_LIBRARIES_CXX})
target_link_libraries(utahlsm_offline ${NETCDF_LIBRARIES_C})

add_executable(utahlsm_offline_f
  utahlsm_offline.f90
)
set_target_properties(utahlsm_offline_f PROPERTIES LINKER_LANGUAGE Fortran) 
enable_language(Fortran)
target_link_libraries(utahlsm_offline_f lsm)
target_link_libraries(utahlsm_offline_f ${LIBS})
target_link_libraries(utahlsm_offline_f ${NETCDF_LIBRARIES_CXX})
target_link_libraries(utahlsm_offline_f ${NETCDF_LIBRARIES_C})
