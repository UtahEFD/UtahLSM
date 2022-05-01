#
# UtahLSM
# 
# Copyright (c) 2017–2022 Jeremy A. Gibbs
# Copyright (c) 2017–2022 Rob Stoll
# Copyright (c) 2017–2022 Eric Pardyjak
# Copyright (c) 2017–2022 Pete Willemsen
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

file( GLOB CASES ${CMAKE_SOURCE_DIR}/cases/* )
set(COUNTER 1)
foreach( case_folder ${CASES} )
  ADD_CUSTOM_TARGET(symlink_c++_${COUNTER} ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/utahlsm_offline ${case_folder}/utahlsm_offline)
  MATH(EXPR COUNTER "${COUNTER}+1")
  endforeach( case_folder ${CASES} )

add_executable(utahlsm_offline_f
  utahlsm_offline.f90
)

set(COUNTER 1)
foreach( case_folder ${CASES} )
  ADD_CUSTOM_TARGET(symlink_f90_${COUNTER} ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/utahlsm_offline_f ${case_folder}/utahlsm_offline_f)
  MATH(EXPR COUNTER "${COUNTER}+1")
endforeach( case_folder ${CASES} )

set_target_properties(utahlsm_offline_f PROPERTIES LINKER_LANGUAGE Fortran) 
enable_language(Fortran)
target_link_libraries(utahlsm_offline_f lsm)
target_link_libraries(utahlsm_offline_f ${LIBS})
target_link_libraries(utahlsm_offline_f ${NETCDF_LIBRARIES_CXX})
target_link_libraries(utahlsm_offline_f ${NETCDF_LIBRARIES_C})
