include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(utahlsm 
  utahlsm.cpp
  )

target_link_libraries(utahlsm lsm)
target_link_libraries(utahlsm ${LIBS})
target_link_libraries(utahlsm ${NETCDF_LIBRARIES_CXX})
target_link_libraries(utahlsm ${NETCDF_LIBRARIES_C})
