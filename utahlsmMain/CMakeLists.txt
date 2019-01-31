include_directories(${CMAKE_SOURCE_DIR}/include)

LINK_DIRECTORIES(/usr/local/lib)

add_executable(utahlsm 
  utahlsm.cpp
  )

target_link_libraries(utahlsm lsm)
target_link_libraries(utahlsm ${LIBS})
target_link_libraries(utahlsm ${NETCDF_LIBRARIES_CXX})
target_link_libraries(utahlsm ${NETCDF_LIBRARIES_C})


