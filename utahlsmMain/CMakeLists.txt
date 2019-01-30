include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(utahlsm 
  utahlsm.cpp
  )

target_link_libraries(utahlsm lsm)
target_link_libraries(utahlsm ${LIBS})