include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(lsm 
  input.cpp
  matrix.cpp
  most.cpp
  output.cpp
  radiation.cpp
  soil.cpp
  utah_lsm.cpp
)