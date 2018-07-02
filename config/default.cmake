# system specific settings for joinwrf
set(ENV{CC} cc)
set(ENV{CXX} c++)

set(USER_CXX_FLAGS "--sysroot -restrict -std=c++11") 
set(USER_CXX_FLAGS_RELEASE "-O3")
set(USER_CXX_FLAGS_DEBUG "-debug -g -check=conversions,stack,uninit -fp-stack-check -fp-trap=common -fp-trap-all=common") 

set(NETCDF_INCLUDE_DIR "/usr/local/netcdf/include/")
set(NETCDF_LIB_C       "/usr/local/netcdf/lib/libnetcdf.dylib")
set(NETCDF_LIB_CPP     "/usr/local/netcdf/lib/libnetcdf-cxx4.1.dylib")
set(HDF5_LIB_1         "/usr/local/hdf5/lib/libhdf5.a")
set(HDF5_LIB_2         "/usr/local/hdf5/lib/libhdf5_hl.a")
set(SZIP_LIB           "/usr/local/szlib/lib/libsz.a")

set(LIBS ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB})
set(INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=restrict)
