# UtahLSM
University of Utah Land Surface Model

## Overview
UtahLSM is a lightweight, fast-response land-surface model suitable for use in large-eddy simulation (LES) codes.

More documentation to come

## Compiling

When compiling on CHPC, NetCDF is installed in different directories.  First, make sure the correct module is loaded.

Then, 

````
mkdir build
cd build
cmake -DNETCDF_DIR=/uufs/chpc.utah.edu/sys/installdir/netcdf-c/4.4.1-c7/include -DNETCDF_CXX_DIR=/uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-5.4.0g/include ..
````

