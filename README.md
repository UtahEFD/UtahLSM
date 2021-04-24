# UtahLSM
University of Utah Land Surface Model



## Overview
UtahLSM is a lightweight, fast-response land-surface model suitable for use in large-eddy simulation (LES) codes.

More documentation to come

## Compiling

Compiling should be relatively quick and painless. To compile:

````
mkdir build
cd build
cmake -DNETCDF_DIR=/path/to/netcdf-c/include -DNETCDF_CXX_DIR=/path/to/netcdf-cxx/include ..
make
````
The executables will be placed in the `bin` folder at the top level of file hierarchy, and be symlinked in each case folder.
