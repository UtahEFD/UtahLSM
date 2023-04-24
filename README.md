# UtahLSM
University of Utah Land Surface Model (UtahLSM), created at the University of Utah and the NOAA National Severe Storms Laboratory

## Overview
UtahLSM is a lightweight, fast-response land-surface model suitable for use in large-eddy simulation (LES) codes.

More documentation to come

## Compiling

Compiling should be relatively quick and painless. To compile:

````shell
mkdir build
cd build
cmake -DNETCDF_DIR=/path/to/netcdf-c/include -DNETCDF_CXX_DIR=/path/to/netcdf-cxx/include ..
make
````
The executables will be placed in the `/bin` folder at the top level of file hierarchy, and be symlinked in each case folder.

## Documentation

UtahLSM supports Doxygen documentation generation. After the `cmake` command above, optionally type:
````shell
make doc
````

The generated html and latex documentation will be placed in the `/doc` folder.
