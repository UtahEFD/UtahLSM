# UtahLSM
University of Utah Land Surface Model (UtahLSM), created at the University of Utah and the NOAA National Severe Storms Laboratory

## Overview
UtahLSM is a lightweight, fast-response land-surface model suitable for use in large-eddy simulation (LES) codes.

There is a C++ version of the model in the `/cpp` folder and a Python version of the model in the `/python` folder.

More documentation to come

## Compiling C++ Version

Compiling should be relatively quick and painless. To compile:

````shell
mkdir build
cd build
cmake -DNETCDF_DIR=/path/to/netcdf-c/include -DNETCDF_CXX_DIR=/path/to/netcdf-cxx/include ..
make
````
The executables will be placed in the `/bin` folder at the top level of `/cpp` file hierarchy.

## Documentation

### C++ Version

UtahLSM supports Doxygen documentation generation. After the `cmake` command above, optionally type:
````shell
make doc
````

The generated html and latex documentation will be placed in the `/doc` folder.

## Running

A series of offline test cases exist in the top-level `/cases` folder. Within each test case, run 
````shell
python makeInput.py
````
to generate the appropriate settings files and input data.

### C++

Running an offline test is easy. From the `/cpp` directory, run
````shell
./bin/utahlsm_offline -c [case name] -o [output filename, optional]
````
where `[case name]` is the folder name of a given case under `/cases`.

### Python

Running an offline test is easy. From the `/python` directory, run
````shell
python utahlsm.py -c [case name] -o [output filename, optional]
````
where `[case name]` is the folder name of a given case under `/cases`.