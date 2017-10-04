# UtahLSM

This is an offline version of the University of Utah Land-Surface Model. 

UtahLSM is modernized version of that presented in Shingelton (2010).

## Building
Right now the makefile is limited to a gfortran version that supports Fortran 2003.

To build: make

## Running
./LSM case_name

Here, case_name is a subfolder in the inputs directory. The folder must contain LSMinputs.txt, .ini files, and a timeseries file.
 
