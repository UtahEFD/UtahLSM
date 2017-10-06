# UtahLSM

This is an offline version of the University of Utah Land-Surface Model. 

UtahLSM version 2 is a modernized version of the original code described in Shingelton (2010).

## Building
Right now the makefile is limited to a gfortran version that supports Fortran 2003.

To build, type *make*

## Running
*./LSM case_name*

Here, *case_name* is a subfolder in the inputs directory. The folder must contain LSMinputs.txt, .ini files, and a timeseries file.
 
