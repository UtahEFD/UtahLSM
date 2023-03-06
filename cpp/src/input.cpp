/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2023 Jeremy A. Gibbs
 * Copyright (c) 2017–2023 Rob Stoll
 * Copyright (c) 2017–2023 Eric Pardyjak
 * Copyright (c) 2017–2023 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "input.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <netcdf>
#include <string>
#include <vector>

using namespace netCDF;
using namespace netCDF::exceptions;

// Constructor for Input class
Input :: Input(std::string input_file) {
    input = new NcFile(input_file, NcFile::read);
}

// Retrieve the requested dimension from the input netcdf file
void Input :: getDim(int& external, std::string name) {
    NcDim dim = input->getDim(name);
    external = dim.getSize();
}

// Retrieve the requested int from the input netcdf file
void Input :: getData(int& external, std::string name) {
    NcVar var = input->getVar(name);
    std::vector<size_t> var_index = {0};
    var.getVar(var_index,&external);
}

// Retrieve the requested double from the input netcdf file
void Input :: getData(double& external, std::string name) {
    NcVar var = input->getVar(name);
    std::vector<size_t> var_index = {0};
    var.getVar(var_index,&external);
}

// Retrieve the requested vector<int> from the input netcdf file
void Input :: getData(std::vector<int>& external, std::string name) {
    NcVar var = input->getVar(name);
    int countp = var.getDim(0).getSize();
    std::vector<size_t> var_index = {0};
    std::vector<size_t> var_size  = {static_cast<unsigned long>(countp)};
    var.getVar(var_index,var_size,&external[0]);
}

// Retrieve the requested vector<double> from the input netcdf file
void Input :: getData(std::vector<double>& external, std::string name) {
    NcVar var = input->getVar(name);
    int countp = var.getDim(0).getSize();
    std::vector<size_t> var_index = {0};
    std::vector<size_t> var_size  = {static_cast<unsigned long>(countp)};
    var.getVar(var_index,var_size,&external[0]);
}

//////////////////////////////////////////////////////////////
// C-style interface for compatibility with other languages //
//////////////////////////////////////////////////////////////

// C-style wrapper for the Input constructor
InputObject GetInput(char* input_file) {
    
    // Convert from char* to std::string
    std::string inputFile(input_file);

    // Remove trailing spaces sent from Fortran
    while(inputFile.size() && isspace(inputFile.back())) 
        inputFile.pop_back();
    
    // Return Input object
    Input *input = new Input(inputFile);

    return (InputObject)input;
}

// C-style wrapper for the getDim function
void GetDim(InputObject input, int* external, char* name) {
   
   // Get Input object
   Input* input_obj = (Input*)input;
   
   // Convert from char* to std::string
   std::string dimName(name);
    
   // Remove trailing spaces sent from Fortran
    while(dimName.size() && isspace(dimName.back())) 
        dimName.pop_back();
   
   // Get Data from Input object
    input_obj->getDim(*external,dimName);
    
}

// C-style wrapper for the getData function for an integer
void GetDataInt(InputObject input, int* external, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get Data from Input object
    input_obj->getData(*external,inputName);
    
    return; 
}

// C-style wrapper for the getData function for a double
void GetDataDbl(InputObject input, double* external, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get Data from Input object
    input_obj->getData(*external,inputName);
    
    return; 
}

// C-style wrapper for the getData function for a vector<double>
void GetDataIntArr(InputObject input, int external[], int* size, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Create a local vector
    std::vector<int> local(*size);
    
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get Data from Input object
    input_obj->getData(local,inputName);
    
    // Copy values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}

// C-style wrapper for the getData function for a vector<double>
void GetDataDblArr(InputObject input, double external[], int* size, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Create a local vector
    std::vector<double> local(*size);
    
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get Data from Input object
    input_obj->getData(local,inputName);
    
    // Copy values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}