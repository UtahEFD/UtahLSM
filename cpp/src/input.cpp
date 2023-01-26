/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2022 Jeremy A. Gibbs
 * Copyright (c) 2017–2022 Rob Stoll
 * Copyright (c) 2017–2022 Eric Pardyjak
 * Copyright (c) 2017–2022 Pete Willemsen
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
void Input :: getItem(int& external, std::string name) {
    NcVar var = input->getVar(name);
    std::vector<size_t> var_index = {0};
    var.getVar(var_index,&external);
}

// Retrieve the requested double from the input netcdf file
void Input :: getItem(double& external, std::string name) {
    NcVar var = input->getVar(name);
    std::vector<size_t> var_index = {0};
    var.getVar(var_index,&external);
}

// Retrieve the requested vector<int> from the input netcdf file
void Input :: getItem(std::vector<int>& external, std::string name) {
   NcVar var = input->getVar(name);
   int countp = var.getDim(0).getSize();
   std::vector<size_t> var_index = {0};
   std::vector<size_t> var_size  = {static_cast<unsigned long>(countp)};
   var.getVar(var_index,var_size,&external[0]);
}

// Retrieve the requested vector<double> from the input netcdf file
void Input :: getItem(std::vector<double>& external, std::string name) {
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
InputObject GetInput(char* input_file, char* input_type) {
    
    // Convert from char* to std::string
    std::string inputFile(input_file);

    // Remove trailing spaces sent from Fortran
    while(inputFile.size() && isspace(inputFile.back())) 
        inputFile.pop_back();
    
    // Return Input object
    Input *input = new Input(inputFile);

    return (InputObject)input;
}

// C-style wrapper for the getItem function for an integer
void GetItemInt(InputObject input, int* external, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Input object
    input_obj->getItem(*external,inputName);
    
    return; 
}

// C-style wrapper for the getItem function for a double
void GetItemDbl(InputObject input, double* external, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Input object
    input_obj->getItem(*external,inputName);
    
    return; 
}

// C-style wrapper for the getItem function for a vector<double>
void GetItemIntArr(InputObject input, int external[], int* size, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Create a local vector
    std::vector<int> local(*size);
    
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Input object
    input_obj->getItem(local,inputName);
    
    // Copy values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}

// C-style wrapper for the getItem function for a vector<double>
void GetItemDblArr(InputObject input, double external[], int* size, char* name) {
    
    // Get Input object
    Input* input_obj = (Input*)input;
        
    // Create a local vector
    std::vector<double> local(*size);
    
    // Convert from char* to std::string
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Input object
    input_obj->getItem(local,inputName);
    
    // Copy values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}