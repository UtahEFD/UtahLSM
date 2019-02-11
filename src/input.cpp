//
//  input.cpp
//  
//  This class handles reading in user options
//  This is modified from version in MicroHH
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#include "input.hpp"
#include "json.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<cstring>

using json = nlohmann::json;

Input :: Input(std::string input_file) {
    
    // define json wrapper
    json input;
    
    // read input
    std::cout<<"[Input] \t Reading "<<input_file<<std::endl;
    readInputFile(input_file);
    
}

// Read in the namelist file
void Input :: readInputFile(std::string input_file) {
    
    // read file and de-serialize
    std::ifstream i(input_file);
    i >> input;
}

// Getter for integer
void Input :: getItem(int& external, std::string section, std::string name) {
    external = input[section][name].get<int>();
}

// Getter for double
void Input :: getItem(double& external, std::string section, std::string name) {
    external = input[section][name].get<double>();
}

// Getter for vector of ints
void Input :: getItem(std::vector<int>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<int>>();
}

// Getter for vector of doubles
void Input :: getItem(std::vector<double>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<double>>();
}

// Getter for vector of strings
void Input :: getItem(std::vector<std::string>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<std::string>>();
}

// C-style functions

// Interface function to return an input object
InputObject GetInput(char* input_file) {
    
    std::string inputFile(input_file);
    
    // remove trailing spaces sent from fortran
    while(inputFile.size() && isspace(inputFile.back())) 
        inputFile.pop_back();
    
    Input *input = new Input(inputFile);
    return (InputObject)input;
}

// Interface function to retrieve an int from an input object
void GetItemInt(InputObject input, int* external,char* section, char* name) {
    
    // Get input object
    Input* input_obj = (Input*)input;
        
    // Convert incoming char to string
    std::string inputSection(section);
    std::string inputName(name);
    
    // remove trailing spaces sent from fortran
    while(inputSection.size() && isspace(inputSection.back())) 
        inputSection.pop_back();
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from input object
    input_obj->getItem(*external,inputSection,inputName);
    
    return; 
}

// Interface function to retrieve a double from an input object
void GetItemDbl(InputObject input, double* external,char* section, char* name) {
    
    // Get input object
    Input* input_obj = (Input*)input;
        
    // Convert incoming char to string
    std::string inputSection(section);
    std::string inputName(name);
    
    // remove trailing spaces sent from fortran
    while(inputSection.size() && isspace(inputSection.back())) 
        inputSection.pop_back();
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from input object
    input_obj->getItem(*external,inputSection,inputName);
    
    return; 
}

// Interface function to retrieve a double array from an input object
void GetItemDblArr(InputObject input, double external[], int* size, char* section, char* name) {
    
    // Get input object
    Input* input_obj = (Input*)input;
        
    // create a local vector
    std::vector<double> local(*size);
    
    // Convert incoming char to string
    std::string inputSection(section);
    std::string inputName(name);
    
    // remove trailing spaces sent from fortran
    while(inputSection.size() && isspace(inputSection.back())) 
        inputSection.pop_back();
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from input object
    input_obj->getItem(local,inputSection,inputName);
    
    // Update values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}