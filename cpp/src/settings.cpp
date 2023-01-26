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

#include "settings.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "json.hpp"

using json = nlohmann::json;

// Constructor for Input class
Settings :: Settings(std::string input_file) {
    
   // json input
   json input;
   readInputFile(input_file);
}

// Read the input file
void Settings :: readInputFile(std::string input_file) {
    
    // read file and de-serialize
    std::ifstream i(input_file);
    i >> input;
}

// Retrieve the requested integer from the input file
void Settings :: getItem(int& external, std::string section, std::string name) {
    external = input[section][name].get<int>();
}

// Retrieve the requested double from the input file
void Settings :: getItem(double& external, std::string section, std::string name) {
    external = input[section][name].get<double>();
}

// Retrieve the requested vector<int> from the input file
void Settings :: getItem(std::vector<int>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<int> >();
}

// Retrieve the requested vector<double> from the input file
void Settings :: getItem(std::vector<double>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<double> >();
}

// Retrieve the requested vector<string> from the input file
void Settings :: getItem(std::vector<std::string>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<std::string> >();
}

//////////////////////////////////////////////////////////////
// C-style interface for compatibility with other languages //
//////////////////////////////////////////////////////////////

// C-style wrapper for the Input constructor
SettingsObject GetInput(char* input_file) {
    
    // Convert from char* to std::string
    std::string inputFile(input_file);

    // Remove trailing spaces sent from Fortran
    while(inputFile.size() && isspace(inputFile.back())) 
        inputFile.pop_back();
    
    // Return Settings object
    Input *input = new Input(inputFile);

    return (SettingsObject)input;
}

// C-style wrapper for the getItem function for an integer
void GetItemInt(SettingsObject input, int* external, char* section, char* name) {
    
    // Get Settings object
    Settings* settings_obj = (Settings*)input;
        
    // Convert from char* to std::string
    std::string inputSection(section);
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputSection.size() && isspace(inputSection.back())) 
        inputSection.pop_back();
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Settings object
    settings_obj->getItem(*external,inputSection,inputName);
    
    return; 
}

// C-style wrapper for the getItem function for a double
void GetItemDbl(SettingsObject input, double* external, char* section, char* name) {
    
    // Get Settings object
    Settings* settings_obj = (Settings*)input;
        
    // Convert from char* to std::string
    std::string inputSection(section);
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputSection.size() && isspace(inputSection.back())) 
        inputSection.pop_back();
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Settings object
    settings_obj->getItem(*external,inputSection,inputName);
    
    return; 
}

// C-style wrapper for the getItem function for a vector<double>
void GetItemDblArr(SettingsObject input, double external[], int* size, char* section, char* name) {
    
    // Get Settings object
    Settings* settings_obj = (Settings*)input;
        
    // Create a local vector
    std::vector<double> local(*size);
    
    // Convert from char* to std::string
    std::string inputSection(section);
    std::string inputName(name);
    
    // Remove trailing spaces sent from Fortran
    while(inputSection.size() && isspace(inputSection.back())) 
        inputSection.pop_back();
    while(inputName.size() && isspace(inputName.back())) 
        inputName.pop_back();
    
    // Get item from Settings object
    settings_obj->getItem(local,inputSection,inputName);
    
    // Copy values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}