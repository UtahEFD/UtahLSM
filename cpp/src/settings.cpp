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

// Constructor for Settings class
Settings :: Settings(std::string settings_file) {
    
   // json settings
   json settings;
   readSettingsFile(settings_file);
}

// Read the settings file
void Settings :: readSettingsFile(std::string settings_file) {
    
    // read file and de-serialize
    std::ifstream i(settings_file);
    i >> settings;
}

// Retrieve the requested integer from the settings file
void Settings :: getItem(int& external, std::string section, std::string name) {
    external = settings[section][name].get<int>();
}

// Retrieve the requested double from the settings file
void Settings :: getItem(double& external, std::string section, std::string name) {
    external = settings[section][name].get<double>();
}

// Retrieve the requested vector<int> from the settings file
void Settings :: getItem(std::vector<int>& external, std::string section, std::string name) {
    external = settings[section][name].get<std::vector<int> >();
}

// Retrieve the requested vector<double> from the settings file
void Settings :: getItem(std::vector<double>& external, std::string section, std::string name) {
    external = settings[section][name].get<std::vector<double> >();
}

// Retrieve the requested vector<string> from the settings file
void Settings :: getItem(std::vector<std::string>& external, std::string section, std::string name) {
    external = settings[section][name].get<std::vector<std::string> >();
}

//////////////////////////////////////////////////////////////
// C-style interface for compatibility with other languages //
//////////////////////////////////////////////////////////////

// C-style wrapper for the Settings constructor
SettingsObject GetSettings(char* settings_file) {
    
    // Convert from char* to std::string
    std::string settingsFile(settings_file);

    // Remove trailing spaces sent from Fortran
    while(settingsFile.size() && isspace(settingsFile.back())) 
        settingsFile.pop_back();
    
    // Return Settings object
    Settings *settings = new Settings(settingsFile);

    return (SettingsObject)settings;
}

// C-style wrapper for the getItem function for an integer
void GetItemInt(SettingsObject settings, int* external, char* section, char* name) {
    
    // Get Settings object
    Settings* settings_obj = (Settings*)settings;
        
    // Convert from char* to std::string
    std::string settingsSection(section);
    std::string settingsName(name);
    
    // Remove trailing spaces sent from Fortran
    while(settingsSection.size() && isspace(settingsSection.back())) 
        settingsSection.pop_back();
    while(settingsName.size() && isspace(settingsName.back())) 
        settingsName.pop_back();
    
    // Get item from Settings object
    settings_obj->getItem(*external,settingsSection,settingsName);
    
    return; 
}

// C-style wrapper for the getItem function for a double
void GetItemDbl(SettingsObject settings, double* external, char* section, char* name) {
    
    // Get Settings object
    Settings* settings_obj = (Settings*)settings;
        
    // Convert from char* to std::string
    std::string settingsSection(section);
    std::string settingsName(name);
    
    // Remove trailing spaces sent from Fortran
    while(settingsSection.size() && isspace(settingsSection.back())) 
        settingsSection.pop_back();
    while(settingsName.size() && isspace(settingsName.back())) 
        settingsName.pop_back();
    
    // Get item from Settings object
    settings_obj->getItem(*external,settingsSection,settingsName);
    
    return; 
}

// C-style wrapper for the getItem function for a vector<double>
void GetItemDblArr(SettingsObject settings, double external[], int* size, char* section, char* name) {
    
    // Get Settings object
    Settings* settings_obj = (Settings*)settings;
        
    // Create a local vector
    std::vector<double> local(*size);
    
    // Convert from char* to std::string
    std::string settingsSection(section);
    std::string settingsName(name);
    
    // Remove trailing spaces sent from Fortran
    while(settingsSection.size() && isspace(settingsSection.back())) 
        settingsSection.pop_back();
    while(settingsName.size() && isspace(settingsName.back())) 
        settingsName.pop_back();
    
    // Get item from Settings object
    settings_obj->getItem(local,settingsSection,settingsName);
    
    // Copy values from local vector to external array
    for (int i=0;i<*size;i++) {
        external[i] = local[i];
    }
    
    return; 
}