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

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <vector>

#include "json.hpp"

using json = nlohmann::json;

/**
 * Class for managing json settings / data files.
 * 
 * This class is responsible for opening a supplied settings file
 * and returning requested fields from that file.
 */
class Settings {
    
    private:

        json settings; ///< json wrapper of settings file
        
        /**
         * Reads a json settings file.
         *
         * @param[in] settings_file name of settings file
         */
        void readSettingsFile(std::string settings_file);
    
    public:

        /**
         * Constructs a Settings object.
         *
         * @param[in] settings_file name of settings file
         */
        Settings(std::string settings_file);        
        
        /**
         * Retrieves the requested integer from the settings file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external integer pointer to fill with requested data
         * @param[in]  section name of the section where the requested item resides in the json file
         * @param[in]  name name of the requested field in the json file
         */
        void getItem(int& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested double from the settings file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external double pointer to fill with requested data
         * @param[in] section name of the section where the requested item resides in the json file
         * @param[in] name name of the requested field in the json file
         */
        void getItem(double& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested vector<int> from the settings file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external vector<int> pointer to fill with requested data
         * @param[in]  section name of the section where the requested item resides in the json file
         * @param[in]  name name of the requested field in the json file
         */
        void getItem(std::vector<int>& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested vector<double> from the settings file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external vector<double> pointer to fill with requested data
         * @param[in]  section name of the section where the requested item resides in the json file
         * @param[in]  name name of the requested field in the json file
         */
        void getItem(std::vector<double>& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested vector<string> from the settings file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external vector<string> pointer to fill with requested data
         * @param[in]  section name of the section where the requested item resides in the json file
         * @param[in]  name name of the requested field in the json file
         */
        void getItem(std::vector<std::string>& external, std::string section, std::string name);
};

typedef void * SettingsObject; ///< Pointer representing Settings object 

/** 
 * C-style interface for compatibility with other languages.
 */
extern "C" {
    /**
     * C-style wrapper for the Settings constructor.
     *
     * @param[in] settings_file name of settings file
     */
    SettingsObject GetSettings(char* settings_file);

    /**
     * C-style wrapper for the getItem function for an integer.
     *
     * @param[in]  settings Settings object
     * @param[out] external external integer pointer to fill with requested data
     * @param[in]  section name of the section where the requested item resides in the json file
     * @param[in]  name name of the requested field in the json file
     */
    void GetItemInt(SettingsObject settings, int* external, char* section, char* name);
    
    /**
     * C-style wrapper for the getItem function for a double.
     *
     * @param[in] settings Settings object
     * @param[out] external external double pointer to fill with requested data
     * @param[in] section name of the section where the requested item resides in the json file
     * @param[in] name name of the requested field in the json file
     */
    void GetItemDbl(SettingsObject settings, double* external, char* section, char* name);
    
    /**
     * C-style wrapper for the getItem function for a vector<double>.
     *
     * @param[in]  settings Settings object
     * @param[out] external external array pointer to fill with requested data
     * @param[in]  size number of elements in external array
     * @param[in]  section name of the section where the requested item resides in the json file
     * @param[in]  name name of the requested field in the json file
     */
    void GetItemIntArr(SettingsObject settings, int external[], int* size, char* section, char* name);
    
    /**
     * C-style wrapper for the getItem function for a vector<double>.
     *
     * @param[in]  settings Settings object
     * @param[out] external external array pointer to fill with requested data
     * @param[in]  size number of elements in external array
     * @param[in]  section name of the section where the requested item resides in the json file
     * @param[in]  name name of the requested field in the json file
     */
    void GetItemDblArr(SettingsObject settings, double external[], int* size, char* section, char* name);
}

#endif