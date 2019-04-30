/*
 * UtahLSM
 * 
 * Copyright (c) 2019 Jeremy A. Gibbs
 * Copyright (c) 2019 Pete Willemsen
 * Copyright (c) 2019 Rob Stoll
 * Copyright (c) 2019 Eric Pardyjak
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include "json.hpp"

using json = nlohmann::json;

/**
 * Class for managing input files.
 * This class is responsible for opening a supplied input file
 * and returning requested fields from that file.
 */
class Input {
    
    private:

        json input; ///< json wrapper of input file
        
        /**
         * Reads a json input file
         *
         * @param input_file name of input file.
         */
        void readInputFile(std::string input_file);
    
    public:

        /**
         * Constructs an Input object.
         *
         * @param input_file name of input file.
         */
        Input(std::string input_file);        
    
        /**
         * Retrieves the requested integer from the input file and 
         * places it in the supplied pointer.
         *
         * @param external external integer pointer to fill with requested data.
         * @param section name of the section where the requested item resides in the json file.
         * @param name name of the requested field in the json file.
         */
        void getItem(int& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested double from the input file and 
         * places it in the supplied pointer.
         *
         * @param external external double pointer to fill with requested data.
         * @param section name of the section where the requested item resides in the json file.
         * @param name name of the requested field in the json file.
         */
        void getItem(double& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested vector<int> from the input file and 
         * places it in the supplied pointer.
         *
         * @param external external vector<int> pointer to fill with requested data.
         * @param section name of the section where the requested item resides in the json file.
         * @param name name of the requested field in the json file.
         */
        void getItem(std::vector<int>& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested vector<double> from the input file and 
         * places it in the supplied pointer.
         *
         * @param external external vector<double> pointer to fill with requested data.
         * @param section name of the section where the requested item resides in the json file.
         * @param name name of the requested field in the json file.
         */
        void getItem(std::vector<double>& external, std::string section, std::string name);
        
        /**
         * Retrieves the requested vector<string> from the input file and 
         * places it in the supplied pointer.
         *
         * @param external external vector<string> pointer to fill with requested data.
         * @param section name of the section where the requested item resides in the json file.
         * @param name name of the requested field in the json file.
         */
        void getItem(std::vector<std::string>& external, std::string section, std::string name);
};

typedef void * InputObject; ///< Pointer representing input object 

/** 
 * C-style interface for compatibility with other languages.
 */
extern "C" {
    /**
     * C-style wrapper for the Input constructor.
     *
     * @param input_file name of input file.
     */
    InputObject GetInput(char* input_file);

    /**
     * C-style wrapper for the getItem function for an integer.
     *
     * @param input Input object.
     * @param external external integer pointer to fill with requested data.
     * @param section name of the section where the requested item resides in the json file.
     * @param name name of the requested field in the json file.
     */
    void GetItemInt(InputObject input, int* external, char* section, char* name);
    
    /**
     * C-style wrapper for the getItem function for a double.
     *
     * @param input Input object.
     * @param external external double pointer to fill with requested data.
     * @param section name of the section where the requested item resides in the json file.
     * @param name name of the requested field in the json file.
     */
    void GetItemDbl(InputObject input, double* external, char* section, char* name);
    
    /**
     * C-style wrapper for the getItem function for a vector<double>.
     *
     * @param input Input object.
     * @param external external array pointer to fill with requested data.
     * @param section name of the section where the requested item resides in the json file.
     * @param name name of the requested field in the json file.
     */
    void GetItemDblArr(InputObject input, double external[], int* size, char* section, char* name);
}

#endif