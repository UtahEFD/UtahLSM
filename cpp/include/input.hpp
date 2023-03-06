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

#ifndef INPUT_HPP
#define INPUT_HPP

#include <netcdf>
#include <vector>

using namespace netCDF;
using namespace netCDF::exceptions;

/**
 * Class for managing input files.
 * 
 * This class is responsible for opening a supplied input file
 * and returning requested fields from that file.
 */
class Input {
    
    private:

        NcFile* input; ///< netcdf input file
    
    public:

        /**
         * Constructs an Input object.
         *
         * @param[in] input_file name of input file
         */
        Input(std::string input_file);        
        
        /**
         * Retrieves the requested integer dimension from the input file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external integer pointer to fill with requested data
         * @param[in]  name name of the requested field in the json file
         */
        void getDim(int& external, std::string name);
        
        /**
         * Retrieves the requested integer from the input file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external integer pointer to fill with requested data
         * @param[in]  name name of the requested field in the json file
         */
        void getData(int& external, std::string name);
        
        /**
         * Retrieves the requested double from the input file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external double pointer to fill with requested data
         * @param[in] name name of the requested field in the netcdf file
         */
        void getData(double& external, std::string name);
        
        /**
         * Retrieves the requested vector<int> from the input file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external vector<int> pointer to fill with requested data
         * @param[in]  name name of the requested field in the json file
         */
        void getData(std::vector<int>& external, std::string name);
        
        /**
         * Retrieves the requested vector<double> from the input file and 
         * places it in the supplied pointer.
         *
         * @param[out] external external vector<double> pointer to fill with requested data
         * @param[in]  name name of the requested field in the netcdf file
         */
        void getData(std::vector<double>& external, std::string name);
};

typedef void * InputObject; ///< Pointer representing Input object 

/** 
 * C-style interface for compatibility with other languages.
 */
extern "C" {
    /**
     * C-style wrapper for the Input constructor.
     *
     * @param[in] input_file name of input file
     * @param[in] input_type filetype of input file
     */
    InputObject GetInput(char* input_file);
    
    /**
     * C-style wrapper for the getDim function
     *
     * @param[in]  input Input object
     * @param[out] external external integer pointer to fill with requested data
     * @param[in]  name name of the requested dimension
     */
    void GetDim(InputObject input, int* external, char* name);
    
    /**
     * C-style wrapper for the getData function for an integer.
     *
     * @param[in]  input Input object
     * @param[out] external external integer pointer to fill with requested data
     * @param[in]  name name of the requested field in the json file
     */
    void GetDataInt(InputObject input, int* external, char* name);
    
    /**
     * C-style wrapper for the getData function for a double.
     *
     * @param[in] input Input object
     * @param[out] external external double pointer to fill with requested data
     * @param[in] name name of the requested field in the json file
     */
    void GetDataDbl(InputObject input, double* external, char* name);
    
    /**
     * C-style wrapper for the getData function for a vector<double>.
     *
     * @param[in]  input Input object
     * @param[out] external external array pointer to fill with requested data
     * @param[in]  size number of elements in external array
     * @param[in]  name name of the requested field in the json file
     */
    void GetDataIntArr(InputObject input, double external[], int* size, char* name);
    
    /**
     * C-style wrapper for the getData function for a vector<double>.
     *
     * @param[in]  input Input object
     * @param[out] external external array pointer to fill with requested data
     * @param[in]  size number of elements in external array
     * @param[in]  name name of the requested field in the json file
     */
    void GetDataDblArr(InputObject input, double external[], int* size, char* name);
}

#endif