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

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <map>
#include <netcdf>
#include <string>
#include <vector>

using namespace netCDF;
using namespace netCDF::exceptions;

/**
 * Class for managing output files.
 * 
 * This class is responsible for creating an output file, 
 * managing which fields to save, and generating their metadata.
 */
class Output {
    
    private:
        
        NcFile* outfile;                    ///< netcdf output file
        std::map<std::string,NcVar> fields; ///< map linking name with variable
    
    public:
    
        /**
         * Constructs an Output object.
         *
         * @param[in] output_file name of output file.
         */
        Output(std::string output_file);
        
        /**
         * Adds a dimension to the output file.
         *
         * @param[in] name name of the dimension
         * @param[in] size optional size of the dimension
         * @return    reference to the created NcDim
         */
        NcDim addDimension(std::string name, int size=0);
        
        /**
         * Adds a field to the output file.
         *
         * @param[in] name name of the field
         * @param[in] units units of the field
         * @param[in] long_name description of the field
         * @param[in] dims dimensions of the field
         */
        void addField(std::string name, std::string units,
                      std::string long_name, std::vector<NcDim> dims);
        
        /**
         * Saves a scalar double field to the output file.
         *
         * @param[in] name name of the field
         * @param[in] index data index of where to save field
         * @param[in] data scalar data to save
         */
        void saveFieldScalar(std::string name, const std::vector<size_t> index,
                             double* data);
        
        
        /**
         * Saves a 1D vector<double> field to the output file.
         *
         * @param[in] name name of the field
         * @param[in] data vector data to save
         */
        void saveFieldVector(std::string name, std::vector<double>& data);

        /**
         * Saves a 2D+ vector<double> field to the output file.
         *
         * @param[in] name name of the field
         * @param[in] index data index of where to save field
         * @param[in] size size in each dimension of the field
         * @param[in] data vector data to save
         */
        void saveFieldVector(std::string name, const std::vector<size_t> index,
                             std::vector<size_t> size, std::vector<double>& data);
};


/** 
 * C-style interface for compatibility with other languages.
 */
typedef void * OutputObject; ///< Pointer representing Output object 

extern "C" {
    /**
     * C-style wrapper for the Output constructor.
     *
     * @param[in] output_file name of output file
     */
    OutputObject GetOutput(char* output_file);
}

#endif