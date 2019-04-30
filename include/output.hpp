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

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <string>
#include <vector>
#include <map>
#include <netcdf>

/**
 * This class handles saving output files.
 */

using namespace netCDF;
using namespace netCDF::exceptions;

class Output {
    
    private:
        // netCDF variables
        NcFile* outfile;
        std::map<std::string,NcVar> fields;
    public:
    
        // initializer
        Output(std::string);
        
        // setter
        NcDim addDimension(std::string, int size=0);
        void addField(std::string, std::string, std::string, std::vector<NcDim>);
        void saveField1D(std::string, const std::vector<size_t>, double*);
        void saveField2D(std::string, std::vector<double>&);
        void saveField2D(std::string, const std::vector<size_t>,
                         std::vector<size_t>, std::vector<double>&);
};


// C-style functions
typedef void * OutputObject;

extern "C" {
   OutputObject GetOutput(char*);
}

#endif
