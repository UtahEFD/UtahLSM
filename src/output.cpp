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

#include "output.hpp"

#include <iostream>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

Output :: Output(std::string output_file) {
    
    // create output file
    std::cout<<"[Output] \t Creating "<<output_file<<std::endl;
    outfile = new NcFile(output_file, NcFile::replace);
}

NcDim Output :: addDimension(std::string name, int size) {
    
    if (size) {
        return outfile->addDim(name, size);
    } else {
        return outfile->addDim(name);
    }
}

void Output :: addField(std::string name, std::string units,
                        std::string long_name, std::vector<NcDim> dims) {
 
    NcVar var;

    var = outfile->addVar(name, ncDouble, dims);
    var.putAtt("units", units);
    var.putAtt("long_name", long_name);
    fields[name] = var;
}

void Output :: saveFieldScalar(std::string name, const std::vector<size_t> index,
                               double* data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(index, data);
}

void Output :: saveFieldVector(std::string name, const std::vector<size_t> index,
                           std::vector<size_t> size, std::vector<double>& data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(index, size, &data[0]);
}

void Output :: saveFieldVector(std::string name, std::vector<double>& data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(&data[0]);
}

//////////////////////////////////////////////////////////////
// C-style interface for compatibility with other languages //
//////////////////////////////////////////////////////////////

// Interface function to return an output object
OutputObject GetOutput(char* output_file) {
    
    std::string outputFile(output_file);
    
    // remove trailing spaces sent from fortran
    while(outputFile.size() && isspace(outputFile.back())) 
        outputFile.pop_back();
    
    Output* output = new Output(outputFile);
    return (OutputObject)output;
}