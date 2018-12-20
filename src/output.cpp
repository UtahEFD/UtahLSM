//
//  input.cpp
//  
//  This class handles reading in user options
//  This is modified from version in MicroHH
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#include "utah_lsm.hpp"
#include "output.hpp"
#include "constants.hpp"
#include <iostream>

using namespace netCDF;
using namespace netCDF::exceptions;

namespace {
    namespace c = Constants;
}

Output :: Output(std::string output_file) {
    
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

void Output :: saveField1D(std::string name, const std::vector<size_t> index,
                           double* data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(index, data);
}

void Output :: saveField2D(std::string name, const std::vector<size_t> index,
                           std::vector<size_t> size, std::vector<double>& data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(index, size, &data[0]);
}

void Output :: saveField2D(std::string name, std::vector<double>& data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(&data[0]);
}
