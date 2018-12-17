//
//  input.cpp
//  
//  This class handles reading in user options
//  This is modified from version in MicroHH
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#include "output.hpp"
#include "utah_lsm.hpp"
#include "constants.hpp"
#include <iostream>

using namespace netCDF;
using namespace netCDF::exceptions;

namespace {
    namespace c = Constants;
}

Output :: Output() {
    
    // create output file
    outFile = new NcFile("lsm.nc", NcFile::replace);
}

NcDim Output :: addDimension(std::string name, int size) {
    
    if (size) {
        return outFile->addDim(name, size);
    } else {
        return outFile->addDim(name);
    }
}

void Output :: addField(std::string name, std::string units,
                        std::string long_name, std::vector<NcDim> dims) {
 
    // add field
    NcVar var = outFile->addVar(name, ncDouble, dims);
    var.putAtt("units", units);
    var.putAtt("long_name", long_name);
    
    // map name to field for later access
    fields[name] = var;
}

void Output :: saveField1D(std::string name, std::vector<double>* data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(data);
}

void Output :: saveField1D(std::string name, const std::vector<size_t> index,
                           double* data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(index, data);
}

void Output :: saveField2D(std::string name, const std::vector<size_t> index,
                           std::vector<size_t> size, std::vector<double>* data) {
    
    // write output data
    NcVar var = fields[name];
    var.putVar(index, size, data);
}
