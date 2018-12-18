//
//  input.hpp
//  
//  This class handles reading in user options
//
//  Created by Jeremy Gibbs on 6/29/17.
//

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
        Output();
    
        // setter
        NcDim addDimension(std::string, int size=0);
        void addField(std::string, std::string, std::string, std::vector<NcDim>);
        void saveField1D(std::string, std::vector<double>&);
        void saveField1D(std::string, const std::vector<size_t>, double*);
        void saveField2D(std::string, const std::vector<size_t>,
                         std::vector<size_t>, std::vector<double>&);
};

#endif
