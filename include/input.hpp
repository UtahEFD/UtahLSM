//
//  input.hpp
//  
//  This class handles reading in user options
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include "json.hpp"

using json = nlohmann::json;

/**
 * This class handles reading input files.
 */

class Input {
    
    private:
        // json wrapper
        json input;
    
        // function to read file
        void readInputFile(std::string);
    
    public:
    
        // initializer
        Input(std::string);
    
        // getters
        void getItem(int&, std::string, std::string);
        void getItem(double&, std::string, std::string);
        void getItem(std::vector<int>&, std::string, std::string);
        void getItem(std::vector<double>&, std::string, std::string);
};

#endif
