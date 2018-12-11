//
//  input.cpp
//  
//  This class handles reading in user options
//  This is modified from version in MicroHH
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#include "input.hpp"
#include "utah_lsm.hpp"
#include "constants.hpp"
#include "json.hpp"
#include<vector>
#include<iostream>
#include<fstream>

using json = nlohmann::json;

namespace {
    namespace c = Constants;
}

Input :: Input(std::string input_file) {
    
    // define json wrapper
    json input;
    
    // read input
    std::cout<<"[Input] \t Reading "<<input_file<<std::endl;
    readInputFile(input_file);
    
}

// Read in the namelist file
void Input :: readInputFile(std::string input_file) {
    
    // read file and de-serialize
    std::ifstream i(input_file);
    i >> input;
}

// Getter for integer
void Input :: getItem(int& external, std::string section, std::string name) {
    external = input[section][name].get<int>();
}

// Getter for double
void Input :: getItem(double& external, std::string section, std::string name) {
    external = input[section][name].get<double>();
}

// Getter for vector of ints
void Input :: getItem(std::vector<int>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<int>>();
}

// Getter for vector of doubles
void Input :: getItem(std::vector<double>& external, std::string section, std::string name) {
    external = input[section][name].get<std::vector<double>>();
}
