//
//  run_lsm.cpp
//  
//
//  Created by Jeremy Gibbs on 10/30/17
//
#include "input.hpp"
#include "utah_lsm.hpp"
#include <iostream>

int main () {
    
    // declare local variables
    int nError = 0;
    
    // namelist length section
    double z_o;
    double z_t;
    double z_m;
    double z_s;
        
    // Print a nice little welcome message
    std::cout << std::endl;
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#                     Welcome to UtahLSM                     #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#   A land surface model created at the University of Utah   #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    
    // Initialize the input class
    Input input;
            
    // read in input section
    nError += input.getItem(&z_o, "length", "z_o", "");
    nError += input.getItem(&z_t, "length", "z_t", "");
    nError += input.getItem(&z_m, "length", "z_m", "");
    nError += input.getItem(&z_s, "length", "z_s", "");
    
    if (nError) throw "There was an error reading the input file";
    
    std::cout<<"##############################################################"<<std::endl;
    
    // Initialize the UtahLSM class
    UtahLSM utahlsm;
    
    return 0;
}