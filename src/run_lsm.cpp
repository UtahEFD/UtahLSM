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
    // Namelist input section
    std::string dirIn;
    int gridID;
    std::string timeBegin;
    std::string hInterval;
    std::string timeFinal;
    
    // Namelist patches section
    int procStartX;
    int procStartY;
    int nProcX;
    int nProcY;
    int nProcXin;
    
    // Namelist output section
    std::string dirOut;
    int deflate;
    std::vector<std::string> inVarList;
    std::vector<std::string> outVarList;
    
    // Aditional variables relating to namelist
    int nTimes;
    std::vector<std::string> dateStrings; 
        
    // Print a nice little welcome message
    std::cout << std::endl;
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#                     Welcome to UtahLSM                     #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#   A land surface model created at the University of Utah   #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    
    // Initialize the UtahLSM class
    UtahLSM utahlsm;
    
    // Initialize the input class
    Input input;
            
    // read in input section
    nError += input.getItem(&dirIn,     "input", "dirIn",     "");
    nError += input.getItem(&gridID,    "input", "gridID",    "");
    nError += input.getItem(&timeBegin, "input", "timeBegin", "");
    nError += input.getItem(&hInterval, "input", "hInterval", "");
    nError += input.getItem(&timeFinal, "input", "timeFinal", "");
    
    // read in patches section
    nError += input.getItem(&procStartX, "patches", "procStartX", "");
    nError += input.getItem(&procStartY, "patches", "procStartY", "");
    nError += input.getItem(&nProcX,     "patches", "nProcX", "");
    nError += input.getItem(&nProcY,     "patches", "nProcY", "");
    nError += input.getItem(&nProcXin,   "patches", "nProcXin", "");
    
    // read in output section
    nError += input.getItem(&dirOut,    "output", "dirOut", "");
    nError += input.getItem(&deflate,   "output", "deflate", "", 2);
    nError += input.getList(&inVarList, "output", "varList", "");
    
    if (nError) throw "There was an error reading the input file";
    
    std::cout<<"##############################################################"<<std::endl;
    
    return 0;
}
