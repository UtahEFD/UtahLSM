//
//  run_lsm.cpp
//  
//
//  Created by Jeremy Gibbs on 10/30/17
//
#include "input.hpp"
#include "utah_lsm.hpp"
#include <iostream>
#include <array>

int main () {
    
    // declare local variables
    int nError            = 0;
    const bool required = false;
    const bool optional = true;
    
    // declare namelist time section
    double dt;
    double utcStart;
    int nSteps;
    
    // namelist space section
    double z_o;
    double z_t;
    double z_m;
    double z_s;
    int nSoilZ;
    
    // namelist radiation section
    double albedo;    
    double emissivity;
    double latitude;  
    double longitude;
    int julian_day;
        
    // print a nice little welcome message
    std::cout << std::endl;
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#                     Welcome to UtahLSM                     #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#   A land surface model created at the University of Utah   #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    
    // initialize the input class
    Input input;
    
    // grab values from length section
    std::cout<<"Processing namelist.ini" << std::endl;
    nError += input.getItem(&dt,       "time", "dt",       "");
    nError += input.getItem(&utcStart, "time", "utcStart", "");
    nError += input.getItem(&nSteps,   "time", "nSteps",   "");
    
    // grab values from space section
    nError += input.getItem(&z_o,    "length", "z_o",    "");
    nError += input.getItem(&z_t,    "length", "z_t",    "");
    nError += input.getItem(&z_m,    "length", "z_m",    "");
    nError += input.getItem(&z_s,    "length", "z_s",    "");
    nError += input.getItem(&nSoilZ, "length", "nSoilZ", "");
    
    // grab values from radiation section
    nError += input.getItem(&albedo,     "radiation", "albedo",     "");
    nError += input.getItem(&emissivity, "radiation", "emissivity", "");
    nError += input.getItem(&latitude,   "radiation", "latitude",   "");
    nError += input.getItem(&longitude,  "radiation", "longitude",  "");
    nError += input.getItem(&julian_day, "radiation", "julian_day", "");
    
    if (nError) throw "There was an error reading the input file";
        
    // read in external soil data
    std::cout<<"Processing inputSoil.dat" << std::endl;
    double* zSoil      = new double[nSoilZ];
    double* tempSoil   = new double[nSoilZ];
    double* moisSoil   = new double[nSoilZ];
    double* porosity   = new double[nSoilZ];
    double* moisPotSat = new double[nSoilZ];
    double* hydCondSat = new double[nSoilZ];
    double* b          = new double[nSoilZ];
    double* Ci         = new double[nSoilZ];
    
    nError += input.getProf(zSoil,      "soil", "zSoil",      nSoilZ);
    nError += input.getProf(tempSoil,   "soil", "tempSoil",   nSoilZ);
    nError += input.getProf(moisSoil,   "soil", "moisSoil",   nSoilZ);
    nError += input.getProf(porosity,   "soil", "porosity",   nSoilZ);
    nError += input.getProf(moisPotSat, "soil", "moisPotSat", nSoilZ);
    nError += input.getProf(hydCondSat, "soil", "hydCondSat", nSoilZ);
    nError += input.getProf(b,          "soil", "b",          nSoilZ);
    nError += input.getProf(Ci,         "soil", "Ci",         nSoilZ);
    
    if (nError) throw "There was an error reading the input file";
        
    // read in external atmospheric data
    std::cout<<"Processing inputMetr.dat" << std::endl;
    double* uComp = new double[nSteps];
    double* vComp = new double[nSteps];
    double* pTemp = new double[nSteps];
    double* wvMix = new double[nSteps];
    double* rNet  = new double[nSteps];
    
    nError += input.getProf(uComp, "metr", "uComp", nSteps);
    nError += input.getProf(vComp, "metr", "vComp", nSteps);
    nError += input.getProf(pTemp, "metr", "pTemp", nSteps);
    nError += input.getProf(wvMix, "metr", "wvMix", nSteps);
    nError += input.getProf(rNet,  "metr", "rNet",  nSteps);
    
    std::cout<<"##############################################################"<<std::endl;
    
    if (nError) throw "There was an error reading the input file";
        
    // Initialize the UtahLSM class
    UtahLSM utahlsm;
    
    return 0;
}