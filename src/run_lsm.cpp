//
//  run_lsm.cpp
//  
//
//  Created by Jeremy Gibbs on 10/30/17
//
#include "input.hpp"
#include "utah_lsm.hpp"
#include "constants.hpp"
#include <iostream>
#include <array>
#include <cmath>

int main () {
    
    // declare local variables
    int nError          = 0;
    const bool required = false;
    const bool optional = true;
    double utc, atm_ws;
    double phiM=0,psiM=0,psiM0=0,phiH=0,psiH=0,psiH0=0;
    double ustar,flux_wT,flux_wq;
    
    // declare namelist time section
    double dt, utc_start;
    int nsteps;
    
    // namelist space section
    double z_o, z_t, z_m, z_s;
    int nsoilz;
    
    //namelist pressure section
    double atm_p;
    
    // namelist radiation section
    double albedo, emissivity, latitude, longitude;
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
    nError += input.getItem(&dt,        "time", "dt",        "");
    nError += input.getItem(&utc_start, "time", "utc_start", "");
    nError += input.getItem(&nsteps,    "time", "nsteps",    "");
    
    // grab values from pressure section
    nError += input.getItem(&atm_p, "pressure", "p_o", "");
    
    // grab values from space section
    nError += input.getItem(&z_o,    "length", "z_o",    "");
    nError += input.getItem(&z_t,    "length", "z_t",    "");
    nError += input.getItem(&z_m,    "length", "z_m",    "");
    nError += input.getItem(&z_s,    "length", "z_s",    "");
    nError += input.getItem(&nsoilz, "length", "nsoilz", "");
    
    // grab values from radiation section
    nError += input.getItem(&albedo,     "radiation", "albedo",     "");
    nError += input.getItem(&emissivity, "radiation", "emissivity", "");
    nError += input.getItem(&latitude,   "radiation", "latitude",   "");
    nError += input.getItem(&longitude,  "radiation", "longitude",  "");
    nError += input.getItem(&julian_day, "radiation", "julian_day", "");
    
    if (nError) throw "There was an error reading the input file";
        
    // read in external soil data
    std::cout<<"Processing inputSoil.dat" << std::endl;
    double* soil_z   = new double[nsoilz];
    double* soil_T   = new double[nsoilz];
    double* soil_q   = new double[nsoilz];
    double* porosity = new double[nsoilz];
    double* psi_nsat = new double[nsoilz];
    double* K_nsat   = new double[nsoilz];
    double* b        = new double[nsoilz];
    double* Ci       = new double[nsoilz];
    
    nError += input.getProf(soil_z,   "soil", "soil_z",   nsoilz);
    nError += input.getProf(soil_T,   "soil", "soil_T",   nsoilz);
    nError += input.getProf(soil_q,   "soil", "soil_q",   nsoilz);
    nError += input.getProf(porosity, "soil", "porosity", nsoilz);
    nError += input.getProf(psi_nsat, "soil", "psi_nsat", nsoilz);
    nError += input.getProf(K_nsat,   "soil", "K_nsat",   nsoilz);
    nError += input.getProf(b,        "soil", "b",        nsoilz);
    nError += input.getProf(Ci,       "soil", "Ci",       nsoilz);
    
    if (nError) throw "There was an error reading the input file";
        
    // read in external atmospheric data
    std::cout<<"Processing inputMetr.dat" << std::endl;
    double* atm_u = new double[nsteps];
    double* atm_v = new double[nsteps];
    double* atm_T = new double[nsteps];
    double* atm_q = new double[nsteps];
    double* R_net = new double[nsteps];
    
    nError += input.getProf(atm_u, "metr", "atm_u", nsteps);
    nError += input.getProf(atm_v, "metr", "atm_v", nsteps);
    nError += input.getProf(atm_T, "metr", "atm_T", nsteps);
    nError += input.getProf(atm_q, "metr", "atm_q", nsteps);
    nError += input.getProf(R_net, "metr", "R_net",nsteps);
        
    std::cout<<"##############################################################"<<std::endl;
    
    if (nError) throw "There was an error reading the input file";
    
    std::cout<<"Running UtahLSM"<<std::endl;;
    std::cout<<"##############################################################"<<std::endl;
    nsteps=2;
    for (int t=0; t<nsteps; ++t) {
        
        utc = utc_start + float(t+1)*dt;
        atm_ws = sqrt( pow(atm_u[t],2) + pow(atm_v[t],2) );
                
        // Initialize the UtahLSM class
        UtahLSM utahlsm(z_o,z_t,z_s,z_m,
                        atm_p,atm_ws,atm_T[t],atm_q[t],
                        nsoilz,soil_z,soil_T,soil_q,
                        porosity,psi_nsat,K_nsat,b,Ci,
                        julian_day,utc,latitude,longitude,
                        albedo,emissivity,R_net[t],
                        &phiM,&psiM,&psiM0,&phiH,&psiH,&psiH0,
                        &ustar,&flux_wT,&flux_wq);
    }
        
    return 0;
}