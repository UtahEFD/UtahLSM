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
#include <cmath>
#include <vector>
#include <functional>

int main () {
    
    // declare local variables
    int n_error         = 0;
    const bool required = false;
    const bool optional = true;
    double utc, atm_ws;
    double zeta_m=0,zeta_s=0,zeta_o=0,zeta_t=0;
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
    n_error += input.getItem(&dt,        "time", "dt",        "");
    n_error += input.getItem(&utc_start, "time", "utc_start", "");
    n_error += input.getItem(&nsteps,    "time", "nsteps",    "");
    
    // grab values from pressure section
    n_error += input.getItem(&atm_p, "pressure", "p_o", "");
    
    // grab values from space section
    n_error += input.getItem(&z_o,    "length", "z_o",    "");
    n_error += input.getItem(&z_t,    "length", "z_t",    "");
    n_error += input.getItem(&z_m,    "length", "z_m",    "");
    n_error += input.getItem(&z_s,    "length", "z_s",    "");
    n_error += input.getItem(&nsoilz, "length", "nsoilz", "");
    
    // grab values from radiation section
    n_error += input.getItem(&albedo,     "radiation", "albedo",     "");
    n_error += input.getItem(&emissivity, "radiation", "emissivity", "");
    n_error += input.getItem(&latitude,   "radiation", "latitude",   "");
    n_error += input.getItem(&longitude,  "radiation", "longitude",  "");
    n_error += input.getItem(&julian_day, "radiation", "julian_day", "");
    
    if (n_error) throw "There was an error reading the input file";
    
    // read in external soil data
    std::cout<<"Processing inputSoil.dat" << std::endl;
    std::vector<double> soil_z;
    std::vector<double> soil_T;
    std::vector<double> soil_q;
    std::vector<double> porosity;
    std::vector<double> psi_nsat;
    std::vector<double> K_nsat;
    std::vector<double> b;
    std::vector<double> Ci;
        
    n_error += input.getProf(&soil_z,   "soil", "soil_z",   nsoilz);
    n_error += input.getProf(&soil_T,   "soil", "soil_T",   nsoilz);
    n_error += input.getProf(&soil_q,   "soil", "soil_q",   nsoilz);
    n_error += input.getProf(&porosity, "soil", "porosity", nsoilz);
    n_error += input.getProf(&psi_nsat, "soil", "psi_nsat", nsoilz);
    n_error += input.getProf(&K_nsat,   "soil", "K_nsat",   nsoilz);
    n_error += input.getProf(&b,        "soil", "b",        nsoilz);
    n_error += input.getProf(&Ci,       "soil", "Ci",       nsoilz);
        
    if (n_error) throw "There was an error reading the input file";
        
    // read in external atmospheric data
    std::cout<<"Processing inputMetr.dat" << std::endl;
    std::vector<double> atm_u;
    std::vector<double> atm_v;
    std::vector<double> atm_T;
    std::vector<double> atm_q;
    std::vector<double> R_net;
    
    n_error += input.getProf(&atm_u, "metr", "atm_u", nsteps);
    n_error += input.getProf(&atm_v, "metr", "atm_v", nsteps);
    n_error += input.getProf(&atm_T, "metr", "atm_T", nsteps);
    n_error += input.getProf(&atm_q, "metr", "atm_q", nsteps);
    n_error += input.getProf(&R_net, "metr", "R_net", nsteps);
    
    // modify soil levels to be negative
    std::transform(soil_z.begin(), soil_z.end(), soil_z.begin(),
          bind2nd(std::multiplies<double>(), -1.0));
    
    std::cout<<"##############################################################"<<std::endl;
    
    if (n_error) throw "There was an error reading the input file";
    
    std::cout<<"Running UtahLSM"<<std::endl;;
    std::cout<<"##############################################################"<<std::endl;
    //nsteps=10;
    for (int t=0; t<nsteps; ++t) {
        
        utc = utc_start + float(t+1)*dt;
        atm_ws = sqrt( pow(atm_u[t],2) + pow(atm_v[t],2) );
        
        std::cout<<"Processing time: "<<utc<<std::endl;
             
        // Initialize the UtahLSM class
        UtahLSM utahlsm(dt,z_o,z_t,z_m,z_s,
                        atm_p,atm_ws,atm_T[t],atm_q[t],
                        nsoilz,soil_z,soil_T,soil_q,
                        porosity,psi_nsat,K_nsat,b,Ci,
                        julian_day,utc,latitude,longitude,
                        albedo,emissivity,R_net[t],
                        zeta_m,zeta_s,zeta_o,zeta_t,
                        ustar,flux_wT,flux_wq);
    } 
    return 0;
}