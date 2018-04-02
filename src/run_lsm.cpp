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
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

namespace {
    namespace c = Constants;
}

int main () {
    
    // declare local variables
    bool first = true;
    int n_error = 0;
    double utc, atm_ws, net_r;
    double zeta_m=0,zeta_s=0,zeta_o=0,zeta_t=0;
    double ustar,flux_wT,flux_wq;
    NcVar t_var, z_var, ustar_var;
    NcVar flux_wT_var, flux_wq_var;
    NcVar soil_T_var, soil_q_var;
    
    // namelist time section
    double dt, utc_start;
    int nsteps;
    
    // namelist space section
    double z_o, z_t, z_m, z_s;
    int nsoilz;
    
    //namelist pressure section
    double atm_p;
    
    // namelist radiation section
    double albedo, emissivity, latitude, longitude;
    int julian_day, comp_rad;
        
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
    n_error += input.getItem(&comp_rad,   "radiation", "comp_rad",   "");
    
    // convert latitude and longitude into radians
    latitude  = latitude * c::pi / 180.0;
    longitude = longitude * c::pi / 180.0; 
    
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
    if (!comp_rad) n_error += input.getProf(&R_net, "metr", "R_net", nsteps);
    
    // modify soil levels to be negative
    std::transform(soil_z.begin(), soil_z.end(), soil_z.begin(),
          bind2nd(std::multiplies<double>(), -1.0));
    
    if (n_error) throw "There was an error reading the input file";
    
    // create output file
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"Creating output file"<<std::endl;
    NcFile outFile("lsm.nc", NcFile::replace);
    
    // define dimensions
    NcDim t_dim = outFile.addDim("t");
    NcDim z_dim = outFile.addDim("z",nsoilz);
    
    // define variables
    std::vector<NcDim> dim_vector;
    dim_vector.push_back(t_dim);
    dim_vector.push_back(z_dim);
    
    t_var       = outFile.addVar("time",    ncInt,   t_dim);
    z_var       = outFile.addVar("soil_z",  ncFloat, z_dim);
    ustar_var   = outFile.addVar("ustar",   ncFloat, t_dim);
    flux_wT_var = outFile.addVar("flux_wT", ncFloat, t_dim);
    flux_wq_var = outFile.addVar("flux_wq", ncFloat, t_dim);
    soil_T_var  = outFile.addVar("soil_T",  ncFloat, dim_vector);
    soil_q_var  = outFile.addVar("soil_q",  ncFloat, dim_vector);
    
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"Running UtahLSM"<<std::endl;;
    std::cout<<"##############################################################"<<std::endl;
    //nsteps = 1;
    for (int t=0; t<nsteps; ++t) {
        
        if (t>0) first = false;
        
        utc = utc_start + float(t+1)*dt;
        atm_ws = sqrt( pow(atm_u[t],2) + pow(atm_v[t],2) );
        
        std::cout<<"\rProcessing time: "<<utc<<std::flush;
        
        if (comp_rad) {
            net_r = 0;
        } else {
            net_r = R_net[t];
        }
        
        // Initialize the UtahLSM class
        UtahLSM utahlsm(first,dt,z_o,z_t,z_m,z_s,
                        atm_p,atm_ws,atm_T[t],atm_q[t],
                        nsoilz,soil_z,soil_T,soil_q,
                        porosity,psi_nsat,K_nsat,b,Ci,
                        julian_day,utc,latitude,longitude,
                        albedo,emissivity,net_r,comp_rad,
                        zeta_m,zeta_s,zeta_o,zeta_t,
                        ustar,flux_wT,flux_wq);
        
        // write output data
        const std::vector<size_t> index = {t};
        const std::vector<size_t> time_height_index = {static_cast<size_t>(t), 0};
        std::vector<size_t> time_height_size  = {1, nsoilz};
        t_var.putVar(index, utc);
        z_var.putVar(&soil_z[0]);
        ustar_var.putVar(index, ustar);
        flux_wT_var.putVar(index, c::rho_air*c::Cp_air*flux_wT);
        flux_wq_var.putVar(index, c::rho_air*c::Lv*flux_wq);
        soil_T_var.putVar(time_height_index, time_height_size, &soil_T[0]);
        soil_q_var.putVar(time_height_index, time_height_size, &soil_q[0]);
    }
    std::cout<<std::endl;
    std::cout<<"Finished!"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    return 0;
}