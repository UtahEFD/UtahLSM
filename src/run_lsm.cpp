//
//  run_lsm.cpp
//  
//  This program is an example that runs
//  UtahLSM in an offline mode
//
//  Created by Jeremy Gibbs on 10/30/17
//
#include "utah_lsm.hpp"
#include "constants.hpp"
#include "json.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <netcdf>
#include <fstream>
#include <iostream>
#include <iomanip>

using json = nlohmann::json;

using namespace netCDF;
using namespace netCDF::exceptions;

namespace {
    namespace c = Constants;
}

int main () {
    
    // declare local variables
    bool first = true;
    double ustar,flux_wT,flux_wq;
    NcVar t_var, z_var, ustar_var;
    NcVar flux_wT_var, flux_wq_var, Rnet_var;
    NcVar soil_T_var, soil_q_var, flux_gr_var;
    
    // input time
    int ntime; 
    double tstep;
    
    // input meteorology
    std::vector<double> atm_U;
    std::vector<double> atm_T;
    std::vector<double> atm_q;
    std::vector<double> atm_p;
    std::vector<double> R_net;
        
    // print a nice little welcome message
    std::cout << std::endl;
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#                     Welcome to UtahLSM                     #"<<std::endl;
    std::cout<<"#   A land surface model created at the University of Utah   #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    
    // read offline inpout files and de-serialize
    Input* inputOffline = new Input("inputOffline.json");

    inputOffline->getItem(ntime,"time","ntime");
    inputOffline->getItem(tstep,"time","tstep");
    inputOffline->getItem(atm_U,"data","atm_U");
    inputOffline->getItem(atm_T,"data","atm_T");
    inputOffline->getItem(atm_q,"data","atm_q");
    inputOffline->getItem(atm_p,"data","atm_p");
    inputOffline->getItem(R_net,"data","R_net");
    
    // initialize an instance of UtahLSM
    Input* inputLSM = new Input("inputLSM.json");
    flux_wq = 0;
    UtahLSM* utahlsm = new UtahLSM(inputLSM,ustar,flux_wT,flux_wq);
    std::cout<<"##############################################################"<<std::endl;
    
    struct timespec start, finish;
    double utc=0,elapsed;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int t=0; t<ntime; ++t) {

        // check if first time through
        if (t>0) first = false;

        // set time
        utc += tstep;
        std::cout<<"\rProcessing time: "<<std::setw(6)<<utc<<std::flush;
        
        // update user-specified fields
        utahlsm->updateFields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t]);
        
        // run the model
        utahlsm->run();
        

//        // write output data
//        const std::vector<size_t> index = {static_cast<unsigned long>(t)};
//        const std::vector<size_t> time_height_index = {static_cast<size_t>(t), 0};
//        std::vector<size_t> time_height_size  = {1, static_cast<unsigned long>(nsoilz)};
//        t_var.putVar(index, utc);
//        z_var.putVar(&soil_z[0]);
//        ustar_var.putVar(index, ustar);
//        flux_wT_var.putVar(index, c::rho_air*c::Cp_air*flux_wT);
//        flux_wq_var.putVar(index, c::rho_air*c::Lv*flux_wq);
//        flux_gr_var.putVar(index, flux_gr);
//        Rnet_var.putVar(index, net_r);
//        soil_T_var.putVar(time_height_index, time_height_size, &soil_T[0]);
//        soil_q_var.putVar(time_height_index, time_height_size, &soil_q[0]);
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout<<std::endl;
    std::cout<<"New flux: "<<flux_wq<<std::endl;
    std::cout<<"Finished in "<<elapsed<<" seconds!"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    

//    // create output file
//    std::cout<<"##############################################################"<<std::endl;
//    std::cout<<"Creating output file"<<std::endl;
//    NcFile outFile("lsm.nc", NcFile::replace);
//    
//    // define dimensions
//    NcDim t_dim = outFile.addDim("t");
//    NcDim z_dim = outFile.addDim("z",nsoilz);
//    
//    // define variables
//    std::vector<NcDim> dim_vector;
//    dim_vector.push_back(t_dim);
//    dim_vector.push_back(z_dim);
//    
//    t_var       = outFile.addVar("time",    ncInt,   t_dim);
//    z_var       = outFile.addVar("soil_z",  ncFloat, z_dim);
//    ustar_var   = outFile.addVar("ustar",   ncFloat, t_dim);
//    flux_wT_var = outFile.addVar("flux_wT", ncFloat, t_dim);
//    flux_wq_var = outFile.addVar("flux_wq", ncFloat, t_dim);
//    flux_gr_var = outFile.addVar("flux_gr", ncFloat, t_dim);
//    Rnet_var    = outFile.addVar("Rnet",    ncFloat, t_dim);
//    soil_T_var  = outFile.addVar("soil_T",  ncFloat, dim_vector);
//    soil_q_var  = outFile.addVar("soil_q",  ncFloat, dim_vector);
//
    return 0;
}
