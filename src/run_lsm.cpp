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
    double ustar,flux_wT,flux_wq;
    
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
    
    // read offline input file
    Input* inputOffline = new Input("inputOffline.json");

    // get offline input data
    inputOffline->getItem(ntime,"time","ntime");
    inputOffline->getItem(tstep,"time","tstep");
    inputOffline->getItem(atm_U,"data","atm_U");
    inputOffline->getItem(atm_T,"data","atm_T");
    inputOffline->getItem(atm_q,"data","atm_q");
    inputOffline->getItem(atm_p,"data","atm_p");
    inputOffline->getItem(R_net,"data","R_net");
    
    // initialize an instance of UtahLSM input
    Input* inputLSM = new Input("inputLSM.json");
    
    // initialize an instance of UtahLSM
    UtahLSM* utahlsm = new UtahLSM(inputLSM,ustar,flux_wT,flux_wq);
    
    struct timespec start, finish;
    double utc=0,elapsed;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int t=0; t<ntime; ++t) {

        // set time
        utc += tstep;
        std::cout<<std::fixed<<"\r[UtahLSM] \t Running for time: "<<std::setw(7)<<utc<<std::flush;
        
        // update user-specified fields
        utahlsm->updateFields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t]);
        
        // run the model
        utahlsm->run();
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout<<std::endl;
    std::cout<<"[UtahLSM] \t Finished in "<<elapsed<<" seconds!"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;

    return 0;
}
