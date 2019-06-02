/*
 * UtahLSM
 * 
 * Copyright (c) 2019 Jeremy A. Gibbs
 * Copyright (c) 2019 Pete Willemsen
 * Copyright (c) 2019 Rob Stoll
 * Copyright (c) 2019 Eric Pardyjak
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <netcdf>
#include <vector>

#include "constants.hpp"
#include "input.hpp"
#include "json.hpp"
#include "output.hpp"
#include "utah_lsm.hpp"

using json = nlohmann::json;

using namespace netCDF;
using namespace netCDF::exceptions;

namespace {
    namespace c = constants;
}

int main () {
    
    // print a nice little welcome message
    std::cout << std::endl;
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#                     Welcome to UtahLSM                     #"<<std::endl;
    std::cout<<"#   A land surface model created at the University of Utah   #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    
    // declare local variables
    std::vector<double> ustar(1,0.0);
    std::vector<double> flux_wT(1,0.0);
    std::vector<double> flux_wq(1,0.0);
    
    // input time
    int ntime; 
    double tstep;
    
    // input grid
    int nx, ny;
    
    // input meteorology
    std::vector<double> atm_U;
    std::vector<double> atm_T;
    std::vector<double> atm_q;
    std::vector<double> atm_p;
    std::vector<double> R_net;
    
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
    
    // initialize an instance of UtahLSM input and output
    Input* inputLSM  = new Input("inputLSM.json");
    Output* outputLSM = new Output("lsm.nc");
    
    // Get grid information
    inputLSM->getItem(nx,"grid","nx");
    inputLSM->getItem(ny,"grid","ny");
    
    // Initialize UtahLSM instance
    UtahLSM* utahlsm = new UtahLSM(inputLSM,outputLSM,ustar,flux_wT,flux_wq,nx,ny);
  
    // set up time information
    struct timespec start, finish;
    double utc=0,elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    for (int t=0; t<ntime; ++t) {

        // set time
        utc += tstep;
        std::cout<<std::fixed<<"\r[UtahLSM] \t Running for time: "<<std::setprecision(2)<<utc<<std::flush;
        
        // update user-specified fields
        utahlsm->updateFields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t]);

        // run the model
        utahlsm->run();

        // save output
        utahlsm->save(outputLSM);
    }

    // compuet run time information
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout<<std::endl;
    std::cout<<"[UtahLSM] \t Finished in "<<elapsed<<" seconds!"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;

    return 0;
}