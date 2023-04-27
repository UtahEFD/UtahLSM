/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2023 Jeremy A. Gibbs
 * Copyright (c) 2017–2023 Rob Stoll
 * Copyright (c) 2017–2023 Eric Pardyjak
 * Copyright (c) 2017–2023 Pete Willemsen
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
#include "logger.hpp"
#include "input.hpp"
#include "json.hpp"
#include "output.hpp"
#include "settings.hpp"
#include "utah_lsm.hpp"

using json = nlohmann::json;

using namespace netCDF;
using namespace netCDF::exceptions;

namespace {
    namespace c = constants;
}

int main (int argc, char *argv[]) {
    
    // print a nice little welcome message
    std::cout << std::endl;
    std::cout<<"##############################################################"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"#                     Welcome to UtahLSM                     #"<<std::endl;
    std::cout<<"#   A land surface model created at the University of Utah   #"<<std::endl;
    std::cout<<"#       and the NOAA National Severe Storms Laboratory       #"<<std::endl;
    std::cout<<"#                                                            #"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;
    
    // process command line arguments
    std::string case_name, outf_name;
    
    if (argc <= 1) {
        throw std::runtime_error("No case name mode specified");
    }
    else {
        // get case name
        case_name = argv[1];
        
        // get the name of the output file
        if (argc > 2) {
            outf_name = argv[2];
        } else {
            outf_name = "lsm_"+case_name+"_cpp.nc";
        }
    }
    std::cout<<"[UtahLSM: Info] \t Running offline for the "+case_name+" case"<<std::endl;
    
    // declare local variables
    double ustar,flux_wT,flux_wq;
    
    // input time
    int ntime; 
    double tstep;
    
    // input grid
    int nx, ny;
    
    // read offline input file
    std::string offline = "../cases/"+case_name+"/lsm_offline.nc";
    Input* inputOffline = new Input(offline);
        
    // get offline input data
    inputOffline->getDim(ntime,"t");
    inputOffline->getData(tstep,"tstep");
    std::vector<double> atm_U(ntime);
    std::vector<double> atm_T(ntime);
    std::vector<double> atm_q(ntime);
    std::vector<double> atm_p(ntime);
    std::vector<double> R_net(ntime);
    inputOffline->getData(atm_U,"atm_U");
    inputOffline->getData(atm_T,"atm_T");
    inputOffline->getData(atm_q,"atm_q");
    inputOffline->getData(atm_p,"atm_p");
    inputOffline->getData(R_net,"R_net");
    
    // initialize an instance of UtahLSM input and output
    std::string settings = "../cases/"+case_name+"/lsm_namelist.json";
    std::string initfile = "../cases/"+case_name+"/lsm_init.nc";

    Input* inputLSM       = new Input(initfile);    
    Settings* settingsLSM = new Settings(settings);
    Output* outputLSM     = new Output(outf_name);
    
    // Get grid information
    settingsLSM->getItem(nx,"grid","nx");
    settingsLSM->getItem(ny,"grid","ny");
    
    // Initialize a vector of UtahLSM instances
    std::vector<UtahLSM*> globalUtahLSM(ny*nx);
    
    // Fill array with LSM instances
    int k = 0;
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            globalUtahLSM[k] = new UtahLSM(settingsLSM,inputLSM,outputLSM,ustar,flux_wT,flux_wq,j,i);
            k++;
        }
    }
    
    // set up time information
    struct timespec start, finish;
    double utc=0,elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    for (int t=0; t<ntime; ++t) {

        // set time
        utc += tstep;
        std::cout<<std::fixed<<"\r[UtahLSM: Run] \t\t Running for time: "<<std::setprecision(2)<<utc<<" of "<<ntime*tstep<<std::flush;
        
        // loop through each LSM instance
        int k = 0;
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                
                // update user-specified fields
                globalUtahLSM[k]->updateFields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t]);
                
                // run the model
                globalUtahLSM[k]->run();
                
                // save output
                globalUtahLSM[k]->save(outputLSM);
                
                k++;
            }
        }
    }
    
    // compute run time information
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout<<std::endl;
    std::cout<<std::setprecision(4);
    std::cout<<"[UtahLSM: Run] \t\t Done! Finished in "<<elapsed<<" seconds!"<<std::endl;
    std::cout<<"##############################################################"<<std::endl;

    return 0;
}