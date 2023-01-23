#!/usr/bin/env python
import argparse
import numpy as np
import os
import sys
import time
from physics import Radiation, Soil, Surface
from util import io, matrix

# custom error message for user case entry
class InvalidCase(Exception):
    pass

# land-surface model class
class UtahLSM:
    """This is the main UtahLSM class
    
    :param inputLSM: A handle to the :class:`util.io.Input`
    :param outputLSM: A handle to the :class:`util.io.Output`
    """
    
    # model class initialization
    def __init__(self,inputLSM, outputLSM, ustar, flux_wT, flux_wq, col_j, col_i):
        """Constructor method
        """
        
        # set the input and output fields
        self.input  = inputLSM
        self.output = outputLSM
        
        print("[UtahLSM: Setup] \t Reading input settings")
        # configuration variables
        nz         = self.input.nsoil
        soil_model = self.input.model
        
        print("[UtahLSM: Setup] \t Reading input data")
        # input fields
        self.soil_column    = np.empty(4).astype(np.object_)
        self.soil_column[0] = self.input.soil_z
        self.soil_column[1] = self.input.soil_T
        self.soil_column[2] = self.input.soil_q
        self.soil_column[3] = self.input.soil_type
        
        print("[UtahLSM: Setup] \t Creating radiation model")
        # choose radiation model
        self.rad = Radiation.get_model(1,self.input)
        
        print("[UtahLSM: Setup] \t Creating soil model")
        # choose surface model
        self.soil = Soil.get_model(soil_model,self.input)
        
        print("[UtahLSM: Setup] \t Creating surface model")
        # choose surface model
        #self.sfc = Surface.get_model(1,self.input)

        print("[UtahLSM: Setup] \t Creating additional fields")
        # initialize flux arrays
        self.obl = np.zeros(1)
        self.ust = np.zeros(1)
        self.shf = np.zeros(1)
        self.lhf = np.zeros(1)
        self.ghf = np.zeros(1)
        
        # Local time data
        self.first=True   # flag whether first time step or not
        self.step_count=0 # number of times the LSM has been called
        self.tstep=0      # current time step
        self.runtime=0    # current elapsed time
        self.utc=0        # current time in UTC

        print("[UtahLSM: Setup] \t Creating output file")    
        # set reference to output dimensions
        self.output_dims = {
            't':0,
            'z':nz
        }
        self.output.set_dims(self.output_dims)

        # set reference to output fields
        self.output_fields = {
            'soil_z':self.soil_column[0],
            'soil_T':self.soil_column[1],
            'soil_q':self.soil_column[2],
            'soil_type':self.soil_column[3],
            'ust':self.ust,
            'obl':self.obl,
            'shf':self.shf,
            'lhf':self.lhf,
            'ghf':self.ghf
        }
        self.output.set_fields(self.output_fields)

        # write initial data
        self.output.save(self.output_fields,0,0,initial=True)
        
    # Update atmospheric quantities prior to solving
    #def update_fields(dt,atm_U, atm_T, atm_q, atm_p, R_net=0):
        
        
        
    # void UtahLSM :: updateFields(double dt,double u,double T,double q,double p,double rad=0) {
    #     
    #     tstep = dt;
    #     atm_U = u;
    #     atm_T = T;
    #     atm_q = q;
    #     atm_p = p;
    #     R_net = rad;
    #     runtime += tstep;
    #     
    #     // Run radiation model and update time/date if needed
    #     if (comp_rad==1) {
    #         utc = std::fmod(runtime,86400);
    #         julian_day += int(runtime/86400);
    #         R_net = radiation->computeNet(julian_day,utc,soil_T[0]);
    #     }
    #     
    #     // Keep winds from being exactly zero
    #     if (atm_U==0) atm_U = 0.1;
    # }

# main program to run the LSM
if __name__ == "__main__":
    
    # let's time this thing
    t1 = time.time()

    # a nice welcome message
    print("##############################################################")
    print("#                                                            #")
    print("#                     Welcome to UtahLSM                     #")
    print("#   A land surface model created at the University of Utah   #")
    print("#       and the NOAA National Severe Storms Laboratory       #")
    print("#                                                            #")
    print("##############################################################")
    
    # get case from user
    parser = argparse.ArgumentParser(description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest='case', required=True,
                        action='store', type=str, help="Case name")
    parser.add_argument("-o", "--output", dest='outfile', 
                        action='store', type=str, help="Output file name")
    args = parser.parse_args()
    case = args.case
    outf = args.outfile
    
    # create Input instance
    try:
        if os.path.exists('../cases/%s/'%case):
            namelist    = '../cases/%s/lsm_namelist.json'%case
            initfile    = '../cases/%s/lsm_init.nc'%case 
            offlinefile = '../cases/%s/lsm_offline.nc'%case  
            inputLSM    = io.Input(namelist,initfile,offlinefile)
        else:
            raise InvalidCase('Error: The folder ../cases/%s does not exist.'%case)
    except InvalidCase as e:
        print(e)
        sys.exit(1)
    
    # create Output instance
    if not outf:
        outf='lsm_%s_py.nc'%case
    outputLSM = io.Output(outf)
    
    # create UtahLSM instance
    lsm = UtahLSM(inputLSM,outputLSM)
    
    # run the model
    #scm.run()
    
    # time info
    t2 = time.time()
    tt = t2 - t1
    print("\n[UtahLSM: Run] \t Done! Completed in %0.2f seconds"%tt)
    print("##############################################################")