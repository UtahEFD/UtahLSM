#!/usr/bin/env python
import argparse
import numpy as np
import os
import sys
import time
from physics import Radiation, Soil, Surface
from util import constants, io, matrix

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
    def __init__(self,inputLSM, outputLSM, ustar, flux_wT, flux_wq, col_j=0, col_i=0):
        """Constructor method
        """
        # TODO: finish writing main init function
        # set the input and output fields
        self.input  = inputLSM
        self.output = outputLSM
        
        print("[UtahLSM: Setup] \t Reading input settings")
        # Input time section
        self.dt_seb = inputLSM.step_seb
        self.dt_dif = inputLSM.step_dif
        
        # Input grid section
        self.nx = inputLSM.nx
        self.ny = inputLSM.ny
        
        # Input length scale section
        self.z_o = self.input.z_o
        self.z_t = self.input.z_t
        self.z_m = self.input.z_m
        self.z_s = self.input.z_s
        
        # Input soil section
        self.nz             = self.input.nsoil
        self.soil_param     = self.input.param
        self.soil_model     = self.input.model
        
        print("[UtahLSM: Setup] \t Reading input data")
        self.soil_z    = self.input.soil_z
        self.soil_T    = self.input.soil_T
        self.soil_q    = self.input.soil_q
        self.soil_type = self.input.soil_type
        
        # Initialize new surface values for first run
        self.sfc_T_new = self.soil_T[0]
        self.sfc_q_new = self.soil_q[0]
        
        # Initialize history arrays for first run
        soil_T_last = self.soil_T
        soil_q_last = self.soil_q
        
        # Modify soil levels to be negative
        self.soil_z = -1*self.soil_z
        
        # Input radiation section
        self.comp_rad = self.input.comp_rad
        if (self.comp_rad==1):
            print("[UtahLSM: Setup] \t Creating radiation model")
            self.albedo     = self.input.albedo
            self.emissivity = self.input.emissivity
            self.latitude   = self.input.latitude
            self.longitude  = self.input.longitude
            self.julian_day = self.input.julian_day
            
            # convert latitude and longitude into radians
            self.latitude  = self.latitude * constants.pi / 180.0
            self.longitude = self.longitude * constants.pi / 180.0
            
            # Create radiation model
            self.rad = Radiation.get_model(1,self.input)
        else:
            print("[UtahLSM: Radiation] \t Using offline data, no model")
        
        # Create soil model
        print("[UtahLSM: Setup] \t Creating soil model")
        self.soil = Soil.get_model(self.soil_model,self.input)
        
        print("[UtahLSM: Setup] \t Creating surface model")
        # choose surface model
        self.sfc = Surface.get_model(1)

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
            'z':self.nz
        }
        self.output.set_dims(self.output_dims)

        # set reference to output fields
        self.output_fields = {
            'soil_z':self.soil_z,
            'soil_T':self.soil_T,
            'soil_q':self.soil_q,
            'soil_type':self.soil_type,
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
    # TODO: write update_fields function
    def update_fields(dt,atm_U, atm_T, atm_q, atm_p, R_net=0): pass
    
    # Run the model
    # TODO: write run function
    def run(): pass
        
    # Save output fields
    # TODO: write save function
    def save(): pass
    
    # Solve the surface energy budget
    # TODO: write solve_seb function
    def solve_seb(): pass
    
    # Solve the surface moisture budget
    # TODO: write solve_smb function
    def solve_smb(): pass
    
    # Solve the diffusion equation for soil heat
    # TODO: write solve_diffusion_heat function
    def solve_diffusion_heat(): pass
    
    # Solve the diffusion equation for soil moisture
    # TODO: write solve_diffusion_mois function
    def solve_diffusion_mois(): pass
    
    # Solve the diffusion equation for soil heat and moisture
    # 1=heat, 2=moisture
    # TODO: write solve_diffusion function
    def solve_diffusion(diff_type): pass
    
    # Compute the surface energy budget
    # TODO: write compute_seb function
    def compute_seb(sfc_T): pass
    
    # Compute the derivative of the surface energy budget
    # TODO: write compute_dseb function
    def compute_dseb(sfc_T): pass

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
    
    print("[UtahLSM: Info] \t Running offline for the %s case"%case)
    
    # create Output instance
    if not outf:
        outf='lsm_%s_py.nc'%case
    outputLSM = io.Output(outf)
    
    # grid information
    nx = inputLSM.nx
    ny = inputLSM.ny
    
    # get offline input data
    ntime = inputLSM.ntime
    tstep = inputLSM.tstep
    atm_U = inputLSM.atm_U
    atm_T = inputLSM.atm_T
    atm_q = inputLSM.atm_q
    atm_p = inputLSM.atm_p
    R_net = inputLSM.r_net
    
    # local fluxes to be modified by lsm
    ustar   = 0.0
    flux_wq = 0.0
    flux_wT = 0.0
    
    # Initialize an array of UtahLSM instances
    global_utah_lsm = np.empty((ny*nx)).astype(np.object_)
    
    # Fill array with lsm instances
    k = 0
    for j in range(0,ny):
        for i in range(0,nx):
            # create UtahLSM instance
            global_utah_lsm[k] = UtahLSM(inputLSM,outputLSM,ustar,flux_wq,flux_wT,j,i)
            k+=1
    
    # Loop through each time
    utc = 0
    for t in range(0,ntime):
        utc += tstep
        sys.stdout.write("\r[UtahLSM: Run] \t\t Running for time %05.2f of %05.2f"%(utc,ntime*tstep))
        sys.stdout.flush()
        
        # loop through each lsm instance
        k = 0
        #for j in range(0,ny):
        #    for i in range(0,nx):
                # update user-specified fields
                #global_utah_lsm[k].update_fields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t])
                
                # run the model
                #global_utah_lsm[k].run()
                
                # save output
                #global_utah_lsm[k].save(outputLSM)
    
    # time info
    t2 = time.time()
    tt = t2 - t1
    print("\n[UtahLSM: Run] \t\t Done! Completed in %0.2f seconds"%tt)
    print("##############################################################")