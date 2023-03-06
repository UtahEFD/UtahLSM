# 
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

class SoilType(object):
    
    def __init__(self):
        
        # TODO: implement check on type
        
        self.b=0.0        # exponent (unitless)
        self.psi_sat=0.0  # saturation moisture potential (m)
        self.porosity=0.0 # saturated soil moisture 
        self.residual=0.0 # residual moisture (volume/volume)
        self.K_sat=0.0    # hydraulic conductivity (m/s)
        self.ci=0.0       # volumetric heat capacity (J/m^3/K)
    
    @staticmethod
    def get_properties(dataset,soil_type):
        if soil_type == 1:
            return Sand(dataset)
        if soil_type == 2:
            return LoamySand(dataset)
        if soil_type == 3:
            return SandyLoam(dataset)
        if soil_type == 4:
            return SiltyLoam(dataset)
        if soil_type == 5:
            return Loam(dataset)
        if soil_type == 6:
            return SandyClayLoam(dataset)
        if soil_type == 7:
            return SiltyClayLoam(dataset)
        if soil_type == 8:
            return ClayLoam(dataset)
        if soil_type == 9:
            return SandyClay(dataset)
        if soil_type == 10:
            return SiltyClay(dataset)
        if soil_type == 11:
            return Clay(dataset)
        if soil_type == 12:
            return Peat(dataset)
        if soil_type == 13:
            return B11(dataset)
        if soil_type == 14:
            return O12(dataset)
        if soil_type == 15:
            return O16(dataset)

# properties for sand (type=1)
class Sand(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 4.05
            self.psi_sat  = -0.121
            self.porosity = 0.395
            self.residual = 0.0000
            self.K_sat    = 1.76e-04
            self.ci       = 1470000
        
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 2.79
            self.psi_sat  = -0.023
            self.porosity = 0.339
            self.residual = 0.0000
            self.K_sat    = 1.60e-05
            self.ci       = 1470000
        
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 1.44
            self.psi_sat  = -0.160
            self.porosity = 0.437
            self.residual = 0.0200
            self.K_sat    = 5.83e-05
            self.ci       = 1470000

# properties for loamy sand (type=2)
class LoamySand(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 4.38
            self.psi_sat  = -0.090
            self.porosity = 0.410
            self.residual = 0.0000
            self.K_sat    = 1.56e-04
            self.ci       = 1410000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 4.26
            self.psi_sat  = -0.018
            self.porosity = 0.421
            self.residual = 0.0000
            self.K_sat    = 9.52e-06
            self.ci       = 1410000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 1.00
            self.psi_sat  = -0.206
            self.porosity = 0.437
            self.residual = 0.0350
            self.K_sat    = 1.70e-05
            self.ci       = 1410000

# properties for sandy loam (type=3)
class SandyLoam(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 4.90
            self.psi_sat  = -0.218
            self.porosity = 0.435
            self.residual = 0.0000
            self.K_sat    = 3.41e-05
            self.ci       = 1340000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 4.74
            self.psi_sat  = -0.032
            self.porosity = 0.434
            self.residual = 0.0000
            self.K_sat    = 6.19e-06
            self.ci       = 1340000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 81.00
            self.psi_sat  = -0.302
            self.porosity = 0.453
            self.residual = 0.0410
            self.K_sat    = 7.19e-06
            self.ci       = 1340000

# properties for silty loam (type=4)
class SiltyLoam(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 5.30
            self.psi_sat  = -0.786
            self.porosity = 0.485
            self.residual = 0.0000
            self.K_sat    = 7.20e-06
            self.ci       = 1270000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 5.33
            self.psi_sat  = -0.066
            self.porosity = 0.476
            self.residual = 0.0000
            self.K_sat    = 4.73e-06
            self.ci       = 1270000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 2.65
            self.psi_sat  = -0.401
            self.porosity = 0.463
            self.residual = 0.0270
            self.K_sat    = 1.89e-06
            self.ci       = 1270000

# properties for loam (type=5)
class Loam(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 5.39
            self.psi_sat  = -0.478
            self.porosity = 0.451
            self.residual = 0.0000
            self.K_sat    = 7.00e-06
            self.ci       = 1210000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 5.25
            self.psi_sat  = -0.047
            self.porosity = 0.439
            self.residual = 0.0000
            self.K_sat    = 5.12e-06
            self.ci       = 1210000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 3.97
            self.psi_sat  = -0.509
            self.porosity = 0.501
            self.residual = 0.0150
            self.K_sat    = 3.67e-06
            self.ci       = 1210000

# properties for sandy clay loam (type=6)
class SandyClayLoam(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 7.12
            self.psi_sat  = -0.299
            self.porosity = 0.420
            self.residual = 0.0000
            self.K_sat    = 6.30e-06
            self.ci       = 1180000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 6.77
            self.psi_sat  = -0.031
            self.porosity = 0.404
            self.residual = 0.0000
            self.K_sat    = 5.78e-06
            self.ci       = 1180000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 4.27
            self.psi_sat  = -0.594
            self.porosity = 0.398
            self.residual = 0.0680
            self.K_sat    = 1.19e-06
            self.ci       = 1180000

# properties for silty clay loam (type=7)
class SiltyClayLoam(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 7.75
            self.psi_sat  = -0.356
            self.porosity = 0.477
            self.residual = 0.0000
            self.K_sat    = 1.70e-06
            self.ci       = 1320000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 8.72
            self.psi_sat  = -0.060
            self.porosity = 0.464
            self.residual = 0.0000
            self.K_sat    = 4.11e-06
            self.ci       = 1320000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 3.13
            self.psi_sat  = -0.564
            self.porosity = 0.464
            self.residual = 0.0750
            self.K_sat    = 6.39e-07
            self.ci       = 1320000

# properties for clay loam (type=8)
class ClayLoam(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 8.52
            self.psi_sat  = -0.630
            self.porosity = 0.476
            self.residual = 0.0000
            self.K_sat    = 2.50e-06
            self.ci       = 1230000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 8.17
            self.psi_sat  = -0.041
            self.porosity = 0.465
            self.residual = 0.0000
            self.K_sat    = 4.45e-06
            self.ci       = 1230000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 4.13
            self.psi_sat  = -0.703
            self.porosity = 0.471
            self.residual = 0.0400
            self.K_sat    = 4.17e-07
            self.ci       = 1230000

# properties for sandy clay (type=9)
class SandyClay(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 10.40
            self.psi_sat  = -0.153
            self.porosity = 0.426
            self.residual = 0.0000
            self.K_sat    = 2.20e-06
            self.ci       = 1180000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 10.73
            self.psi_sat  = -0.027
            self.porosity = 0.406
            self.residual = 0.0000
            self.K_sat    = 7.12e-05
            self.ci       = 1180000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 5.65
            self.psi_sat  = -0.795
            self.porosity = 0.430
            self.residual = 0.1090
            self.K_sat    = 3.33e-07
            self.ci       = 1180000

# properties for silty clay (type=10)
class SiltyClay(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 10.40
            self.psi_sat  = -0.490
            self.porosity = 0.492
            self.residual = 0.0000
            self.K_sat    = 1.00e-06
            self.ci       = 1150000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 10.39
            self.psi_sat  = -0.045
            self.porosity = 0.468
            self.residual = 0.0000
            self.K_sat    = 3.43e-06
            self.ci       = 1150000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 4.48
            self.psi_sat  = -0.765
            self.porosity = 0.479
            self.residual = 0.0560
            self.K_sat    = 2.50e-07
            self.ci       = 1150000

# properties for clay (type=11)
class Clay(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 11.40
            self.psi_sat  = -0.405
            self.porosity = 0.482
            self.residual = 0.0000
            self.K_sat    = 1.30e-06
            self.ci       = 1090000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 10.55
            self.psi_sat  = -0.053
            self.porosity = 0.468
            self.residual = 0.0000
            self.K_sat    = 2.99e-06
            self.ci       = 1090000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 6.67
            self.psi_sat  = -0.856
            self.porosity = 0.475
            self.residual = 0.0900
            self.K_sat    = 1.67e-07
            self.ci       = 1090000

# properties for clay (type=12)
class Peat(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Clapp and Hornberger (1974)
        if (dataset==1):
            self.b        = 7.75
            self.psi_sat  = -0.356
            self.porosity = 0.863
            self.residual = 0.0000
            self.K_sat    = 8.00e-06
            self.ci       = 840000
       
        # Cosby et al. (1984)
        if (dataset==2):
            self.b        = 7.75
            self.psi_sat  = -0.356
            self.porosity = 0.863
            self.residual = 0.0000
            self.K_sat    = 8.00e-06
            self.ci       = 840000
       
        # Rawls and Brakensiek (1982)
        if (dataset==3):
            self.b        = 6.06
            self.psi_sat  = -0.356
            self.porosity = 0.863
            self.residual = 0.1763
            self.K_sat    = 8.00e-06
            self.ci       = 840000

# properties for B11 (type=13)
class B11(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Heinen, Bakker, Wosten (Cabauw-specific)
        if (dataset==4):
            self.b        = 9.35
            self.psi_sat  = -0.463
            self.porosity = 0.591
            self.residual = 0.01
            self.K_sat    = 7.30e-07
            self.ci       = 1090000

# properties for O12 (type=14)
class O12(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Heinen, Bakker, Wosten (Cabauw-specific)
        if (dataset==4):
            self.b        = 6.33
            self.psi_sat  = -1.14
            self.porosity = 0.561
            self.residual = 0.01
            self.K_sat    = 1.25e-07
            self.ci       = 1090000

# properties for O16 (type=15)
class O16(SoilType):
    
    def __init__(self,dataset):
        
        # initialize parent class
        super().__init__()
        
        # Heinen, Bakker, Wosten (Cabauw-specific)
        if (dataset==4):
            self.b        = 2.75
            self.psi_sat  = -1.03
            self.porosity = 0.889
            self.residual = 0.00
            self.K_sat    = 1.69e-07
            self.ci       = 1090000