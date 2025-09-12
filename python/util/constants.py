# 
# UtahLSM
# 
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

vonk    = 0.4              # von Karman constant []
grav    = 9.81             # acceleration due to gravity [m/s^2]
pi      = 3.14159265358979 # pi
sb      = 5.6697e-8        # Stefan-Boltzmann constant [W/m^2-K^4]
sc      = 1367             # solar constant [K-m/s]
rho_air = 1.204            # density of air [kg/m^3]
rho_wat = 1000.0           # density of water [kg/m^3]
Rd      = 287              # gas constant for dry air [J/kg-K]
Rv      = 461.4            # gas constant for water vapor [J/kg-K]
epsilon = Rd/Rv            # ratio of mass of dry air to moist air
Lv      = 2.45e6           # latent heat of vaporization [J/kg]
Ci_wat  = 4186000.0        # volumetric heat capacity of water [J/m^3-K]
Cp_air  = 1004.0           # specific heat of air [J/kg-K]