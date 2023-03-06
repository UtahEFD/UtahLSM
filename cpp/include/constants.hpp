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

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/**
 * Namespace for storing commonly used constants.
 */
namespace constants {
    const double vonk    = 0.4;              ///< von Karman constant
    const double grav    = 9.81;             ///< acceleration due to gravity [m/s^2]
    const double pi      = 3.14159265358979; ///< pi
    const double sb      = 5.6697e-8;        ///< Stefan-Boltzmann constant [W/m^2-K^4]
    const double sc      = 1367;             ///< solar constant [K-m/s]
    const double rho_air = 1.204;            ///< density of air [kg/m^3]
    const double rho_wat = 1000.0;           ///< density of water [kg/m^3]
    const double Rv      = 461.4;            ///< gas constant for water vapor [J/kg-K]
    const double Lv      = 2.45e6;           ///< latent heat of vaporization [J/kg]
    const double Ci_wat  = 4186000.0;        ///< volumetric heat capacity of water [J/m^3-K]
    const double Cp_air  = 1004.0;           ///< specific heat of air [J/kg-K]
};

#endif