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

#ifndef MOST_HPP
#define MOST_HPP

/**
 * Namespace for computing functions relevant to
 * Monin-Obukhov similarity theory.
 */
namespace most {
    
    /**
     * Computes gradient momentum function.
     *
     * @param[in] z height of the measurement
     * @param[in] obukL Obukhov length
     * @return the solution to gradient momentum function
     */
    double phim(const double z, const double obukL);
    
    /**
     * Helper function for phim when zeta >= 0 (neutral or stable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to phim under neutral or stable conditions
     */
    double phimStable(const double zeta);
    
    /**
     * Helper function for phim when zeta < 0 (unstable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to phim under unstable conditions
     */
    double phimUnstable(const double zeta);
    
    /**
     * Computes gradient scalar function.
     *
     * @param[in] z height of the measurement
     * @param[in] obukL Obukhov length
     * @return the solution to gradient scalar function
     */
    double phih(const double z, const double obukL);
    
    /**
     * Helper function for phih when zeta >= 0 (neutral or stable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to phih under neutral or stable conditions
     */
    double phihStable(const double zeta);
    
    /**
     * Helper function for phih when zeta < 0 (unstable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to phih under unstable conditions
     */
    double phihUnstable(const double zeta);
    
    /**
     * Computes the integral stability correction for momentum.
     *
     * @param[in] z height of the measurement
     * @param[in] obukL Obukhov length
     * @return the solution to integral stability correction for momentum
     */
    double psim(const double z, const double obukL);
    
    /**
     * Helper function for psim when zeta >= 0 (neutral or stable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to psim under neutral or stable conditions
     */
    double psimStable(const double zeta);
    
    /**
     * Helper function for psim when zeta < 0 (unstable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to psim under unstable conditions
     */
    double psimUnstable(const double zeta);
    
    /**
     * Computes the integral stability correction for scalars.
     *
     * @param[in] z height of the measurement
     * @param[in] obukL Obukhov length
     * @return the solution to integral stability correction for scalars
     */
    double psih(const double z, const double obukL);
    
    /**
     * Helper function for psih when zeta >= 0 (neutral or stable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to psih under neutral or stable conditions
     */
    double psihStable(const double zeta);
    
    /**
     * Helper function for psih when zeta < 0 (unstable) 
     *
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to psih under unstable conditions
     */
    double psihUnstable(const double zeta);
    
    /**
     * Computes common log-law function for momentum. Here, 
     * zm is measurement height of momentum and zo is roughness 
     * height for momentum.
     *
     * @param[in] z1 height of upper measurement
     * @param[in] z0 height of lower measurement
     * @param[in] obukL Obukhov length
     * @return the solution to log-law function for momentum
     */
    double fm(const double z1, const double z0, const double obukL);
    
    /**
     * Computes common log-law function for scalars. Here, 
     * zs is measurement height of scalars and zt is roughness 
     * height for scalars.
     *
     * @param[in] z1 height of upper measurement
     * @param[in] z0 height of lower measurement
     * @param[in] obukL Obukhov length
     * @return the solution to log-law function for scalars
     */
    double fh(const double z1, const double z0, const double obukL);
};

#endif
