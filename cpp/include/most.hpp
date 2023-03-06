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
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to gradient momentum function
     */
    double phim(const double zeta);
    
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
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to gradient scalar function
     */
    double phih(const double zeta);
    
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
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to integral stability correction for momentum
     */
    double psim(const double zeta);
    
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
     * @param[in] zeta ratio zeta = z/L
     * @return the solution to integral stability correction for scalars
     */
    double psih(const double zeta);
    
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
     * @param[in] zm_over_zo ratio zm_over_zo = zm / zo
     * @param[in] zeta_m ratio zeta_m = zm/L
     * @param[in] zeta_o ratio zeta_o = zo/L
     * @return the solution to log-law function for momentum
     */
    double fm(const double zm_over_zo, const double zeta_m, const double zeta_o);
    
    /**
     * Computes common log-law function for scalars. Here, 
     * zs is measurement height of scalars and zt is roughness 
     * height for scalars.
     *
     * @param[in] zs_over_zt ratio zs_over_zt = zs / zt
     * @param[in] zeta_s ratio zeta_s = zs/L
     * @param[in] zeta_t ratio zeta_t = zt/L
     * @return the solution to log-law function for scalars
     */
    double fh(const double zs_over_zt, const double zeta_s, const double zeta_t);
};

#endif
