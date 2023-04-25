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

#ifndef SFC_HPP
#define SFC_HPP

#include <stdexcept>
#include <vector>

/**
 * Class for managing the surface model
 * 
 * This class computes various terms using a surface model.
 */
class Surface {
	
    public:
    
        /**
         * Constructs a Surface object.
         *
         * @param[in] sfc_model surface model
         */
        Surface(const int sfc_model);
        
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
        virtual double fm(const double z1, const double z0, const double obukL)=0;
        
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
        virtual double fh(const double z1, const double z0, const double obukL)=0;

        /**
         * Factory method that returns the correct surface model.
         * 
         * @param[in] model surface model
         * @return pointer to proper surface model
         */
        static Surface* getModel(const int sfc_model);
};
#endif