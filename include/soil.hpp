/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2022 Jeremy A. Gibbs
 * Copyright (c) 2017–2022 Rob Stoll
 * Copyright (c) 2017–2022 Eric Pardyjak
 * Copyright (c) 2017–2022 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#ifndef SOIL_HPP
#define SOIL_HPP

#include <stdexcept>
#include <vector>

#include "soil_properties.hpp"

/**
 * Class for managing the soil model
 * 
 * This class computes various terms using a soil model.
 */
class Soil {
	
    public:

        std::vector<SoilType*> properties; ///< vector of soil properties
    
        /**
         * Constructs a Soil object.
         *
         * @param[in] soil_type soil type at each depth
         * @param[in] soil_param soil properties dataset
         * @param[in] soil_model soil model
         * @param[in] depth number of soil levels
         */
        Soil(const std::vector<int>& soil_type, const int soil_param, 
             const int soil_model, const int levels);

        /**
         * Computes heat capacity of soil.
         *
         * @param[in] porosity saturated soil moisture
         * @param[in] Ci dry volumetric heat capacity
         * @param[in] soil_q soil moisture
         * @return heat capacity of soil
         */
        double heatCapacity(const double soil_q, const int level);

        /**
         * Computes surface mixing ratio.
         *
         * @param[in] sfc_T surface moisture
         * @param[in] sfc_q surface moisture
         * @param[in] atm_p atmospheric pressure
         * @return surface mixing ratio
         */
        double surfaceMixingRatio(const double sfc_T, const double sfc_q, const double atm_p);

        /**
         * Computes soil surface moisture.
         *
         * @param[in] psi soil moisture potential
         * @return soil surface moisture
         */
        virtual double surfaceWaterContent(const double psi)=0;

        /**
         * Estimate soil surface moisture from surface mixing ratio.
         *
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] b b exponent 
         * @param[in] sfc_T surface moisture
         * @param[in] sfc_q surface moisture
         * @param[in] atm_p atmospheric pressure
         * @param[in] model soil model
         * @return soil surface moisture
         */
        virtual double surfaceWaterContentEstimate(const double sfc_T, const double sfc_q,
                                                   const double atm_p)=0;

        /**
         * Computes soil water potential (single level).
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil water potential
         */
        virtual double waterPotential(const double soil_q, const int level)=0;
        
        /**
         * Computes soil thermal conductivity.
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil thermal conductivity
         */
        double conductivityThermal(const double soil_q, const int level);

        /**
         * Computes soil thermal diffusivity.
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil thermal diffusivity
         */
        double diffusivityThermal(const double soil_q, const int level);

        /**
         * Computes soil moisture conductivity.
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil moisture conductivity/diffusivity
         */
        virtual double conductivityMoisture(const double soil_q, const int level)=0;

        /**
         * Computes soil moisture diffusivity.
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil moisture conductivity/diffusivity
         */
        virtual double diffusivityMoisture(const double soil_q, const int level)=0;

        /**
         * Factory method that returns the correct soil model.
         * 
         * @param[in] model soil model
         * @return pointer to proper soil model
         */
        static Soil* getModel(const std::vector<int>& soil_type, const int soil_param, 
                              const int soil_model, const int levels);
};
#endif