/*
 * UtahLSM
 * 
 * Copyright (c) 2021 Jeremy A. Gibbs
 * Copyright (c) 2021 Rob Stoll
 * Copyright (c) 2021 Eric Pardyjak
 * Copyright (c) 2021 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#ifndef BROOKSCOREY_HPP
#define BROOKSCOREY_HPP

#include "soil.hpp"

#include <stdexcept>
#include <vector>

#include "soil_properties.hpp"

/**
 * Class for managing the soil model
 * 
 * This class computes various terms using a soil model.
 */
class BrooksCorey : public Soil {
	
    public:

        /**
         * Constructs a Soil object.
         *
         * @param[in] soil_type soil type at each depth
         * @param[in] soil_param soil properties dataset
         * @param[in] soil_model soil model
         * @param[in] depth number of soil levels
         */
        BrooksCorey(const std::vector<int>& soil_type, const int soil_param, 
                    const int soil_model, const int levels) : 
                    Soil(soil_type,soil_param,soil_model,levels){}

        /**
         * Computes soil surface moisture.
         *
         * @param[in] psi soil moisture potential
         * @return soil surface moisture
         */
        virtual double surfaceWaterContent(const double psi);

        /**
         * Estimate soil surface moisture from surface mixing ratio.
         *
         * @param[in] sfc_T surface moisture
         * @param[in] sfc_q surface moisture
         * @param[in] atm_p atmospheric pressure
         * @return soil surface moisture
         */
        virtual double surfaceWaterContentEstimate(const double sfc_T, const double sfc_q,
                                                   const double atm_p);

        /**
         * Computes soil water potential (single level).
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil water potential
         */
        virtual double waterPotential(const double soil_q, const int level);

        /**
         * Computes soil moisture conductivity.
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil moisture conductivity/diffusivity
         */
        virtual double conductivityMoisture(const double soil_q, const int level);

        /**
         * Computes soil moisture diffusivity.
         *
         * @param[in] soil_q soil moisture
         * @param[in] level soil level for computation
         * @return soil moisture conductivity/diffusivity
         */
        virtual double diffusivityMoisture(const double soil_q, const int level);
};

#endif