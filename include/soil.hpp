/*
 * UtahLSM
 * 
 * Copyright (c) 2019 Jeremy A. Gibbs
 * Copyright (c) 2019 Pete Willemsen
 * Copyright (c) 2019 Rob Stoll
 * Copyright (c) 2019 Eric Pardyjak
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

/**
 * Namespace for managing soil properties
 */
class Soil {
	
    public:
        /**
         * Struct to hold soil moisture transfer coefficients.
         */
	    struct MoistureTransfer {
            std::vector<double> d; ///< soil moisture diffusivity
            std::vector<double> k; ///< soil moisture conductivity
        };

        /**
         * Struct to hold soil thermal transfer coefficients.
         */
        struct ThermalTransfer {
            std::vector<double> d; ///< soil moisture diffusivity
            std::vector<double> k; ///< soil moisture conductivity
        };

        /**
         * Struct to hold soil properties.
         */
        struct Properties {
	        std::vector<double> b;        ///< b exponent 
	        std::vector<double> psi_sat;  ///< saturated soil moisture potential
	        std::vector<double> porosity; ///< saturated soil moisture 
            std::vector<double> residual; ///< minimum soil moisture
	        std::vector<double> K_sat;    ///< saturated hydraulic conductivity
	        std::vector<double> Ci;       ///< dry volumetric heat capacity
        };

        /**
         * Computes heat capacity of soil.
         *
         * @param[in] porosity saturated soil moisture
         * @param[in] Ci dry volumetric heat capacity
         * @param[in] soil_q soil moisture
         * @return heat capacity of soil
         */
        double heatCapacity(const double porosity, const double Ci, const double soil_q);

        /**
         * Computes surface mixing ratio.
         *
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] b b exponent 
         * @param[in] sfc_T surface moisture
         * @param[in] sfc_q surface moisture
         * @param[in] atm_p atmospheric pressure
         * @param[in] model soil model
         * @return surface mixing ratio
         */
        double surfaceMixingRatio(const double psi_sat, const double porosity,
                                  const double residual,  const double b,
                                  const double sfc_T, const double sfc_q,
                                  const double atm_p, const int model);

        /**
         * Computes soil surface moisture.
         *
         * @param[in] psi soil moisture potential
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] b b exponent
         * @param[in] model soil model
         * @return soil surface moisture
         */
        double surfaceWaterContent(const double psi, const double psi_sat,
                                   const double porosity, const double residual,
                                   const double b, const int model);

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
        double surfaceWaterContentEstimate(const double psi_sat, const double porosity,
                                           const double residual, const double b,
                                           const double sfc_T, const double sfc_q,
                                           const double atm_p, const int model);

        /**
         * Computes soil water potential (single level).
         *
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] soil_q soil moisture
         * @param[in] b b exponent 
         * @param[in] model soil model
         * @return soil water potential
         */
        double waterPotential(const double psi_sat, const double porosity,
                              const double residual, const double soil_q,
                              const double b, const int model);

        /**
         * Computes soil water potential (full column).
         *
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] soil_q soil moisture
         * @param[in] b b exponent 
         * @param[in] depth number of soil levels
         * @param[in] model soil model
         * @return soil water potential
         */
        std::vector<double> waterPotential(const std::vector<double>& psi_sat,
                                           const std::vector<double>& porosity,
                                           const std::vector<double>& residual,
                                           const std::vector<double>& soil_q,
                                           const std::vector<double>& b,
                                           const int depth, const int model);

        /**
         * Computes soil thermal conductivity/diffusivity.
         *
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] soil_q soil moisture
         * @param[in] b b exponent
         * @param[in] Ci dry volumetric heat capacity
         * @param[in] depth number of soil levels
         * @param[in] model soil model
         * @return soil thermal conductivity/diffusivity
         */
        ThermalTransfer thermalTransfer(const std::vector<double> &psi_sat,
                                        const std::vector<double> &porosity,
                                        const std::vector<double> &residual,
                                        const std::vector<double> &soil_q,
                                        const std::vector<double> &b,
                                        const std::vector<double> &Ci,
                                        const int depth, const int model);

        /**
         * Computes soil moisture conductivity/diffusivity.
         *
         * @param[in] psi_sat saturated soil moisture potential
         * @param[in] K_sat saturated hydraulic conductivity
         * @param[in] porosity saturated soil moisture
         * @param[in] residual minimum soil moisture
         * @param[in] soil_q soil moisture
         * @param[in] b b exponent
         * @param[in] depth number of soil levels
         * @param[in] model soil model
         * @return soil moisture conductivity/diffusivity
         */
        MoistureTransfer moistureTransfer(const std::vector<double>& psi_sat,
                                          const std::vector<double>& K_sat,
                                          const std::vector<double>& porosity,
                                          const std::vector<double>& residual,
                                          const std::vector<double>& soil_q,
                                          const std::vector<double>& b,
                                          const int depth, const int model);

        /**
         * Set soil type properties at each soil level.
         *
         * @param[in] soil_type soil type at each soil level
         * @param[in] depth number of soil levels
         * @param[in] src soil parameter set
         * @return soil type properties
         */
        Properties properties(const std::vector<int>& soil_type, const int depth, const int src);
};

#endif