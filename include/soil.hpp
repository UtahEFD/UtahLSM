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

#include <vector>
#include <stdexcept>

/**
 * This class handles tasks associated with soil.
 */

namespace soil {
	
	struct MoistureTransfer {
        std::vector<double> d; 
        std::vector<double> k;
    };
    
    struct ThermalTransfer {
        std::vector<double> d;
        std::vector<double> k;
    };
    
    struct Properties {
	    std::vector<double> b;
	    std::vector<double> psi_sat;
	    std::vector<double> porosity;
        std::vector<double> residual;
	    std::vector<double> K_sat;
	    std::vector<double> Ci;
    };
    
    double heatCapacity(const double, const double, const double);
    
    double surfaceMixingRatio(const double,const double, const double,
                              const double,const double, const double,
                              const double, const int);
    
    double surfaceWaterContent(const double, const double, const double,
                               const double, const double, const int);
    
    double surfaceWaterContentEstimate(const double,const double,const double,
                                       const double,const double,const double,
                                       const double,const int);
    
    double heatCapacity(const double, const double, const double);
    
    double waterPotential(const double, const double, const double,
                          const double, const double, const int);
    
    std::vector<double> waterPotential(const std::vector<double>&, const std::vector<double>&,
                                       const std::vector<double>&, const std::vector<double>&,
                                       const std::vector<double>&, const int, const int);
    
    ThermalTransfer thermalTransfer(const std::vector<double>&, const std::vector<double>&,
                                    const std::vector<double>&, const std::vector<double>&,
                                    const std::vector<double>&, const std::vector<double>&,
                                    const int, const int);
    
    MoistureTransfer moistureTransfer(const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const std::vector<double>&,
                                      const int, const int);
    
    Properties properties(const std::vector<int>&, const int, const int);
};

#endif
