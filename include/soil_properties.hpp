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

#ifndef SOILPROPERTIES_HPP
#define SOILPROPERTIES_HPP

#include <iostream>
#include <vector>

/**
 * Base struct for managing soil properties.
 * 
 * This strcut is responsible for holding the appropriate 
 * soil properties for a given dataset and soil type.
 * 
 * Datasets:
 * 1 = Clapp and Hornberger (1974)
 * 2 = Cosby et al. (1984)
 * 3 = Rawls and Brakensiek (1982)
 * 
 * Soil Types:
 * 01 = sand
 * 02 = loamy sand
 * 03 = sandy loam
 * 04 = silty loam
 * 05 = loam
 * 06 = sandy clay loam
 * 07 = silty clay loam
 * 08 = clay loam
 * 09 = sandy clay
 * 10 = silty clay
 * 11 = clay
 * 12 = peat
 * 
 */
struct SoilType {

    public:
        
        double b=0;        ///< exponent (unitless)
        double psi_sat=0;  ///< saturation moisture potential (m)
        double porosity=0; ///< saturated soil moisture 
        double residual=0; ///< residual moisture (volume/volume)
        double K_sat=0;    ///< hydraulic conductivity (m/s)
        double ci=0;       ///< volumetric heat capacity (J/m^3/K)

        /**
         * Factory method that returns the correct soil properties.
         * 
         * @param[in] type the soil type
         * @param[in] dataset soil properties dataset
         */
        static SoilType* getProperties(int type, int dataset);
};

/**
 * Derived struct for managing properties for sand (type=1).
 */
struct Sand : public SoilType {
    Sand(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 4.05;
            psi_sat  = -0.121;
            porosity = 0.395;
            residual = 0.0000;
            K_sat    = 1.76e-04;
            ci       = 1470000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 2.79;
            psi_sat  = -0.023;
            porosity = 0.339;
            residual = 0.0000;
            K_sat    = 1.60e-05;
            ci       = 1470000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 1.44;
            psi_sat  = -0.160;
            porosity = 0.437;
            residual = 0.0200;
            K_sat    = 5.83e-05;
            ci       = 1470000;
        }
    }   
};

/**
 * Derived struct for managing properties for loamy sand (type=2).
 */
struct LoamySand : public SoilType {
    LoamySand(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 4.38;
            psi_sat  = -0.090;
            porosity = 0.410;
            residual = 0.0000;
            K_sat    = 1.56e-04;
            ci       = 1410000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 4.26;
            psi_sat  = -0.018;
            porosity = 0.421;
            residual = 0.0000;
            K_sat    = 9.52e-06;
            ci       = 1410000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 1.00;
            psi_sat  = -0.206;
            porosity = 0.437;
            residual = 0.0350;
            K_sat    = 1.70e-05;
            ci       = 1410000;
        }
    }   
};

/**
 * Derived struct for managing properties for sandy loam (type=3).
 */
struct SandyLoam : public SoilType {
    SandyLoam(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 4.90;
            psi_sat  = -0.218;
            porosity = 0.435;
            residual = 0.0000;
            K_sat    = 3.41e-05;
            ci       = 1340000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 4.74;
            psi_sat  = -0.032;
            porosity = 0.434;
            residual = 0.0000;
            K_sat    = 6.19e-06;
            ci       = 1340000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 81.00;
            psi_sat  = -0.302;
            porosity = 0.453;
            residual = 0.0410;
            K_sat    = 7.19e-06;
            ci       = 1340000;
        }
    }   
};

/**
 * Derived struct for managing properties for silty loam (type=4).
 */
struct SiltyLoam : public SoilType {
    SiltyLoam(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 5.30;
            psi_sat  = -0.786;
            porosity = 0.485;
            residual = 0.0000;
            K_sat    = 7.20e-06;
            ci       = 1270000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 5.33;
            psi_sat  = -0.066;
            porosity = 0.476;
            residual = 0.0000;
            K_sat    = 4.73e-06;
            ci       = 1270000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 2.65;
            psi_sat  = -0.401;
            porosity = 0.463;
            residual = 0.0270;
            K_sat    = 1.89e-06;
            ci       = 1270000;
        }
    }   
};

/**
 * Derived struct for managing properties for loam (type=5).
 */
struct Loam : public SoilType {
    Loam(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 5.39;
            psi_sat  = -0.478;
            porosity = 0.451;
            residual = 0.0000;
            K_sat    = 7.00e-06;
            ci       = 1210000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 5.25;
            psi_sat  = -0.047;
            porosity = 0.439;
            residual = 0.0000;
            K_sat    = 5.12e-06;
            ci       = 1210000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 3.97;
            psi_sat  = -0.509;
            porosity = 0.501;
            residual = 0.0150;
            K_sat    = 3.67e-06;
            ci       = 1210000;
        }
    }   
};

/**
 * Derived struct for managing properties for sandy clay loam (type=6).
 */
struct SandyClayLoam : public SoilType {
    SandyClayLoam(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 7.12;
            psi_sat  = -0.299;
            porosity = 0.420;
            residual = 0.0000;
            K_sat    = 6.30e-06;
            ci       = 1180000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 6.77;
            psi_sat  = -0.031;
            porosity = 0.404;
            residual = 0.0000;
            K_sat    = 5.78e-06;
            ci       = 1180000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 4.27;
            psi_sat  = -0.594;
            porosity = 0.398;
            residual = 0.0680;
            K_sat    = 1.19e-06;
            ci       = 1180000;
        }
    }   
};

/**
 * Derived struct for managing properties for silty clay loam (type=7).
 */
struct SiltyClayLoam : public SoilType {
    SiltyClayLoam(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 7.75;
            psi_sat  = -0.356;
            porosity = 0.477;
            residual = 0.0000;
            K_sat    = 1.70e-06;
            ci       = 1320000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 8.72;
            psi_sat  = -0.060;
            porosity = 0.464;
            residual = 0.0000;
            K_sat    = 4.11e-06;
            ci       = 1320000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 3.13;
            psi_sat  = -0.564;
            porosity = 0.464;
            residual = 0.0750;
            K_sat    = 6.39e-07;
            ci       = 1320000;
        }
    }   
};

/**
 * Derived struct for managing properties for clay loam (type=8).
 */
struct ClayLoam : public SoilType {
    ClayLoam(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 8.52;
            psi_sat  = -0.630;
            porosity = 0.476;
            residual = 0.0000;
            K_sat    = 2.50e-06;
            ci       = 1230000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 8.17;
            psi_sat  = -0.041;
            porosity = 0.465;
            residual = 0.0000;
            K_sat    = 4.45e-06;
            ci       = 1230000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 4.13;
            psi_sat  = -0.703;
            porosity = 0.471;
            residual = 0.0400;
            K_sat    = 4.17e-07;
            ci       = 1230000;
        }
    }   
};

/**
 * Derived struct for managing properties for sandy clay (type=9).
 */
struct SandyClay : public SoilType {
    SandyClay(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 10.40;
            psi_sat  = -0.153;
            porosity = 0.426;
            residual = 0.0000;
            K_sat    = 2.20e-06;
            ci       = 1180000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 10.73;
            psi_sat  = -0.027;
            porosity = 0.406;
            residual = 0.0000;
            K_sat    = 7.12e-05;
            ci       = 1180000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 5.65;
            psi_sat  = -0.795;
            porosity = 0.430;
            residual = 0.1090;
            K_sat    = 3.33e-07;
            ci       = 1180000;
        }
    }   
};

/**
 * Derived struct for managing properties for silty clay (type=10).
 */
struct SiltyClay : public SoilType {
    SiltyClay(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 10.40;
            psi_sat  = -0.490;
            porosity = 0.492;
            residual = 0.0000;
            K_sat    = 1.00e-06;
            ci       = 1150000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 10.39;
            psi_sat  = -0.045;
            porosity = 0.468;
            residual = 0.0000;
            K_sat    = 3.43e-06;
            ci       = 1150000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 4.48;
            psi_sat  = -0.765;
            porosity = 0.479;
            residual = 0.0560;
            K_sat    = 2.50e-07;
            ci       = 1150000;
        }
    }   
};

/**
 * Derived struct for managing properties for clay (type=11).
 */
struct Clay : public SoilType {
    Clay(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 11.40;
            psi_sat  = -0.405;
            porosity = 0.482;
            residual = 0.0000;
            K_sat    = 1.30e-06;
            ci       = 1090000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 10.55;
            psi_sat  = -0.053;
            porosity = 0.468;
            residual = 0.0000;
            K_sat    = 2.99e-06;
            ci       = 1090000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 6.67;
            psi_sat  = -0.856;
            porosity = 0.475;
            residual = 0.0900;
            K_sat    = 1.67e-07;
            ci       = 1090000;
        }
    }   
};

/**
 * Derived struct for managing properties for peat (type=12).
 */
struct Peat : public SoilType {
    Peat(int dataset) {
        // Clapp and Hornberger (1974)
        if (dataset==1) {
            b        = 7.75;
            psi_sat  = -0.356;
            porosity = 0.863;
            residual = 0.0000;
            K_sat    = 8.00e-06;
            ci       = 840000;
        }
        // Cosby et al. (1984)
        if (dataset==2) {
            b        = 7.75;
            psi_sat  = -0.356;
            porosity = 0.863;
            residual = 0.0000;
            K_sat    = 8.00e-06;
            ci       = 840000;
        }
        // Rawls and Brakensiek (1982)
        if (dataset==3) {
            b        = 6.06;
            psi_sat  = -0.356;
            porosity = 0.863;
            residual = 0.1763;
            K_sat    = 8.00e-06;
            ci       = 840000;
        }
    }   
};

/**
 * Derived struct for managing properties for B11 (type=13).
 */
struct B11 : public SoilType {
    B11(int dataset) {
        // Heinen, Bakker, Wosten (Cabauw-specific)
        if (dataset==4) {
            b        = 9.35;
            psi_sat  = -0.463;
            porosity = 0.591;
            residual = 0.01;
            K_sat    = 7.30e-07;
            ci       = 1090000;
        } else {
            std::cout<<"Invalid soil property dataset: this is specific to Cabauw and must equal 4."<<std::endl;
            throw(1);
        }
    }   
};

/**
 * Derived struct for managing properties for O12 (type=14).
 */
struct O12 : public SoilType {
    O12(int dataset) {
        // Heinen, Bakker, Wosten (Cabauw-specific)
        if (dataset==4) {
            b        = 6.33;
            psi_sat  = -1.14;
            porosity = 0.561;
            residual = 0.01;
            K_sat    = 1.25e-07;
            ci       = 1090000;
        } else {
            std::cout<<"Invalid soil property dataset: this is specific to Cabauw and must equal 4."<<std::endl;
            throw(1);
        }
    }   
};


/**
 * Derived struct for managing properties for O16 (type=15).
 */
struct O16 : public SoilType {
    O16(int dataset) {
        // Heinen, Bakker, Wosten (Cabauw-specific)
        if (dataset==4) {
            b        = 2.75;
            psi_sat  = -1.03;
            porosity = 0.889;
            residual = 0.00;
            K_sat    = 1.69e-07;
            ci       = 1090000;
        } else {
            std::cout<<"Invalid soil property dataset: this is specific to Cabauw and must equal 4."<<std::endl;
            throw(1);
        }
    }   
};

#endif