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

#ifndef RADIATION_HPP
#define RADIATION_HPP

/**
 * Class for managing radiation.
 * 
 * This class computes the surface radiation budget.
 */
class Radiation {
    
    public:
        
        /**
         * Constructs a Radiation object.
         *
         * @param[in] latitude latitude of the site
         * @param[in] longitude longitude of the site
         * @param[in] albedo surface albedo of the site
         * @param[in] emissivity surface emissivity of the site
         */
        Radiation(const double latitude, const double longitude,
                  const double albedo, const double emissivity);
        
        /**
         * Computes the surface net radiation.
         *
         * @param[in] julian_day julian day of the year
         * @param[in] time_utc time of the day in UTC
         * @param[in] sfc_T surface temperature
         */
        double computeNet(const double julian_day, const double time_utc, const double sfc_T);

    private:
        
        // Local copies
        double latitude;    ///< site latitude
        double longitude;   ///< site longitude
        double albedo;      ///< surface albedo
        double emissivity;  ///< surface emissivity
        
        /**
         * Computes the downward longwave radiation at the surface.
         *
         * @return downward longwave radiation 
         */
        double longwaveIn();

        /**
         * Computes the upward longwave radiation at the surface.
         *
         * @param[in] emissivity surface emissivity
         * @param[in] sfc_T surface temperature
         * @return    upward longwave radiation
         */
        double longwaveOut(const double emissivity, const double sfc_T);

        /**
         * Computes the downward shortwave radiation at the surface.
         *
         * @param[in] julian_day julian day of the year
         * @param[in] time_utc time of the day in UTC
         * @param[in] latitude latitude of the site
         * @param[in] longitude of the site
         * @return    downward shortwave radiation
         */
        double shortwaveIn(const double julian_day,const double time_utc,
                           const double latitude,const double longitude);
        
        /**
         * Computes the upward shortwave radiation at the surface.
         *
         * @param[in] albedo surface albedo of the site
         * @param[in] shortwave_in downward shortwave radiation at the surface
         * @return    upward shortwave radiation
         */
        double shortwaveOut(const double albedo, const double shortwave_in);  
};

#endif