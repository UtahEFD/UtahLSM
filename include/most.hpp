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

#ifndef MOST_HPP
#define MOST_HPP

/**
 * This class computes functions relevant to
 * Monin-Obukhov similarity theory
 */

namespace most {
    
    // gradient functions
    double phim(const double);
    double phimStable(const double);
    double phimUnstable(const double);
    double phih(const double);
    double phihStable(const double);
    double phihUnstable(const double);
    
    // integral functions
    double psim(const double);
    double psimStable(const double);
    double psimUnstable(const double);
    double psih(const double);
    double psihStable(const double);
    double psihUnstable(const double);
    
    // log-law functions
    double fm(const double, const double, const double);
    double fh(const double, const double, const double);
};

#endif
