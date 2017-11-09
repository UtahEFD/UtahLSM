//
//  most.hpp
//  
//  This namespace handles Monin-Obukhov Similarity Theory
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#ifndef MOST_HPP
#define MOST_HPP

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