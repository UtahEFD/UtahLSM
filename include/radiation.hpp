//
//  radiation.hpp
//  
//  This namespace handles functions related to soil
//
//  Created by Jeremy Gibbs on 11/6/17.
//

#ifndef RADIATION_HPP
#define RADIATION_HPP

namespace radiation {
    
    double longwaveIn();
    double longwaveOut(const double,const double);
    double shortwaveIn(const double,const double,const double,const double);
    double shortwaveOut(const double,const double);
    double net(const double,const double,const double,
                        const double,const double,const double,const double);
};

#endif