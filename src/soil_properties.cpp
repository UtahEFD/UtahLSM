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

#include "soil_properties.hpp"

// Factory method that returns the correct soil properties
SoilType* SoilType::getProperties(int type, int source) {

    if (type==1) {
        return new Sand(source);
    } else if (type==2) {
        return new LoamySand(source);
    } else if (type==3) {
        return new SandyLoam(source);
    } else if (type==4) {
        return new SiltyLoam(source);
    } else if (type==5) {
        return new Loam(source);
    } else if (type==6) {
        return new SandyClayLoam(source);
    } else if (type==7) {
        return new SiltyClayLoam(source);
    } else if (type==8) {
        return new ClayLoam(source);
    } else if (type==9) {
        return new SandyClay(source);
    } else if (type==10) {
        return new SiltyClay(source);
    } else if (type==11) {
        return new Clay(source);
    } else if (type==12) {
        return new Peat(source);
    } else {
        std::cout<<"Invalid soil property dataset: must be an integer between 1-12"<<std::endl;
        throw(1);
    }
}