/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2023 Jeremy A. Gibbs
 * Copyright (c) 2017–2023 Rob Stoll
 * Copyright (c) 2017–2023 Eric Pardyjak
 * Copyright (c) 2017–2023 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "sfc.hpp"
#include "sfc_most.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

Surface::Surface(const int sfc_model) {};

// Factory method that returns the correct surface model
Surface* Surface::getModel(const int sfc_model) {

    if (sfc_model==1) {
        std::cout<<"[UtahLSM: Surface] \t --- using the MOST model"<<std::endl;
        return new MOST(sfc_model);
    } else {
        std::cout<<"[UtahLSM: Surface] \t Invalid surface model: must be =1"<<std::endl;
        throw(1);
    }
}
