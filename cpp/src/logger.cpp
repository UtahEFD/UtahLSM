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

#include "logger.hpp"

#include <iostream>
#include <iomanip>

void Logger :: print_double(double x, std::string label) {
    std::cout<<std::defaultfloat;
    std::cout<<std::setprecision(17);
    std::cout<<"LOG -> "<<label<<" "<<x<<std::endl;
}

void Logger :: print_double(double x) {
    std::cout<<std::defaultfloat;
    std::cout<<std::setprecision(17);
	std::cout<<"LOG -> "<<x<<std::endl;
}

void Logger :: print_hex(double x, std::string label) {
    std::cout<<std::hexfloat;
    std::cout<<"LOG -> "<<label<<" "<<x<<std::endl;
}

void Logger :: print_hex(double x) {
    std::cout<<std::hexfloat;
    std::cout<<"LOG -> "<<x<<std::endl;
}