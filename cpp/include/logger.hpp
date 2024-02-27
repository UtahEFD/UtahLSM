/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2024 Jeremy A. Gibbs
 * Copyright (c) 2017–2024 Rob Stoll
 * Copyright (c) 2017–2024 Eric Pardyjak
 * Copyright (c) 2017–2024 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <vector>

class Logger {
    public:
        void print_double(double x, std::string section, std::string label);
        void print_double(double x);
        void print_hex(double x, std::string section, std::string label);
        void print_hex(double x);
};

#endif