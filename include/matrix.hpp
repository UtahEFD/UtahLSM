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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

/**
 * This class solves a tridiagonal matrix using
 * the Thomas algorithm
 */

namespace matrix {

    void tridiagonal(const std::vector<double>&,const std::vector<double>&,
                     const std::vector<double>&,const std::vector<double>&,
                     std::vector<double>&);
};

#endif
