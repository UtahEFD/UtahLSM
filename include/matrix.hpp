/*
 * UtahLSM
 * 
 * Copyright (c) 2021 Jeremy A. Gibbs
 * Copyright (c) 2021 Rob Stoll
 * Copyright (c) 2021 Eric Pardyjak
 * Copyright (c) 2021 Pete Willemsen
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
 * Class for managing matrix algebra.
 * 
 * This class currently solves a tridiagonal matrix using
 * the Thomas algorithm.
 */
namespace matrix {

    /**
     * Solves a tridiagonal matrix using the Thomas algorithm.
     *
     * @param[in]  a first off-diagonal band
     * @param[in]  b diagonal band
     * @param[in]  c second off-diagonal band
     * @param[in]  r right-hand column vector
     * @param[out] u solution vector
     */
    void tridiagonal(const std::vector<double>& a,const std::vector<double>& b,
                     const std::vector<double>& c,const std::vector<double>& r,
                     std::vector<double> &u);
};

#endif
