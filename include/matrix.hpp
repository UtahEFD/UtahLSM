//
//  matrix.hpp
//  
//  This namespace includes a simple tridiagonal solver
//
//  Created by Jeremy Gibbs on 12/1/17.
//

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
