//
//  matrix.hpp
//  
//  This namespace includes a simple tridiagonal solver
//
//  Created by Jeremy Gibbs on 12/1/17.
//

#include <iostream>
#include <vector>

namespace matrix {
    
    // solve tridiagonal matrix
    void tridiagonal(const std::vector<double>& a,const std::vector<double>& b,
                     const std::vector<double>& c,const std::vector<double>& r,
                     std::vector<double> &u,int n) {
        
        // local variables
        double bet;
        double gam[n];
        
        // solve system
        bet  = b[0];
        u[1] = r[0] / bet;
        for (int j=1; j<n; j++) {
            gam[j] = c[j-1]/bet;
            bet    = b[j]-a[j]*gam[j];
            u[j+1]   = (r[j]-a[j]*u[j])/bet;
        }
        for (int j=n-1; j>=0; j--) {
            u[j+1]=u[j+1]-gam[j+1]*u[j+2];
        }        
    }
};