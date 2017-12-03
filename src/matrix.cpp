//
//  matrix.hpp
//  
//  This namespace includes a simple tridiagonal solver
//
//  Created by Jeremy Gibbs on 12/1/17.
//

#include <iostream>

namespace matrix {
    
    // solve tridiagonal matrix
    void tridiagonal(double* a,double* b,double* c,double* r,double* u,int n) {
        
        // local variables
        double bet;
        double gam[n];
        
        // solve system
        bet  = b[0];
        u[0] = r[0] / bet;
        for (int j=1; j<n; j++) {
            gam[j] = c[j-1]/bet;
            bet    = b[j]-a[j]*gam[j];
            u[j]   = (r[j]-a[j]*u[j-1])/bet;
        }
        for (int j=n-1; j>=0; j--) {
            u[j]=u[j]-gam[j+1]*u[j+1];
        }        
    }
};