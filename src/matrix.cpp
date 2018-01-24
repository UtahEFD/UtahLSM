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
                     std::vector<double> &u) {
        
        // local variables
        int j,n=a.size();
        double bet;
        std::vector<double> gam(n);
        
        // solve system
        if (b[0] == 0.0) throw(std::string("Error 1 in tridag"));
        u[1]=r[0]/(bet=b[0]);
        for (j=1;j<n;j++) {
            gam[j]=c[j-1]/bet;
            bet=b[j]-a[j]*gam[j];
            if (bet == 0.0) throw(std::string("Error 2 in tridag"));
            u[j+1]=(r[j]-a[j]*u[j])/bet;
        }
        for (j=(n-2);j>=0;j--)
            u[j+1] -= gam[j+1]*u[j+2];
    }
};