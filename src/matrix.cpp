/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2022 Jeremy A. Gibbs
 * Copyright (c) 2017–2022 Rob Stoll
 * Copyright (c) 2017–2022 Eric Pardyjak
 * Copyright (c) 2017–2022 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include <iostream>
#include <vector>
#include <span>

namespace matrix {
    
    // Solve tridiagonal matrix using the Thomas algorithm
    void tridiagonal(const std::vector<double>& a,const std::vector<double>& b,
                     const std::vector<double>& c,const std::vector<double>& r,
                     std::span<double> u) {
        
        // Local variables
        int j,n=int(a.size());
        double bet;
        std::vector<double> gam(n);
        
        // Make sure diagonal band is not zero
        if (b[0] == 0.0) throw(std::string("Error 1 in tridag"));

        // Initialize first element of solution vector
        u[0] = r[0]/(bet=b[0]);
        
        // Forward sweep 
        for (j=1;j<n;j++) {
            
            gam[j] = c[j-1]/bet;
            bet=b[j]-a[j]*gam[j];

            // Error check on bet
            if (bet == 0.0) throw(std::string("Error 2 in tridag"));

            u[j]=(r[j]-a[j]*u[j-1])/bet;
        }

        // Backward sweep
        for (j=(n-2);j>=0;j--)
            u[j] -= gam[j+1]*u[j+1];
    }
};