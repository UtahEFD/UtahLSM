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

#include <iostream>
#include <iomanip>
#include <vector>
#include <span>

namespace matrix {
    
    // Solve tridiagonal matrix using the Thomas algorithm
    void tridiagonal(const std::vector<double>& a,const std::vector<double>& b,
                     const std::vector<double>& c,const std::vector<double>& r,
                     std::span<double>& u) {
        
        // Local variables
        int j;
        int n=int(a.size());
        double bet=b[0];
        std::vector<double> gam(n);
        std::cout<<std::endl;
        std::cout << std::fixed;
        std::cout<<"tridiag---------"<<std::endl;
        std::cout<<"n: "<<n<<std::endl;
           
        // Make sure diagonal band is not zero
        if (b[0] == 0.0) throw(std::string("Error 1 in tridag"));

        // Initialize first element of solution vector
        u[0] = r[0]/(bet);
        std::cout<<"----------------"<<std::endl;
        std::cout<<std::setprecision(17)<<"uj (0): "<<u[0]<<std::endl;
        std::cout<<std::setprecision(17)<<"rj (0): "<<r[0]<<std::endl;
        std::cout<<std::setprecision(17)<<"bj (0): "<<bet<<std::endl;
        
        // Forward sweep
        for (j=1;j<n;j++) {
            std::cout<<"----------------"<<std::endl;
            gam[j] = c[j-1]/bet;
            bet=b[j]-a[j]*gam[j];

            // Error check on bet
            if (bet == 0.0) throw(std::string("Error 2 in tridag"));

            u[j]=(r[j]-a[j]*u[j-1])/bet;
            std::cout<<"aj ("<<j<<"): "<<std::setprecision(17)<<a[j]<<std::endl;
            std::cout<<"bj ("<<j<<"): "<<std::setprecision(17)<<b[j]<<std::endl;
            std::cout<<"cj ("<<j<<"): "<<std::setprecision(17)<<c[j-1]<<std::endl;
            std::cout<<"bt ("<<j<<"): "<<std::setprecision(17)<<bet<<std::endl;
            std::cout<<"gm ("<<j<<"): "<<std::setprecision(17)<<gam[j]<<std::endl;
            std::cout<<"uj ("<<j<<"): "<<std::setprecision(17)<<u[j]<<std::endl;
        }
        std::cout<<"----------------"<<std::endl;

        // Backward sweep
        for (j=(n-2);j>=0;j--) {
            u[j] -= gam[j+1]*u[j+1];
            // std::cout<<"gj ("<<j<<"): "<<std::setprecision(17)<<gam[j+1]<<std::endl;
            // std::cout<<"uj ("<<j<<"): "<<std::setprecision(17)<<gam[j+1]<<std::endl;
        }
    }
};