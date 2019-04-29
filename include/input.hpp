//
//  input.hpp
//  
//  This class is an input manager
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include "json.hpp"

using json = nlohmann::json;

/**
 * Class for managing input files.
 * This class is responsible for opening a supplied input file
 * and returning requested fields from that file.
 */

class Input {
    
    private:

        json input;                      ///< json wrapper of input file
        void readInputFile(std::string); ///< Read the input file
    
    public:
    
        Input(std::string);        
    
        // getters
        void getItem(int&, std::string, std::string);
        void getItem(double&, std::string, std::string);
        void getItem(std::vector<int>&, std::string, std::string);
        void getItem(std::vector<double>&, std::string, std::string);
        void getItem(std::vector<std::string>&, std::string, std::string);
};

// C-style functions
typedef void * InputObject;

extern "C" {
   InputObject GetInput(char*);
   void GetItemInt(InputObject,int*,char*,char*);
   void GetItemDbl(InputObject,double*,char*,char*);
   void GetItemDblArr(InputObject,double[],int*,char*,char*);
}

#endif
