//
//  input.hpp
//  
//  This class handles reading in user options
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#ifndef INPUT_HPP
#define INPUT_HPP

#include <map>
#include <string>
#include <vector>

typedef std::map<std::string, std::vector<double> > dataMap;

/**
 * This class handles reading input files.
 */

class Input {

    private:
        int readNamelist();
        int readDataFile(dataMap*, std::string, bool);
        int checkItemExists(std::string, std::string, std::string el="default");
        
        template <class valuetype>
        int parseItem(valuetype*,   std::string, std::string, std::string, bool, valuetype);
        int checkItem(int*,         std::string, std::string, std::string el="default");
        int checkItem(double*     , std::string, std::string, std::string el="default");
        int checkItem(std::string*, std::string, std::string, std::string el="default");
        
        template <class valuetype>
        int parseList(std::vector<valuetype> *, std::string, std::string, std::string);
        int checkList(std::vector<std::string>*, std::string, std::string, std::string el="default");
        
        struct inputType {
            std::string data;
            bool isused;
        };
        
        typedef std::map<std::string, inputType > inputMap1d;
        typedef std::map<std::string, inputMap1d> inputMap2d;
        typedef std::map<std::string, inputMap2d> inputMap;
    
        inputMap inputList;
        dataMap inputMetr;
        dataMap inputSoil;
    
    public:
        Input();
        
        // Item retrieval functions
        int getItem(int*,         std::string, std::string, std::string);
        int getItem(int*,         std::string, std::string, std::string, int);
        int getItem(double*,      std::string, std::string, std::string);
        int getItem(double*,      std::string, std::string, std::string, double);
        int getItem(std::string*, std::string, std::string, std::string);
        int getItem(std::string*, std::string, std::string, std::string, std::string);
                
        // List retrieval function
        int getList(std::vector<std::string> *, std::string, std::string, std::string);
        int getProf(double*, std::string, std::string, int size);
        int getProf(std::vector<int>*, std::string, std::string, int size);
        int getProf(std::vector<double>*, std::string, std::string, int size);        
};

#endif
