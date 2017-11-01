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

class Input {

    private:
        int readNamelist();
        int checkItemExists(std::string, std::string, std::string el="default");
        
        template <class valuetype>
        int parseItem(valuetype*, std::string, std::string, std::string, bool, valuetype);
        int checkItem(int*, std::string, std::string, std::string el="default");
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
    
    public:
        Input();
        
        // Item retrieval functions
        int getItem(int*, std::string, std::string, std::string);
        int getItem(int*, std::string, std::string, std::string, int);
        int getItem(std::string*, std::string, std::string, std::string);
        int getItem(std::string*, std::string, std::string, std::string, std::string);
        
        // List retrieval function
        int getList(std::vector<std::string> *, std::string, std::string, std::string); 
};

#endif