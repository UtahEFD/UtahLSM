//
//  input.cpp
//  
//  This class handles reading in user options
//  This is modified from version in MicroHH
//
//  Created by Jeremy Gibbs on 6/29/17.
//

#include "input.hpp"
#include "utah_lsm.hpp"
#include <string>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include<vector>

Input :: Input() {
    
    int nError          = 0;
    const bool required = false;
    
    std::cout<<"Reading namelist.ini"<<std::endl;
    nError += readNamelist();
    std::cout<<"Reading inputMetr.dat"<<std::endl;
    nError += readDataFile(&inputMetr, + "inputMetr.dat", required);
    std::cout<<"Reading inputSoil.dat"<<std::endl;
    nError += readDataFile(&inputSoil, + "inputSoil.dat", required);
    std::cout<<"##############################################################"<<std::endl;
    
    if (nError) throw 1;
}

// Read in the namelist file
int Input :: readNamelist() {
    int nerror = 0;
    char inputline[256], temp1[256], block[256], lhs[256], rhs[256], dummy[256], element[256];

    // read the input file
    FILE *inputfile = 0;
    std::string inputfilename = "namelist.ini";

    inputfile = fopen(inputfilename.c_str(), "r");
    if (inputfile == NULL) {
        std::printf("ERROR \"%s\" does not exist\n", inputfilename.c_str());
        ++nerror;
    }

    int n;
    bool blockset = false;
    int nErrors = 0;
    int nLines  = 0;
    int nLine;

    while (std::fgets(inputline, 256, inputfile) != NULL) {
        nLines++;
    }
    rewind(inputfile);

     // check the cases: comments, empty line, block, value, rubbish
    for (int nn=0; nn<nLines; nn++) {
        nLine = nn+1;

        // fetch a line and broadcast it
        std::fgets(inputline, 256, inputfile);

        // check for empty line
        n = std::sscanf(inputline, " %s ", temp1);
        if (n == 0)
            continue;

        // check for comments
        n = std::sscanf(inputline, " #%[^\n]", temp1);
        if (n > 0)
            continue;

        n = std::sscanf(inputline, " [%[^]]] ", temp1);
        if (n == 1) {
            n = std::sscanf(temp1, "%s %s", block, dummy);
            if (n == 1) {
                blockset = true;
            }
            else {
                std::printf("ERROR line %d: illegal block specification [%s]\n", nLine, temp1);
                return 1;
            }
            continue;
        }

        // read items
        n = std::sscanf(inputline, "%[^=] = %[^\n]", temp1, rhs);
        if (n == 2) {

            n = std::sscanf(temp1, " %[a-zA-Z0-9_()][%[^]]] %s", lhs, element, dummy);
            if (n <= 2) {
                if (!blockset) {
                    std::printf("ERROR line %d: illegal item [?][%s] = \"%s\"\n", nLine, lhs, rhs);
                    nErrors++;
                    return 1;
                }

                if (n ==1) {
                    std::strcpy(element,"default");
                }
                std::string blockstring(block);
                std::string itemstring(lhs);
                std::string elementstring(element);
                std::string valuestring(rhs);
                if (checkItemExists(blockstring, itemstring, elementstring)) {
                    inputList[blockstring][itemstring][elementstring].data   = valuestring;
                    inputList[blockstring][itemstring][elementstring].isused = false;
                }
                else {
                    std::printf("ERROR line %d: Item [%s][%s][%s] defined for the second time\n", nLine, block, lhs, element);
                    return 1;
                }
            }
            else {
                n = std::sscanf(inputline, "%[^=]", temp1);
                std::printf("ERROR line %d: illegal item  [%s][%s]\n", nLine, block, temp1);
                nErrors++;
            }
        }

        // throw exception
        else {
            n = std::sscanf(inputline, "%[^\n]", temp1);
            if (n > 0) {
                std::printf("ERROR line %d: \"%s\" is illegal input\n", nLine, temp1);
                nErrors++;
            }
        }
    }

    fclose(inputfile);

    return nErrors;

}

int Input::checkItemExists(std::string cat, std::string item, std::string el) {
    inputMap::const_iterator it1 = inputList.find(cat);

    bool readerror = false;

    if (it1 != inputList.end()) {
        inputMap2d::const_iterator it2 = it1->second.find(item);

        if (it2 != it1->second.end()) {
            inputMap1d::const_iterator it3 = it2->second.find(el);
            if (it3 == it2->second.end())
                readerror = true;
        }
        else
            readerror = true;
    }
    else
        readerror = true;

    if (readerror)
        return 1;

    return 0;
}

// int functions
int Input::getItem(int* value, std::string cat, std::string item, std::string el) {
    bool optional = false;
    int dummy = 0;

    if (parseItem(value, cat, item, el, optional, dummy))
        return 1;

    return 0;
}

int Input::getItem(int* value, std::string cat, std::string item, std::string el, int def) {
    bool optional = true;

    if (parseItem(value, cat, item, el, optional, def))
        return 1;

    return 0;
}

template <class valuetype>
int Input::parseItem(valuetype* value, std::string cat, std::string item, std::string el, bool optional, valuetype def) {
    std::string itemout, itemtype;
    itemout = "[" + cat + "][" + item + "]";

    if (!el.empty()) {
        itemout += "[" + el + "]";
        if (!checkItemExists(cat, item, el)) {
            if (checkItem(value, cat, item, el))
                return 1;
            //itemtype = "(element specific)";
        }
    }
    if (itemtype.empty()) {
        if (checkItemExists(cat, item)) {
            if (optional) {
                *value = def;
                //itemtype = "(default)";
            }
            else {
                std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
                return 1;
            }
        }
        else {
            if (checkItem(value, cat, item))
                return 1;
        }
    }
    //std::cout << std::left  << std::setw(30) << itemout << "= "
    //    << std::right << std::setw(30) << std::setprecision(5) << std::boolalpha << *value
    //    << "   " << std::endl;

    return 0;
}

int Input::checkItem(int* value, std::string cat, std::string item, std::string el) {
    char inputstring[256], temp[256];
    std::strcpy(inputstring, inputList[cat][item][el].data.c_str());

    int inputint;
    int n = std::sscanf(inputstring, " %d %[^\n] ", &inputint, temp);

    if (n == 1)
        *value = inputint;
    else {
        if (std::strcmp(inputstring,"")) {
            if (el == "default") {
                std::printf("ERROR [%s][%s] = \"%s\" is not of type INT\n", cat.c_str(), item.c_str(), inputstring);
            }
            else {
                std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type INT\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
            }
            return 1;
        }
    }
    inputList[cat][item][el].isused = true;

    return 0;
}

// double functions
int Input::getItem(double* value, std::string cat, std::string item, std::string el) {
    bool optional = false;
    double dummy = 0.;

    if (parseItem(value, cat, item, el, optional, dummy))
        return 1;

    return 0;
}

int Input::getItem(double* value, std::string cat, std::string item, std::string el, double def) {
    bool optional = true;

    if (parseItem(value, cat, item, el, optional, def))
        return 1;

    return 0;
}

int Input::checkItem(double* value, std::string cat, std::string item, std::string el) {
    char inputstring[256], temp[256];
    std::strcpy(inputstring, inputList[cat][item][el].data.c_str());

    double inputdouble;
    int n = std::sscanf(inputstring, " %lf %[^\n] ", &inputdouble, temp);
    // catch the situation where a double is closed with a ".", which is not read by sscanf's %f
    if (n == 1 || (n == 2 && !std::strcmp(".", temp)))
        *value = inputdouble;
    else {
        if (std::strcmp(inputstring,"")) {
            if (el == "default") {
                std::printf("ERROR [%s][%s] = \"%s\" is not of type DOUBLE\n", cat.c_str(), item.c_str(), inputstring);
            }
            else {
                std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type DOUBLE\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
            }
            return 1;
        }
    }
    inputList[cat][item][el].isused = true;

    return 0;
}

// strings
int Input::getItem(std::string* value, std::string cat, std::string item, std::string el) {
    bool optional = false;
    std::string dummy = "";

    if (parseItem(value, cat, item, el, optional, dummy))
        return 1;

    return 0;
}

int Input::getItem(std::string* value, std::string cat, std::string item, std::string el, std::string def) {
    bool optional = true;

    if (parseItem(value, cat, item, el, optional, def))
        return 1;

    return 0;
}

int Input::checkItem(std::string* value, std::string cat, std::string item, std::string el) {
    char inputstring[256], stringval[256], dummy[256];
    std::strcpy(inputstring, inputList[cat][item][el].data.c_str());

    int n = std::sscanf(inputstring, " %s %[^\n] ", stringval, dummy);

    if (n == 1)
        *value = stringval;
    else {
        if (std::strcmp(inputstring,"")) {
            if (el == "default") {
                std::printf("ERROR [%s][%s] = \"%s\" is not of type STRING\n", cat.c_str(), item.c_str(), inputstring);
            }
            else {
                std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type STRING\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
            }
            return 1;
        }
    }
    inputList[cat][item][el].isused = true;

    return 0;
}

// list retrieval function
int Input::getList(std::vector<std::string>* value, std::string cat, std::string item, std::string el) {
    if (parseList(value, cat, item, el))
        return 1;

    return 0;
}

template <class valuetype>
int Input::parseList(std::vector<valuetype>* value, std::string cat, std::string item, std::string el) {
    std::string itemout, listout;
    std::stringstream liststream;

    itemout = "[" + cat + "][" + item + "]";
    if (checkItemExists(cat, item)) {
        std::cout << std::left  << std::setw(30) << itemout << "= "
            << std::right << std::setw(30) << "EMPTY LIST" << std::endl;
    }
    else {
        if (checkList(value, cat, item))
            return 1;
        typedef typename std::vector<valuetype>::iterator itertype;
        for (itertype it = value->begin(); it !=value->end()-1; ++it) {
            liststream << *it << ", ";
        }
        liststream << *(value->end()-1);
        std::cout << std::left  << std::setw(30) << itemout << "= "
            << std::right << std::setw(30) << liststream.str() << std::endl;
    }

    return 0;
}

int Input::checkList(std::vector<std::string>* value, std::string cat, std::string item, std::string el) {
    char inputstring[256], dummy[256];
    std::strcpy(inputstring, inputList[cat][item][el].data.c_str());

    char temp1[256];
    char* temp2;
    std::strcpy(temp1, inputstring);

    // first, split string on the delimiter
    temp2 = std::strtok(temp1, ",");

    while (temp2 != NULL) {
        // read in the string part in temp1
        int n = std::sscanf(temp2, "%s %s", temp1, dummy);

        // store the contents in the vector, or throw exception
        if (n == 1)
            value->push_back(temp1);
        else {
            if (std::strcmp(inputstring,"")) {
                if (el == "default") {
                    std::printf("ERROR [%s][%s] = \"%s\" is not a list of type STRING\n", cat.c_str(), item.c_str(), inputstring);
                }
                else {
                    std::printf("ERROR [%s][%s][%s] = \"%s\" is not a list of type STRING\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
                }
                // empty the vector
                value->clear();
                return 1;
            }
        }

        // retrieve the next raw substring
        temp2 = std::strtok(NULL, ",");
    }
    inputList[cat][item][el].isused = true;

    return 0;
}

// read timeseries data
int Input::readDataFile(dataMap* series, std::string inputname, bool optional) {
    int nerror = 0;
    char inputline[256], temp1[256];
    char* substring;
    int n;

    // read the input file
    FILE* inputfile = 0;
    std::string inputfilename = inputname;

    int doreturn = 0;
    inputfile = fopen(inputfilename.c_str(), "r");
    if (inputfile == NULL) {
        if (optional)
            doreturn = true;
        else {
            std::printf("ERROR \"%s\" does not exist\n", inputfilename.c_str());
            nerror++;
        }
    }

    // broadcast the error count
    if (nerror) return 1;
    if (doreturn) return 0;

    int nlines = 0;
    int nline;
    int nvar = 0;
    std::vector<std::string> varnames;

    while (std::fgets(inputline, 256, inputfile) != NULL)
        nlines++;
    rewind(inputfile);
    int nn;

    // first find the header
    for (nn=0; nn<nlines; nn++) {
        nline = nn+1;
        // fetch a line and broadcast it
        std::fgets(inputline, 256, inputfile);

        // check for empty line
        n = std::sscanf(inputline, " %s ", temp1);
        if (n == 0) continue;

        // check for comments
        n = std::sscanf(inputline, " #%[^\n]", temp1);
        if (n > 0) continue;

        // read the header
        // read the first substring
        substring = std::strtok(inputline, " ,;\t\n");
        while (substring != NULL) {
            nvar++;

            // temporarily store the variable name
            varnames.push_back(std::string(substring));

            // read the next substring
            substring = std::strtok(NULL, " ,;\t\n");
        }

        if (nvar == 0) {
            std::printf("ERROR no variable names in header\n");
            fclose(inputfile);
            return 1;
        }

        // step out of the fgets loop
        break;
    }

    // second read the data, continue reading
    int ncols;
    double datavalue;

    std::vector<double> varvalues;

    // continue the loop from the exit value of nn
    for (nn++; nn<nlines; nn++) {
        nline = nn+1;
        
        // fetch a line and broadcast it
        std::fgets(inputline, 256, inputfile);

        // check for empty line
        n = std::sscanf(inputline, " %s ", temp1);
        if (n == 0) continue;

        // check for comments
        n = std::sscanf(inputline, " #%[^\n]", temp1);
        if (n > 0) continue;

        // read the data
        ncols = 0;
        varvalues.clear();
        
        // read the first substring
        substring = std::strtok(inputline, " ,;\t\n");
        while (substring != NULL) {
            ncols++;

            // scan the line, while checking that the whole string has been read
            n = std::sscanf(substring, " %lf %[^\n]", &datavalue, temp1);

            if (n != 1) {
                std::printf("ERROR line %d: \"%s\" is not a correct data value\n", nline, substring);
                fclose(inputfile);
                return 1;
            }

            // temporarily store the data
            varvalues.push_back(datavalue);

            // read the next substring
            substring = std::strtok(NULL, " ,;\t\n");
        }

        if (ncols != nvar) {
            std::printf("ERROR line %d: %d data columns, but %d defined variables\n", nline, ncols, nvar);
            fclose(inputfile);
            return 1;
        }

        // store the data
        for (n=0; n<nvar; n++)
            (*series)[varnames[n]].push_back(varvalues[n]);
    }

    fclose(inputfile);
    return 0;
}

// read input data
int Input::getProf(double* data, std::string inputType, std::string varname, int kmaxin) {
    
    dataMap::const_iterator it;
    dataMap inputData;
    
    if (inputType=="soil") {
        it = inputSoil.find(varname);
        inputData = inputSoil;
    }
    if (inputType=="metr") {
        it = inputMetr.find(varname);
        inputData = inputMetr;
    }
        
    if (it != inputData.end()) {
        int profsize = int(inputData[varname].size());
        if (profsize < kmaxin) {
            std::printf("ERROR only %d of %d levels can be read for variable \"%s\"\n", profsize, kmaxin, varname.c_str());
            return 1;
        }
        if (profsize > kmaxin)
           std::printf("WARNING %d is larger than the number of grid points %d for variable \"%s\"\n", profsize, kmaxin, varname.c_str());

        for (int k=0; k<kmaxin; k++)
            data[k] = inputData[varname][k];
    }
    else {
        std::printf("WARNING no profile data for variable \"%s\", values set to zero\n", varname.c_str());
        for (int k=0; k<kmaxin; k++)
            data[k] = 0.;
    }

    return 0;
}

int Input::getProf(std::vector<double>* data, std::string inputType, std::string varname, int kmaxin) {
    
    dataMap::const_iterator it;
    dataMap inputData;
    
    if (inputType=="soil") {
        it = inputSoil.find(varname);
        inputData = inputSoil;
    }
    if (inputType=="metr") {
        it = inputMetr.find(varname);
        inputData = inputMetr;
    }
        
    if (it != inputData.end()) {
        int profsize = int(inputData[varname].size());
        if (profsize < kmaxin) {
            std::printf("ERROR only %d of %d levels can be read for variable \"%s\"\n", profsize, kmaxin, varname.c_str());
            return 1;
        }
        if (profsize > kmaxin)
           std::printf("WARNING %d is larger than the number of grid points %d for variable \"%s\"\n", profsize, kmaxin, varname.c_str());

        for (int k=0; k<kmaxin; k++)
            data->push_back(inputData[varname][k]);
        data->resize(kmaxin);
    }
    else {
        std::printf("WARNING no profile data for variable \"%s\", values set to zero\n", varname.c_str());
        for (int k=0; k<kmaxin; k++)
            data->push_back(0.);
        data->resize(kmaxin);
    }

    return 0;
}

int Input::getProf(std::vector<int>* data, std::string inputType, std::string varname, int kmaxin) {
    
    dataMap::const_iterator it;
    dataMap inputData;
    
    if (inputType=="soil") {
        it = inputSoil.find(varname);
        inputData = inputSoil;
    }
    if (inputType=="metr") {
        it = inputMetr.find(varname);
        inputData = inputMetr;
    }
        
    if (it != inputData.end()) {
        int profsize = int(inputData[varname].size());
        if (profsize < kmaxin) {
            std::printf("ERROR only %d of %d levels can be read for variable \"%s\"\n", profsize, kmaxin, varname.c_str());
            return 1;
        }
        if (profsize > kmaxin)
           std::printf("WARNING %d is larger than the number of grid points %d for variable \"%s\"\n", profsize, kmaxin, varname.c_str());

        for (int k=0; k<kmaxin; k++)
            data->push_back(inputData[varname][k]);
        data->resize(kmaxin);
    }
    else {
        std::printf("WARNING no profile data for variable \"%s\", values set to zero\n", varname.c_str());
        for (int k=0; k<kmaxin; k++)
            data->push_back(0.);
        data->resize(kmaxin);
    }

    return 0;
}
