//
//  input.cpp
//  
//  This class handles reading in user options
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

Input :: Input() {
    
    int nError = 0;
    nError += readNamelist();

    if (nError)
        throw 1;
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

    std::cout<<"Processing ini file "<< inputfilename << std::endl;
    std::cout<<"##############################################################"<<std::endl;
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
    std::cout << std::left  << std::setw(30) << itemout << "= "
        << std::right << std::setw(30) << std::setprecision(5) << std::boolalpha << *value
        << "   " << std::endl;

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