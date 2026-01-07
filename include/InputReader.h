#ifndef __INPUT_READER_H__
#define __INPUT_READER_H__

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <istream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <map>
#include <list>
#include <unordered_map>

#include "StringReader.h"

typedef std::vector<double> darray;
typedef std::vector<unsigned> uiarray;
typedef unsigned uint;

class InputReader
{
private:
    friend class Counter;
    bool work;
    uint numberLine;
    std::string error_message;

    uint precision;

    darray zArray;
    darray rArray;
    darray ni;
    uint nz;
    uint nr;
    uint nParticles;
    double sigma;
    double theta;
    std::pair<double, double> position;


    darray sArray;
    std::vector <std::pair<uint, uint>> index;
    uint ns = 0;

    double normaDensity;

    uint countSpace(std::string line) const;

    std::istream & getline(std::istream &in, std::string &line, bool formatLine=false, bool ignoreEqual=false) 
    {
        numberLine++;
        std::getline(in, line);
        if (formatLine)
            line = StringReader::formatLine(line);
        if (ignoreEqual)
            std::replace(line.begin(), line.end(), '=', ' ');
        return in;
    }
    void skip(std::istream &in, std::string &line, bool ignoreEqual=false);
    bool isComment(const std::string &line) const 
    {
        if (line.empty())
            return true;

        uint i = 0;
        while (line[i] == ' ' || line[i] == '\t') 
        {
            i++;
            if (i == line.size())
                return false;
        }
        return line[i] == '#'; 
    }

    std::string readWord(std::string line) 
    {
        std::istringstream iss(line);
        std::string word;
        iss >> word;
        return word;
    }

    void errorMessage(std::string error);
    void errorConfigConstNumberPar(std::string part1, const std::vector<std::string> PAR_NAMES, const bool *array, const uint N_STEP);

    bool isLine(const std::string &line, const std::string &name) const { return line.find(name) != std::string::npos && line.find(name + " end") == std::string::npos; } 

    bool readPosition(std::istream &in, std::pair<double, double> &p);
    bool readAxis(std::istream &in, darray &axis, uint &size, const std::string &name);
    bool readMesh(std::istream &in);
    bool readCount(std::istream &in);

    bool generateInjectionLine();

    bool checkArray(bool *array, const uint N_PAR)
    {

        for (uint i = 0; i < N_PAR; i++)
        {
            if (!array[i])
                return false;
        }

        return true;
    }

    void arrayBit(bool &array, bool test) const {
        array = test | array;
    }

public:
    
    InputReader(std::istream &in=std::cin);

    bool isWork() const { return work; }
    const std::string getError() const { return error_message; } 
    uint getPrecision() const { return precision; }

};

#endif