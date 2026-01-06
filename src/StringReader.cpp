#include "StringReader.h"
#include <stdexcept>
#include <iostream>
#include <sstream>

bool StringReader::findParameterInLine(const std::string &line, const std::string &parameterName)
{
    if (line.empty())
        return false;

    size_t index = line.find(parameterName);
    if (index != std::string::npos && (index == 0 || line[index-1] == ' ')) 
    {
        return true;
    }
    
    return false;
}

std::string StringReader::formatLine(const std::string &line)
{
    std::string newLine;
    newLine.reserve(line.size());  // Reserve space to avoid reallocations
    unsigned numberSpace = 0;
    
    // First pass: Remove extra spaces and tabs
    for (size_t i = 0; i < line.size(); i++) 
    {
        if (line[i] == ' ' || line[i] == '\t') 
        {
            numberSpace++;
            if (numberSpace == 1)
                newLine += ' ';
        }
        else 
        {
            newLine += line[i];
            numberSpace = 0;
        }
    }

    // Second pass: Remove spaces around '='
    std::string finalLine;
    finalLine.reserve(newLine.size());

    for (size_t i = 0; i < newLine.size(); i++) 
    {
        if (newLine[i] == ' ' && i + 1 < newLine.size() && newLine[i + 1] == '=')
        {
            continue;  // Skip space before '='
        }
        else if (newLine[i] == '=' && i > 0 && finalLine.back() == ' ')
        {
            finalLine.back() = '=';  // Replace previous space with '='
        }
        else if (newLine[i] == '=' && i+1 < newLine.size() && newLine[i+1] == ' ') 
        {
            finalLine += '=';
            i++;
        }
        else
        {
            finalLine += newLine[i];
        }
    }

    return finalLine;
}

bool StringReader::getDoubleParameter(const std::string &line, std::string parameterName, double &val)
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = std::stod(line.substr(line.find(parameterName) + parameterName.size()));
        }catch (const std::invalid_argument & e) {   
            return false;                                 
        } catch (const std::out_of_range & e) {      
            return false;                                   
        }        
    else
        return false;
    return true;
}

bool StringReader::getIntParameter(const std::string &line, std::string parameterName, int &val) 
{
    if (findParameterInLine(line, parameterName))
        try {
            val = std::stoi(line.substr(line.find(parameterName) + parameterName.size()));
        }catch (const std::invalid_argument & e) { 
            return false;                                
        } catch (const std::out_of_range & e) {   
            return false;                                       
        }     
    else
        return false;
    return true;
}

bool StringReader::getUnsignedLLIntParameter(const std::string &line, std::string parameterName, unsigned long long int &val) 
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = std::stoull(line.substr(line.find(parameterName) + parameterName.size()));
        }catch (const std::invalid_argument & e) { 
            return false;                                    
        } catch (const std::out_of_range & e) {    
            return false;                                  
        }     
    else
        return false;
    return true;    
}

bool StringReader::getUnsignedParameter(const std::string &line, std::string parameterName, unsigned &val) 
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = 0;
            unsigned long int valL = 0;
            valL = std::stoul(line.substr(line.find(parameterName) + parameterName.size()));
            val = static_cast <unsigned> (valL);
        }catch (const std::invalid_argument & e) {    
            return false;                                
        } catch (const std::out_of_range & e) {   
            return false;                                  
        }     
    else
        return false;
    return true;
}

bool StringReader::getLineParameter(const std::string &line, std::string parameterName, std::string &val)
{
    if(findParameterInLine(line, parameterName)) 
    {
        try {
            val = "";
            val = line.substr(line.find(parameterName) + parameterName.size());
        } catch (const std::invalid_argument & e) {
            return false;
        } catch (const std::out_of_range & e) {
            return false;
        }
    }
    else
        return false;
    return true;
}
