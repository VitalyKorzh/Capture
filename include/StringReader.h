#ifndef _STRING_READER_H_
#define _STRING_READER_H_

#include <string>

class StringReader 
{
private:
    static bool findParameterInLine(const std::string &line, const std::string &parameterName);
public:
    static std::string formatLine(const std::string &line); //убрать лишнии пробелы и табуляции в строке, заменить табуляции на пробелы, убрать пробелы между знаком =
    static bool getDoubleParameter(const std::string &line, std::string parameterName, double &val);
    static bool getIntParameter(const std::string &line, std::string parameterName, int &val);
    static bool getUnsignedLLIntParameter(const std::string &line, std::string parameterName, unsigned long long int &val);
    static bool getUnsignedParameter(const std::string &line, std::string parameterName, unsigned &val);
    static bool getLineParameter(const std::string &line, std::string parameterName, std::string &val);
    virtual ~StringReader()=0;
};

#endif