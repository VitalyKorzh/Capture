#ifndef __COMMAND_H__
#define __COMMAND_H__

#include <string>


void help();
bool commands(std::string &inputFile, std::string &outputFile, int argc, char**argv);


#endif