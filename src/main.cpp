#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "Counter.h"
#include "Command.h"

int main(int argc, char** argv)
{
    std::string inputFile = "../test.in";
    std::string outputFile = "../test.out";


    bool work = commands(inputFile, outputFile, argc, argv);

    if (!work)
        return 1;

    std::ifstream fin(inputFile);
    std::ofstream fout(outputFile);

    if (fin.is_open() && fout.is_open())
    {
        Counter counter(fin, fout);
        if (!counter.isReadSuccess()) {
            std::cerr << counter.getReader().getError();
            fin.close();
            fout.close();
            return 1;
        }
        counter.printStartInfo();
        counter.count();
        counter.printResult();
    }
    fin.close();
    fout.close();

    return 0;

}

