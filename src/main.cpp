#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "Counter.h"

int main(int argc, char** argv)
{

    std::ifstream fin("../test.in");
    std::ofstream fout("../test.txt");

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

