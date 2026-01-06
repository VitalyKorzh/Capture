#ifndef __COUNTER_H__
#define __COUNTER_H__

#include <vector>

#include "InputReader.h"

typedef std::vector<double> darray;
typedef std::vector<unsigned> uiarray;
typedef unsigned uint;


class Counter {
private:

    const InputReader reader;
    std::ostream &os;
    const uint nz;
    const uint nr;
    const darray &ni;
    const darray &zArray;
    const darray &rArray;
    const uint nParticles;
    const double sigma;
    const double theta;

    const darray &sArray;
    const uint ns;

    const std::pair<double, double> &position;
    uiarray nCap;
    uint nFlyby;

    void clearPrevious();

public:
    Counter(std::istream &in=std::cin, std::ostream &os=std::cout);
    void count();
    
    double getSigma() const { return sigma; }
    uint getNParticles() const { return nParticles; }
    uint getNz() const { return nz; }
    uint getNFlyply() const { return nFlyby; }
    const darray & getZArray() const { return zArray; }
    const darray & getNi() const { return ni; }
    const uiarray & getNCap() const { return nCap; }


    double getnCap(uint index) const { return ((double) nCap[index]) / nParticles; }
    double getnFlyby() const { return ((double) nFlyby) / nParticles; }

    bool isReadSuccess() const { return reader.work; }
    const InputReader getReader() const { return reader; }

    void printStartInfo() const;
    void printResult() const;

    ~Counter();

};

#endif