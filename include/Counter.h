#ifndef __COUNTER_H__
#define __COUNTER_H__

#include <vector>

#include "InputReader.h"
#include "Line.h"
#include "Injector.h"

typedef std::vector<double> darray;
typedef std::vector<unsigned> uiarray;
typedef unsigned uint;


class Counter {
private:

    const InputReader reader;
    std::ostream &os;
    const uint nz;
    const uint nr;
    const uint nphi;
    const darray &ni;
    const darray &zArray;
    const darray &rArray;
    const darray &phiArray;
    const std::vector <Injector> &injectors;


    uint nParticles;
    barray intersectionCell;
    uiarray nCap;
    uint nFlyby;

    void clearPrevious();

    uint getIndex(uint iz, uint ir, uint iphi) const { return ir + nr*iz + nr*nz*iphi; }

    void process(const Line &line, const Injector &injector, double gamma);

public:
    Counter(std::istream &in=std::cin, std::ostream &os=std::cout);
    void count();
    
    uint getNz() const { return nz; }
    uint getNr() const { return nr; }
    uint getNPhi() const { return nphi; }
    uint getNFlyply() const { return nFlyby; }
    const darray & getZArray() const { return zArray; }
    const darray & getRArray() const { return rArray; }
    const darray & getPhiArray() const { return phiArray; }
    const darray & getNi() const { return ni; }
    const uiarray & getNCap() const { return nCap; }
    const barray & getIntersectionCell() const { return intersectionCell; }


    double getnCap(uint index) const { return ((double) nCap[index]) / nParticles; }
    double getnFlyby() const { return ((double) nFlyby) / nParticles; }

    bool isReadSuccess() const { return reader.work; }
    const InputReader getReader() const { return reader; }

    void printStartInfo() const;
    void printResult() const;

    ~Counter();

};

#endif