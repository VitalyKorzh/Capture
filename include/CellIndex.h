#ifndef __CELL_INDEX_H__
#define __CELL_INDEX_H__

typedef unsigned uint;

struct CellIndex
{
    uint iZ;
    uint iR;
    uint iPhi;

    CellIndex(uint iZ=-1, uint iR=-1, uint iPhi=-1) : iZ(iZ), iR(iR), iPhi(iPhi)
    {}

    uint getIndex(uint nz, uint nr) const { return iR + nr*iZ + nr*nz*iPhi; }

    bool errorIndex() const { return iZ == (uint) -1 || iR == (uint) -1 || iPhi == (uint)-1; }
};


#endif