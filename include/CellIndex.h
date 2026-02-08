#ifndef __CELL_INDEX_H__
#define __CELL_INDEX_H__

typedef unsigned uint;

struct CellIndex
{
    uint iZ;
    uint iR;
    uint iPhi;

    CellIndex(uint iZ, uint iR, uint iPhi) : iZ(iZ), iR(iR), iPhi(iPhi)
    {}

    uint getIndex(uint nz, uint nr) const { return iR + nr*iZ + nr*nz*iPhi; }
};


#endif