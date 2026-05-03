#include "Line.h"
#include <cmath>
#include <algorithm>
#include <iostream>

inline double Line::getRPoint(double t, int l) const
{ 
    double delta = (d - t * sinTheta);
    double r = rho + l*sign*rhoLarmor;
    return sqrt(r*r + delta*delta); 
}

inline double Line::getRLarmorPoint(double r) const
{
    double r0 = rho + sign*rhoLarmor;
    return sqrt(r0*r0+r*r-rho*rho);
}

inline double Line::getPhiPoint(double t, int l) const
{
    double y = getYPoint(t, l);
    double x = getXPoint(t, l);

    if (x > 0 && y >= 0)
        return atan(y/x);
    else if (x > 0 && y < 0)
        return 2.*M_PI + atan(y/x);
    else if (x < 0)
        return M_PI + atan(y/x);
    else if (x == 0 && y > 0)
        return M_PI_2;
    else if (x == 0 && y < 0)
        return 3. * M_PI_2;
    else
        return 0.0;
}

Line::Line(
    double rho, double theta, double phi0_rho, double z0_rho,
    uint nz, uint nr, uint nphi,
    const darray &zArray, const darray &rArray, const darray &phiArray, double t_epsilon, 
    double t_epsilon_first, bool plusDirection, double rhoLarmor, double rMax, bool count_centers) :
                                                            rho(rho),
                                                            theta(theta), phi0_rho(phi0_rho), z0_rho(z0_rho),
                                                            sinTheta(sin(theta)), cosTheta(cos(theta)), sinPhi0(sin(phi0_rho)),
                                                            cosPhi0(cos(phi0_rho)), sign(plusDirection ? 1 : -1), sx(-sinTheta * sinPhi0*sign),
                                                            sy(sinTheta * cosPhi0*sign), sz(cosTheta), 
                                                            //intersectionPoints(nz * nr * nphi, false),
                                                            //intescetionPointsCenter(nz*nr*nphi, false),
                                                            nz(nz), nr(nr), nphi(nphi), crossAxis(false), t_epsilon(t_epsilon), 
                                                            t_epsilon_first(t_epsilon_first), rhoLarmor(rhoLarmor), count_centers(count_centers),
                                                            lineWork(false)
{
    if (nz == 0 || nz == 0 || nphi == 0)
        return;

    if (theta <= 0. || rho < 0)
        return;

    Rmax = rArray.back();
    iR_crit = nr-1;
    if (rMax > 0. && rMax < Rmax)
    {
        iR_crit = findIndex(rArray, nr, rMax, false);
        Rmax = rArray[iR_crit+1];
    }

    zmin = zArray.front();
    zmax = zArray.back();

    if (Rmax <= rho) {
        std::cerr << "вышла за радиальную границу Rmax <= rho\n";
        return;
    }

    d = sqrt(Rmax*Rmax - rho*rho);
    t_crit = d / sinTheta;

    x00 = rho*cosPhi0 + d * sinPhi0*sign;
    y00 = rho*sinPhi0 - d * cosPhi0*sign;
    z00 = z0_rho - d * cosTheta / sinTheta;

    if (z00 < zmin || z00 >= zmax) {
        std::cerr << "вышли за левую границу оси z\n";
        return;
    }

    if (z0_rho + d * cosTheta / sinTheta > zmax) {
        std::cerr << "вышли за правую границу оси z\n";
        return;
    }

    lineWork = createDataArray(nz, nr, nphi, zArray, rArray, phiArray);

    // if (lineWork)
    //     for (const LineData &id : data) // заполнить точки пересечения
    //     {
    //         intersectionPoints[id.index.getIndex(nz, nr)] = true;
    //         if (!id.indexCenter.errorIndex())
    //             intescetionPointsCenter[id.indexCenter.getIndex(nz, nr)] = true;
    //     }
}

inline double Line::getTR(double r, double tPrevious, int l) const
{
    double t1 = -1.;
    double t2 = -1.;

    double rho_l = rho + l*sign*rhoLarmor;

    if (r*r >= rho_l*rho_l) {
        double sq = sqrt(r*r-rho_l*rho_l);
        t1 = 1./sinTheta * (d - sq);
        t2 = 1./sinTheta * (d + sq);
    }

    return t_crit > tPrevious ? t1 : t2;
}

inline double Line::getTZ(double z, int l) const
{
    if (std::abs(cosTheta) > 1e-12)
        return (z - z00) / cosTheta;
    else
        return -1.; // очень большой луч не пересикает
}

inline double Line::getTPhi(double phi, int l) const
{
    double tanPhi = tan(phi);

    double x00_l = getX0L(l);
    double y00_l = getY0L(l);

    if (phi != M_PI_2 && phi != 3. * M_PI_2)
        return (tanPhi*x00_l-y00_l) / (sy - sx*tanPhi);
    else if (sx != 0.)
        return - x00_l / sx;
    else
        return -1.; //если nan то ломается сортировка
}

void Line::getNewIndex(const Boundary &boundary, CellIndex &index)
{
    if (boundary.type == Boundary::IntersectionType::Z)
        index.iZ = boundary.index;
    else if (boundary.type == Boundary::IntersectionType::R)
    {
        index.iR = boundary.index;
    }
    else if (!crossAxis)
    {
        index.iPhi = boundary.index;
        if (index.iR == nr+1)
            index.iPhi = nphi-1;
    }

}

bool Line::checkBoundary(const std::vector<Boundary> &boundaryArray, const Boundary &boundary, uint index) const
{
    for (uint ii = index; ii < boundaryArray.size(); ii++)
    {
        if (boundary == boundaryArray[ii])
            return true;
    }

    return false;
}

bool Line::checkBoundaryR(const std::vector<Boundary> &boundaryArray, Boundary::IntersectionDirection direction, bool isCenter, uint index) const
{
    for (uint ii = index; ii < boundaryArray.size(); ii++)
    {
        if (boundaryArray[ii].type == Boundary::IntersectionType::R && boundaryArray[ii].direction == direction && boundaryArray[ii].isCenter == isCenter)
            return true;
    }
    return false;
}

LineData Line::traceLine(uint nz, uint nr, uint nphi, const darray &zArray, const darray &rArray, const darray &phiArray, double &tPrevious,
                         CellIndex &index, std::vector<Boundary> &bPrevious)
{
    //double z_current = z00 + tPrevious * sz;
    //double r_current = getRPoint(tPrevious);
    double phi_current = getPhiPoint(tPrevious);

    double s = 0;

    double r1 = rArray[index.iR];
    double r2 = rArray[index.iR+1];
    double z1 = zArray[index.iZ];
    double z2 = zArray[index.iZ+1];
    double phi1 = phiArray[index.iPhi];
    double phi2 = phiArray[index.iPhi+1];

    CellIndex indexPrev = index;

    double t_r1 = getTR(r1, tPrevious);
    double t_r2 = getTR(r2, t_r1 >= 0. ? tPrevious : t_crit+1.);

    std::vector<Boundary> boundaries= 
    {
        Boundary(getTZ(z1), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MIN, index.iZ > 0 ? index.iZ-1 : nz),
        Boundary(getTZ(z2), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MAX, index.iZ+1),
        Boundary(t_r1, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MIN, index.iR > 0 ? index.iR-1: nr+1),
        Boundary(t_r2, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MAX, index.iR+1),
        Boundary(getTPhi(phi1), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MIN, index.iPhi > 0 ? (index.iPhi-1) % nphi : nphi-1),
        Boundary(getTPhi(phi2), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MAX, (index.iPhi+1) % nphi )
    };

    if (t_r1 < 0. && std::abs(tPrevious - t_crit) <= t_epsilon_first && checkBoundaryR(bPrevious, Boundary::IntersectionDirection::MIN))
      index.iR++;
    
    std::sort(boundaries.begin(), boundaries.end(), [](Boundary a, Boundary b) 
        {
            return a.t_boundary < b.t_boundary;
        }
    );

    bool findBoundary = false;

    uint index_boundaries = 0;

    for (Boundary boundary : boundaries)
    {
        double t = boundary.t_boundary;

        if (t >= tPrevious && !checkBoundary(bPrevious, boundary, index_boundaries))
        {
            if (t > tPrevious + t_epsilon_first && !findBoundary)
            {
                s = (t - tPrevious);
                tPrevious = t;
                findBoundary = true;
                getNewIndex(boundary, index);
    
                bPrevious[index_boundaries] = boundary;
                index_boundaries++;
            }
            else if ((t - tPrevious) <= t_epsilon && findBoundary)
            {
                getNewIndex(boundary, index);
                s += (t-tPrevious);
                tPrevious = t;
                
                bPrevious[index_boundaries] = boundary;
                index_boundaries++;
            }
        }
    }

    for (uint ii = index_boundaries; ii < bPrevious.size(); ii++)
        bPrevious[ii].type = Boundary::IntersectionType::None;

    LineData data(indexPrev, s);

    if (index.iR == nr+1) // луч прошел через ось
    {
        crossAxis = true;
        index.iR = 0;
        double phi_new = phi_current+M_PI >= 2. * M_PI ? phi_current - M_PI : phi_current + M_PI;
        index.iPhi = findIndex(phiArray, nphi, phi_new);
    }
    else
        crossAxis = false;

    return data;
}

LineData Line::traceLineCenter(uint nz, uint nr, uint nphi, 
                        const darray &zArray, const darray &rArray, const darray &phiArray, 
                        double &tPrevious, CellIndex &index, CellIndex &indexCenter, 
                        std::vector <Boundary> &bPrevious)
{
    //double r_current = getRPoint(tPrevious);
    //double r_current_c = getRPoint(tPrevious, 1);
    //double z_current = getZPoint(tPrevious);
    double phi_current = getPhiPoint(tPrevious);
    double phi_current_c = getPhiPoint(tPrevious, 1);
    double s = 0.;
    
    double r1 = rArray[index.iR];
    double r2 = rArray[index.iR+1];
    double z1 = zArray[index.iZ];
    double z2 = zArray[index.iZ+1];
    double phi1 = phiArray[index.iPhi];
    double phi2 = phiArray[index.iPhi+1];
    
    double r1_c = rArray[indexCenter.iR];
    double r2_c = rArray[indexCenter.iR+1];
    double z1_c = zArray[indexCenter.iZ];
    double z2_c = zArray[indexCenter.iZ+1];
    double phi1_c = phiArray[indexCenter.iPhi];
    double phi2_c = phiArray[indexCenter.iPhi+1];
    
    CellIndex indexPrev = index;
    CellIndex indexCenterPrev = indexCenter;

    double t_r1 = getTR(r1, tPrevious);
    double t_r2 = getTR(r2, t_r1 >= 0. ? tPrevious : t_crit+1.);
    double t_r1_c = getTR(r1_c, tPrevious, 1);
    double t_r2_c = getTR(r2_c, t_r1_c >= 0. ? tPrevious : t_crit+1., 1);

    std::vector<Boundary> boundaries= 
    {
        Boundary(getTZ(z1), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MIN, index.iZ > 0 ? index.iZ-1 : nz),
        Boundary(getTZ(z2), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MAX, index.iZ+1),
        Boundary(t_r1, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MIN, index.iR > 0 ? index.iR-1: nr+1),
        Boundary(t_r2, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MAX, index.iR+1),
        Boundary(getTPhi(phi1), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MIN, index.iPhi > 0 ? (index.iPhi-1) % nphi : nphi-1),
        Boundary(getTPhi(phi2), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MAX, (index.iPhi+1) % nphi ),

        Boundary(getTZ(z1_c, 1), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MIN, indexCenter.iZ > 0 ? indexCenter.iZ-1 : nz, true),
        Boundary(getTZ(z2_c, 1), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MAX, indexCenter.iZ+1, true),
        Boundary(t_r1_c, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MIN, indexCenter.iR > 0 ? indexCenter.iR-1: nr+1, true),
        Boundary(t_r2_c, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MAX, indexCenter.iR+1, true),
        Boundary(getTPhi(phi1_c, 1), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MIN, 
        indexCenter.iPhi > 0 ? (indexCenter.iPhi-1) % nphi : nphi-1, true),
        Boundary(getTPhi(phi2_c, 1), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MAX, (indexCenter.iPhi+1) % nphi, true)
    };

    { //нужно акуратно развернуть луч
        if (t_r1 < 0. && std::abs(tPrevious - t_crit) <= t_epsilon_first && checkBoundaryR(bPrevious, Boundary::IntersectionDirection::MIN))
            index.iR++;
        if (t_r1_c < 0. && std::abs(tPrevious - t_crit) <= t_epsilon_first && checkBoundaryR(bPrevious, Boundary::IntersectionDirection::MIN, true))
            indexCenter.iR++;
    }
    
    std::sort(boundaries.begin(), boundaries.end(), [](Boundary a, Boundary b) 
        {
            return a.t_boundary < b.t_boundary;
        }
    );

    bool findBoundary = false;
    uint index_boundaries = 0;

    for (Boundary boundary : boundaries)
    {
        double t = boundary.t_boundary;

        if (t >= tPrevious && !checkBoundary(bPrevious, boundary, index_boundaries))
        {
            if (t > tPrevious + t_epsilon_first && !findBoundary)
            {
                s = (t - tPrevious);
                tPrevious = t;
                findBoundary = true;
                getNewIndex(boundary, boundary.isCenter ? indexCenter : index);
    
                bPrevious[index_boundaries] = boundary;
                index_boundaries++;
            }
            else if ((t - tPrevious) <= t_epsilon && findBoundary)
            {
                getNewIndex(boundary, boundary.isCenter ? indexCenter : index);
                s += (t-tPrevious);
                tPrevious = t;
                
                bPrevious[index_boundaries] = boundary;
                index_boundaries++;
            }
        }
    }

    for (uint ii = index_boundaries; ii < bPrevious.size(); ii++)
        bPrevious[ii].type = Boundary::IntersectionType::None;

    LineData data(indexPrev, s, indexCenterPrev);

    if (index.iR == nr+1) // луч прошел через ось
    {
        crossAxis = true;
        index.iR = 0;
        double phi_new = phi_current+M_PI >= 2. * M_PI ? phi_current - M_PI : phi_current + M_PI;
        index.iPhi = findIndex(phiArray, nphi, phi_new);
    }
    else if (indexCenter.iR == nr+1)
    {
        crossAxis = true;
        indexCenter.iR = 0;
        double phi_new = phi_current_c+M_PI >= 2. * M_PI ? phi_current_c - M_PI : phi_current_c + M_PI;
        indexCenter.iPhi = findIndex(phiArray, nphi, phi_new);
    }
    else
        crossAxis = false;

    return data;
}

bool Line::createDataArray(uint nz, uint nr, uint nphi, 
                            const darray &zArray, const darray &rArray, 
                            const darray &phiArray)
{
    double tPrevious = 0.;
    crossAxis = false;
    CellIndex index(0, iR_crit, 0);
    CellIndex indexCenter(nz, nr, nphi);

    { //поиск начального положения луча
        index.iZ = findIndex(zArray, nz, z00);
        index.iPhi = findIndex(phiArray, nphi, getPhiPoint(0.));
        index.iR = findIndex(rArray, nr, Rmax, false);
    }

    if (count_centers) // поиск начального положения луча центров
    {
        indexCenter.iZ = findIndex(zArray, nz, z00);
        indexCenter.iPhi = findIndex(phiArray, nphi, getPhiPoint(0., 1));
        indexCenter.iR = findIndex(rArray, nr, getRLarmorPoint(Rmax), false);

        if (indexCenter.errorIndex()) // пока не будет работать если вышла за границу
        {
            std::cerr << "ларморовский центр не попал в ячейку\n";
            return false;
        }
    }

    ns = 0;

    const uint N_PREV_BOUNDARIES = 4;
    uint size_coeff = count_centers ? 2 : 1;
    std::vector <Boundary> previousBoundaries(size_coeff*N_PREV_BOUNDARIES, Boundary());
    data.reserve(size_coeff*4*nr);
    while (index.iZ != nz && index.iR != iR_crit+1) {
        LineData data0 = count_centers ?
            traceLineCenter(nz, nr, nphi, zArray, rArray, phiArray, tPrevious, index, indexCenter, previousBoundaries) :
            traceLine(nz, nr, nphi, zArray, rArray, phiArray, tPrevious, index, previousBoundaries);
        data.push_back(data0);
        ns++;
        if (ns > nr*nz*nphi*size_coeff) {
            std::cerr << "вышли за предел числа итераций!\n";
            return false;
        }
    }

    return true;
}

uint Line::findIndex(const darray &array, uint n, double value, bool left) const
{
    if (left)
        for (uint i = 0; i < n; i++)
        {
            if (value >= array[i] && value < array[i+1])
                return i;
        }
    else
        for (uint i = 0; i < n; i++)
        {
            if (value > array[i] && value <= array[i+1])
                return i;
        }
    return -1;
}
