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
    const darray &zArray, const darray &rArray, const darray &phiArray, double t_epsilon, double t_epsilon_first, bool plusDirection,
                                                                double rhoLarmor, double rMax) : rho(rho),
                                                                          theta(theta), phi0_rho(phi0_rho), z0_rho(z0_rho),
                                                                          sinTheta(sin(theta)), cosTheta(cos(theta)), sinPhi0(sin(phi0_rho)),
                                                                          cosPhi0(cos(phi0_rho)), sign(plusDirection ? 1 : -1), sx(-sinTheta * sinPhi0*sign),
                                                                          sy(sinTheta * cosPhi0*sign), sz(cosTheta), intersectionPoints(nz * nr * nphi, false),
                                                                          nz(nz), nr(nr), nphi(nphi), crossAxis(false), t_epsilon(t_epsilon), 
                                                                          t_epsilon_first(t_epsilon_first), rhoLarmor(rhoLarmor), lineWork(false)
                                                                        

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

    if (Rmax <= rho)
        return;

    d = sqrt(Rmax*Rmax - rho*rho);
    t_crit = d / sinTheta;

    x00 = rho*cosPhi0 + d * sinPhi0*sign;
    y00 = rho*sinPhi0 - d * cosPhi0*sign;
    z00 = z0_rho - d * cosTheta / sinTheta;

    if (z00 < zmin || z00 >= zmax)
        return;

    if (z0_rho + d * cosTheta / sinTheta > zmax)
        return;

    lineWork = createDataArray(nz, nr, nphi, zArray, rArray, phiArray);

    if (lineWork)
        for (LineData id : data) // заполнить точки пересечения
            intersectionPoints[id.iR + nr*id.iZ + nz*nr*id.iPhi] = true;

}

inline double Line::getTR(double r, double tPrevious, int l) const
{
    double t1 = -1.;
    double t2 = -1.;

    double rho_l = rho + l*sign*rhoLarmor;

    if (r >= rho_l) {
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
    else
        return - x00_l / sx;
}

void Line::getNewIndex(const Boundary &boundary, uint &iZ, uint &iR, uint &iPhi)
{
    if (boundary.type == Boundary::IntersectionType::Z)
        iZ = boundary.index;
    else if (boundary.type == Boundary::IntersectionType::R)
    {
        iR = boundary.index;
    }
    else if (!crossAxis)
    {
        iPhi = boundary.index;
        if (iR == nr+1)
            iPhi = nphi-1;
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

bool Line::checkBoundaryR(const std::vector<Boundary> &boundaryArray, Boundary::IntersectionDirection direction, uint index) const
{
    for (uint ii = index; ii < boundaryArray.size(); ii++)
    {
        if (boundaryArray[ii].type == Boundary::IntersectionType::R && boundaryArray[ii].direction == direction)
            return true;
    }
    return false;
}

LineData Line::traceLine(uint nz, uint nr, uint nphi, const darray &zArray, const darray &rArray, const darray &phiArray, double &tPrevious,
                         uint &iZ, uint &iR, uint &iPhi, std::vector<Boundary> &bPrevious)
{

    // double x_current = x00 + tPrevious * sx;
    // double y_current = y00 + tPrevious * sy;

    double z_current = z00 + tPrevious * sz;
    double r_current = getRPoint(tPrevious);
    double phi_current = getPhiPoint(tPrevious);

    double s = 0;

    double r1 = rArray[iR];
    double r2 = rArray[iR+1];
    double z1 = zArray[iZ];
    double z2 = zArray[iZ+1];
    double phi1 = phiArray[iPhi];
    double phi2 = phiArray[iPhi+1];

    const uint iZprev = iZ;
    const uint iRprev = iR;
    const uint iPhiprev = iPhi;

    double t_r1 = getTR(r1, tPrevious);
    double t_r2 = getTR(r2, t_r1 >= 0. ? tPrevious : t_crit+1.);

    std::vector<Boundary> boundaries= 
    {
        Boundary(getTZ(z1), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MIN, iZ > 0 ? iZ-1 : nz),
        Boundary(getTZ(z2), Boundary::IntersectionType::Z, Boundary::IntersectionDirection::MAX, iZ+1),
        Boundary(t_r1, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MIN, iR > 0 ? iR-1: nr+1),
        Boundary(t_r2, Boundary::IntersectionType::R, Boundary::IntersectionDirection::MAX, iR+1),
        Boundary(getTPhi(phi1), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MIN, iPhi > 0 ? (iPhi-1) % nphi : nphi-1),
        Boundary(getTPhi(phi2), Boundary::IntersectionType::Phi, Boundary::IntersectionDirection::MAX, (iPhi+1) % nphi )
    };

    if (t_r1 < 0. && std::abs(tPrevious - t_crit) <= t_epsilon_first && checkBoundaryR(bPrevious, Boundary::IntersectionDirection::MIN))
      iR++;
    
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
                getNewIndex(boundary, iZ, iR, iPhi);
    
                bPrevious[index_boundaries] = boundary;
                index_boundaries++;
            }
            else if ((t - tPrevious) <= t_epsilon && findBoundary)
            {
                getNewIndex(boundary, iZ, iR, iPhi);
                s += (t-tPrevious);
                tPrevious = t;
                
                bPrevious[index_boundaries] = boundary;
                index_boundaries++;
            }
        }
    }

    for (uint ii = index_boundaries; ii < bPrevious.size(); ii++)
        bPrevious[ii].type = Boundary::IntersectionType::None;

    LineData data(iZprev, iRprev, iPhiprev, s);

    if (iR == nr+1) // луч прошел через ось
    {
        crossAxis = true;
        iR = 0;
        double phi_new = phi_current+M_PI >= 2. * M_PI ? phi_current - M_PI : phi_current + M_PI;
        iPhi = findIndex(phiArray, nphi, phi_new);
        //for (uint iphi = 0; iphi < nphi; iphi++)
        //   if (phi_new >= phiArray[iphi] && phi_new < phiArray[iphi+1])
        //       iPhi = iphi;
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
    uint iR = iR_crit;
    uint iZ = 0;
    uint iPhi = 0;

    {
        iZ = findIndex(zArray, nz, z00);
        iPhi = findIndex(phiArray, nphi, getPhiPoint(0.));
        iR = findIndex(rArray, nr, Rmax, false);
        // for (uint iz = 0; iz < nz; iz++)
        // {
        //     if (z00 >= zArray[iz] && z00 < zArray[iz+1])
        //         iZ = iz;
        // }
            
            
        // double phi00 = getPhiPoint(0.);
        // for (uint iphi = 0; iphi < nphi; iphi++)
        // {
        //     if (phi00 >= phiArray[iphi] && phi00 < phiArray[iphi+1])
        //         iPhi = iphi;
        // }
    }

    ns = 0;

    const uint N_PREV_BOUNDARIES = 4;
    std::vector <Boundary> previousBoundaries(N_PREV_BOUNDARIES, Boundary());
    data.reserve(4*nr);
    while (iZ != nz && iR != iR_crit+1) {
        //std::cout << "iz: " << iZ << " ir: " << iR << " iphi: " << iPhi << "\n";
        data.push_back(traceLine(nz, nr, nphi, zArray, rArray, phiArray, tPrevious, iZ, iR, iPhi, previousBoundaries));
        ns++;
        if (ns > nr*nz*nphi)
            return false;
    }

    //std::cout << "t_start: " << 0 <<  " t_end: " << tPrevious << "\n";

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
