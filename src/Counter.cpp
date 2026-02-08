#include "Counter.h"
#include <cmath>
#include <random>
#include <iostream>

#include "TimeProfiler.h"
#include "PhysicValues.h"
#include "Line.h"

void Counter::clearPrevious()
{
    nFlyby = 0;
    nParticles = 0;
    for (uint i = 0; i < nCap.size(); i++)
        nCap[i] = 0;

    for (uint i = 0; i < intersectionCell.size(); i++)
        intersectionCell[i] = false;
}

Counter::Counter(std::istream &in, std::ostream &os) : reader(in), os(os), nz(reader.nz), nr(reader.nr), nphi(reader.nphi),
                                                       ni(reader.ni), zArray(reader.zArray), rArray(reader.rArray),
                                                       phiArray(reader.phiArray), injectors(reader.injectors), nParticles(0), intersectionCell(nr * nz * nphi, false),
                                                       nCap(nz * nr * nphi), nFlyby(0)
{
    os.precision(reader.precision);
    os << std::scientific;
}

void Counter::process(const Line &line, const Injector &injector, double gamma, darray *sArray)
{
    double P1 = 0.;
    double P2 = 0.;
    bool found = false;
    double integral = 0.;
    uint ns = line.getNs();
    for (uint is = 0; is < ns; is++)
    {
        const LineData &data = line.getData()[is];
        uint ir = data.index.iR;
        uint iz = data.index.iZ;
        uint iphi = data.index.iPhi;
        intersectionCell[getIndex(iz, ir, iphi)] = true;
        integral += data.s*ni[getIndex(iz, ir, 0)]*injector.sigma*reader.normaDensity;
        P2 = 1. - exp(-integral);

        if (P1 < gamma && P2 > gamma)
        {
            nCap[getIndex(iz, ir, iphi)]++;

            if (sArray != nullptr)
                (*sArray)[1+2*ns+is]++;

            found = true;
            break;
        }

        P1 = P2;
    }

    if (!found)
        nFlyby++;

}

darray Counter::fillSArray(const Line &line) const
{
    uint ns = line.getNs();

    darray sArray(ns*3+1, 0.);

    for (uint i = 0; i < ns; i++)
    {
        sArray[1+i] = sArray[i] + line.getData()[i].s;
        sArray[1+ns+i] = ni[getIndex(line.getData()[i].index.iZ, line.getData()[i].index.iR, 0)];
    }

    return sArray;
}

void Counter::printSArray(const darray &sArray, uint i) const
{
    uint ns = (sArray.size()-1) / 3;
    os << "# result injector " << i << "\n";
    os << "# ns=" << ns << "\n";
    double full = 0;
    for (uint i = 0; i < ns; i++)
        full += sArray[i+1+2*ns];

    for (uint i = 0; i < ns; i++)
    {
        os << sArray[i] << 
        " " << sArray[i+1]
        << " " << sArray[i+1+ns] <<
        " " << sArray[i+1+2*ns] / full << "\n";
    }
    os << "#" << std::endl;
}

double Counter::cellVolume(uint iz, uint ir, uint iphi) const
{
    double r1 = rArray[ir];
    double r2 = rArray[ir+1];

    double z1 = zArray[iz];
    double z2 = zArray[iz+1];

    double phi1 = phiArray[iphi];
    double phi2 = phiArray[iphi+1];

    return (phi2-phi1)*(r2-r1)*(z2-z1)*(r1+r2)/2.;
}

void Counter::countIonN()
{
    nF.resize(nr, 0.);

    for (uint ir = 0; ir < nr; ir++)
    {
        double r1 = rArray[ir];
        double r2 = rArray[ir+1];
        uint Nions = 0;
        double volume = M_PI*(zArray.back()-zArray.front())*(r2*r2-r1*r1);

        for (uint iphi = 0; iphi < nphi; iphi++) {
            for (uint iz = 0; iz < nz; iz++)
                Nions += nCap[getIndex(iz, ir, iphi)];
        }

        nF[ir] = Nions / (volume * (nParticles-nFlyby));
    }

}

void Counter::count()
{
    TimeProfiler t_cout("time count full");
    if (!reader.work)
        return;
    clearPrevious();
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution <> distGamma(0., 1.);

    uint index = 0;

    uint extra_flyby = 0;
    for (Injector injector : injectors)
    {
        double rho = injector.rho;
        double phi = injector.phi;
        double r0 = injector.r0;
        double dtheta = injector.deltaTheta;

        darray sArray;

        if (r0 == 0. && dtheta == 0.) {
            Line line(rho, injector.theta, phi,  injector.z, nz, nr, nphi, zArray, rArray,  phiArray, reader.t_epsilon, reader.t_epsilon_first, 
                        injector.plusDirection, injector.getRhoLarmor(reader.Bcenter), injector.Rinjection);

            sArray = fillSArray(line);

            for (uint it = 0; it < injector.particles; it++)
                if (line.getLineWork())
                    process(line, injector, distGamma(gen), &sArray);
                else
                {
                    std::cerr << "не удалось построить линию!\n";
                    break;
                }

            printSArray(sArray, index);
        }
        else
        {
            std::normal_distribution<double> dist_x(rho * cos(phi), r0);
            std::normal_distribution<double> dist_y(rho * sin(phi), r0);
            std::normal_distribution<double> dist_theta(injector.theta, dtheta);

            double cos0 = cos(phi);
            double sin0 = sin(phi);

            
            for (uint it = 0; it < injector.particles; it++)
            {
                bool direction = injector.plusDirection;
                double x = r0 == 0. ? rho*cos(phi) : dist_x(gen);
                double y = r0 == 0. ? rho*sin(phi) : dist_y(gen);
                double theta_new = dtheta == 0. ? injector.theta : dist_theta(gen);
                double rho_new = x*cos0 + y*sin0;
                double t_new = (y*cos0 - x*sin0)/(injector.sign*sin(theta_new));
                double z_new = injector.z - t_new*cos(theta_new);
                double phi_new = phi;
                if (rho_new < 0)
                {
                    rho_new *= -1.;
                    phi_new += M_PI;
                    direction = !direction;
                }

                if (rho_new >= rArray.back()) // не расчитвать такой луч
                {
                    nFlyby++;
                    extra_flyby++;
                    continue;
                }


                Line line(rho_new, theta_new, phi_new,  z_new, nz, nr, nphi, zArray, rArray,  phiArray, reader.t_epsilon, reader.t_epsilon_first, 
                                direction, injector.getRhoLarmor(reader.Bcenter), injector.Rinjection);
                if (line.getLineWork())
                    process(line, injector, distGamma(gen));
                else {
                    std::cerr << "не удалось построить линию!\n";
                    break;
                }
            }
        }

        nParticles += injector.particles;
        index++;
    }

    if (extra_flyby != 0)
        std::cerr << "extra-flyby: " << extra_flyby << "\n"; 
    countIonN();
}

void Counter::printStartInfo() const 
{
    os << "# precision=" << reader.precision << "\n";
    os << "# normaN=" << reader.normaDensity << "\n";
    os << "#\n";
    os << "# mesh\n";
    os << "# \tz-axis\n# \t\tn " << nz << "\n";
    for (const double & z0 : zArray)
        os << "# \t\t\t" << z0 << "\n";
    os << "# \tr-axis\n# \t\tn " << nr << "\n";
    for (const double & r0 : rArray)
        os << "# \t\t\t" << r0 << "\n";
    os << "# \tphi-axis\n# \t\tn " << nphi << "\n";
    for (const double & phi0 : phiArray)
        os << "# \t\t\t" << phi0 << "\n";
    os << "# \tni\n";
    for (uint iz = 0; iz < nz; iz++)
    {
        for (uint ir = 0; ir < nr; ir++)
            os << ni[ir+nr*iz] << " ";
        os << "\n";
    }
    os << "#\n";

    os << "# count\n";
    os << "# \tt-epsilon=" << reader.t_epsilon << "\n";
    os << "# \tt-epsilon-first=" << reader.t_epsilon_first << "\n";
    os << "# \tBcenter=" << reader.Bcenter << "\n";
    for (Injector injector : injectors)
    {
        os << "# \tinjector\n";
        os << "# \t\tparticles=" << injector.particles << "\n";
        os << "# \t\ttheta=" << injector.theta*180./M_PI << "\n";
        os << "# \t\tsigma=" << injector.sigma << "\n";
        os << "# \t\tr0=" << injector.r0 << "\n";
        os << "# \t\tdtheta=" << injector.deltaTheta*180./M_PI << "\n";
        os << "# \t\tr-injection=" << (injector.Rinjection < 0 ? rArray.back() : injector.Rinjection) << "\n";
        os << "# \t\tposition\n";
        os << "# \t\t\trho=" << injector.rho << "\n";
        os << "# \t\t\tz=" << injector.z << "\n";
        os << "# \t\t\tphi=" << injector.phi*180./M_PI << "\n";
        os << "# \t\tparticle\n";
        os << "# \t\t\tE=" << injector.E << "\n";
        os << "# \t\t\tZ=" << injector.Z << "\n";
        os << "# \t\t\tM=" << injector.M << "\n";
        if (injector.plusDirection)
            os << "# \t\t[plus-direction]\n";
        else
            os << "# \t\t[minus-direction]\n";
        os << "# \tinjector end\n#\n";
    }

    os << "# count end\n";
    os << "#" << std::endl;
}

void Counter::printResult() const
{
    os << "# result:\n";
    os << "#\n";

    os << "# nFlyby=" << getnFlyby()*100. << " %\n";
    os << "# nCapture=" << 100. * (1. - getnFlyby()) << " %\n";

    for  (uint iphi = 0; iphi < nphi; iphi++)
    {
        os << "iphi=" << iphi << "\n";
        for (uint iz = 0; iz < nz; iz++)
        {
            for (uint ir=0; ir < nr; ir++)
                os << getnCap(getIndex(iz, ir, iphi)) << " " << (intersectionCell[getIndex(iz, ir, iphi)] ? 1 : 0) << " ";
            os << "\n";
        }
    }

    os << "#\n";
    os << "#n from r:\n";
    for (uint ir = 0; ir < nr; ir++)
        os << nF[ir] << "\n";
    os << "\n";
}

Counter::~Counter()
{
    TimeProfiler::print(os);
}