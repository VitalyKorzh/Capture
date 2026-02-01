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

void Counter::process(const Line &line, const Injector &injector, double gamma)
{
    double P1 = 0.;
    double P2 = 0.;
    bool found = false;
    double integral = 0.;
    for (uint is = 0; is < line.getNs(); is++)
    {
        const LineData &data = line.getData()[is];
        uint ir = data.iR;
        uint iz = data.iZ;
        uint iphi = data.iPhi;
        intersectionCell[getIndex(iz, ir, iphi)] = true;
        integral += data.s*ni[getIndex(iz, ir, 0)]*injector.sigma*reader.normaDensity;
        P2 = 1. - exp(-integral);

        if (P1 < gamma && P2 > gamma)
        {
            nCap[getIndex(iz, ir, iphi)]++;
            found = true;
            break;
        }

        P1 = P2;
    }

    if (!found)
        nFlyby++;

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

    for (Injector injector : injectors)
    {
        double rho = injector.rho;
        double phi = injector.phi;
        double r0 = injector.r0;

        if (r0 == 0) {
            Line line(rho, injector.theta, phi,  injector.z, nz, nr, nphi, zArray, rArray,  phiArray, reader.t_epsilon);
            for (uint it = 0; it < injector.particles; it++)
                if (line.getLineWork())
                    process(line, injector, distGamma(gen));
                else
                {
                    std::cerr << "не удалось построить линию!\n";
                    break;
                }
        }
        else
        {
            std::normal_distribution<double> dist_x(rho * cos(phi), r0);
            std::normal_distribution<double> dist_y(rho * sin(phi), r0);

            for (uint it = 0; it < injector.particles; it++)
            {
                double x = dist_x(gen);
                double y = dist_y(gen);

                double rho_new = sqrt(x*x+y*y);
                double phi_new = atan2(y, x);

                if (phi_new < 0) phi_new += 2 * M_PI;

                Line line(rho_new, injector.theta, phi_new,  injector.z, nz, nr, nphi, zArray, rArray,  phiArray, reader.t_epsilon);
                if (line.getLineWork())
                    process(line, injector, distGamma(gen));
                else {
                    std::cerr << "не удалось построить линию!\n";
                    break;
                }
            }
        }

        nParticles += injector.particles;
    }
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
    for (Injector injector : injectors)
    {
        os << "# \tinjector\n";
        os << "# \t\tparticles=" << injector.particles << "\n";
        os << "# \t\ttheta=" << injector.theta*180./M_PI << "\n";
        os << "# \t\tsigma=" << injector.sigma << "\n";
        os << "# \t\tr0=" << injector.r0 << "\n";
        os << "# \t\tposition\n";
        os << "# \t\t\trho=" << injector.rho << "\n";
        os << "# \t\t\tz=" << injector.z << "\n";
        os << "# \t\t\tphi=" << injector.phi*180./M_PI << "\n";
        os << "# \tinjector end\n";
    }

    os << "# count end\n";
    os << "#\n";
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
}

Counter::~Counter()
{
    TimeProfiler::print(os);
}