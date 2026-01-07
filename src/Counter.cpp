#include "Counter.h"
#include <cmath>
#include <random>
#include <iostream>

#include "TimeProfiler.h"
#include "PhysicValues.h"

void Counter::clearPrevious()
{
    nFlyby = 0;
    for (uint & n0 : nCap)
        n0 = 0;
}

Counter::Counter(std::istream &in, std::ostream &os) : reader(in), os(os), nz(reader.nz), nr(reader.nr),
                                                        ni(reader.ni), zArray(reader.zArray), rArray(reader.rArray),
                                                        nParticles(reader.nParticles), sigma(reader.sigma), theta(reader.theta), 
                                                        sArray(reader.sArray), ns(reader.ns),
                                                        position(reader.position), nCap(nz*nr) 
{
    os.precision(reader.precision);
    os << std::scientific;
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

    for (uint it = 0; it < nParticles; it++)
    {
        double P1 = 0;
        double P2 = 0;
        bool found = false;
        double gamma = distGamma(gen);
        double integral = 0.;
        for (uint is = 0; is < ns; is++)
        {
            uint iz = reader.index[is].first;
            uint ir = reader.index[is].second;
            integral += sArray[is]*ni[iz]*sigma*reader.normaDensity;
            P2 = 1. - exp(-integral);

            if ( (P1 < gamma && P2 > gamma))
            {
                nCap[iz*nr+ir]++;
                found = true;
                break;
            }

            P1 = P2;
        }

        if (!found)
            nFlyby++;
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
    os << "# \tni\n";
    for (const double & ni0 : ni)
        os << "# \t\t" << ni0 << "\n";

    os << "#\n";

    os << "# count\n";
    os << "# \tparticles=" << nParticles << "\n";
    os << "# \tsigma=" << sigma << "\n";
    os << "# \ttheta=" << theta << "\n";
    os << "# \tposition\n";
    os << "# \t\tz " << position.first << "\n# \t\tr " << position.second << "\n";
    os << "#\n";
}

void Counter::printResult() const
{
    os << "# result:\n";
    os << "# " << "nFlyby=" << getnFlyby()*100. << "%" << "\n";
    os << "#\n";
    for (uint iz = 0; iz < nz; iz++)
    {
        for (uint ir = 0; ir < nr; ir++)
            os << getnCap(iz*nr+ir) << " ";
        os << "\n";
    }
}

Counter::~Counter()
{
    TimeProfiler::print(os);
}