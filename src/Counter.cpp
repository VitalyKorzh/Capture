#include "Counter.h"
#include <cmath>
#include <random>
#include <iostream>
#include <forward_list>
#include "TimeProfiler.h"
#include "PhysicValues.h"
#include "Line.h"

void Counter::clearPrevious()
{
    nFlyby = 0;
    nParticles = 0;

    for (uint i = 0; i < nz*nr*nphi; i++)
    {
        intersectionCell[i] = false;
        intersectionCellCenter[i] = false;
        nCap[i] = 0;
        nCapCenter[i] = 0;
    }
}

Counter::Counter(std::istream &in, std::ostream &os) : reader(in), os(os), nz(reader.nz), nr(reader.nr), nphi(reader.nphi),
                                                       ni(reader.ni), zArray(reader.zArray), rArray(reader.rArray),
                                                       phiArray(reader.phiArray), injectors(reader.injectors), nParticles(0), intersectionCell(nr * nz * nphi, false),
                                                       intersectionCellCenter(nr*nz*nphi, false),
                                                       nCap(nz * nr * nphi, 0), nCapCenter(nz*nr*nphi, 0), nFlyby(0)
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
        uint index = getIndex(data.index);
        uint indexCenter = getIndex(data.indexCenter);
        bool isCenter = reader.count_centers && (data.indexCenter.iR < nr) && (data.indexCenter.iZ < nz) && (data.indexCenter.iPhi < nphi); //заполнять ли центры
        intersectionCell[index] = true;
        if (isCenter)
            intersectionCellCenter[indexCenter] = true;
        integral += data.s*ni[getIndex(iz, ir, 0)]*injector.sigma*reader.normaDensity;
        P2 = 1. - exp(-integral);

        if (P1 < gamma && P2 > gamma)
        {
            nCap[index]++;
            if  (isCenter)
                nCapCenter[indexCenter]++;

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

void Counter::countIonN(const uiarray &array, darray &result)
{
    result.resize(nr, 0.);
    uint nCaptureParticles = 0.;
    for (uint ir = 0; ir < nr; ir++)
    {
        double r1 = rArray[ir];
        double r2 = rArray[ir+1];
        uint Nions = 0;
        double volume = M_PI*(r2*r2-r1*r1);

        for (uint iphi = 0; iphi < nphi; iphi++) {
            for (uint iz = 0; iz < nz; iz++)
                Nions += array[getIndex(iz, ir, iphi)];
        }

        result[ir] = Nions / (volume * (nParticles-nFlyby));
        nCaptureParticles += Nions;
    }

    //std::cout << nParticles-nFlyby << " " << nCaptureParticles << "\n";

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


    //std::ofstream fout;
    //fout.open("points.txt");

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
                        injector.plusDirection, injector.getRhoLarmor(reader.Bcenter), injector.Rinjection, reader.count_centers);

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

            std::forward_list <Line> lines;

            uint N_LINES_GENERATE = injector.nLines;
            const uint part_in_line = injector.particles / N_LINES_GENERATE;
            bool lineError = false;
            for (uint l = 0; l < N_LINES_GENERATE; l++)
            {
                bool direction = injector.plusDirection;
                double x = r0 == 0. ? rho*cos(phi) : dist_x(gen);
                double y = r0 == 0. ? rho*sin(phi) : dist_y(gen);
                double theta_new = dtheta == 0. ? injector.theta : dist_theta(gen);
                if (theta_new > M_PI || theta_new < 0) {
                    nFlyby+=part_in_line;
                    extra_flyby+=part_in_line;
                    continue;
                }

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

                double rInj = injector.Rinjection < 0 ? rArray.back() : injector.Rinjection;

                if (rho_new >= std::min(rArray.back(), rInj)) // не расчитвать такой луч
                {
                    nFlyby+=part_in_line;
                    extra_flyby+=part_in_line;
                    continue;
                }


                lines.emplace_front(rho_new, theta_new, phi_new,  z_new, nz, nr, nphi, zArray, rArray,  phiArray, reader.t_epsilon, reader.t_epsilon_first, 
                direction, injector.getRhoLarmor(reader.Bcenter), injector.Rinjection, reader.count_centers);

                if (!lines.begin()->getLineWork())
                {
                    std::cerr << "не удалось построить линию, поменяйти сетку!\n";
                    lineError = true;
                    break;
                }

            }

            auto iline = lines.begin();
            
            for (uint it = 0; it < injector.particles-extra_flyby; it++)
            {

                Line &line = *iline;

                process(line, injector, distGamma(gen));

                iline++;
                if (iline == lines.end())
                    iline = lines.begin();
                // bool direction = injector.plusDirection;
                // double x = r0 == 0. ? rho*cos(phi) : dist_x(gen);
                // double y = r0 == 0. ? rho*sin(phi) : dist_y(gen);
                // double theta_new = dtheta == 0. ? injector.theta : dist_theta(gen);
                // if (theta_new > M_PI || theta_new < 0) {
                //     nFlyby++;
                //     extra_flyby++;
                //     continue;
                // }

                // double rho_new = x*cos0 + y*sin0;
                // double t_new = (y*cos0 - x*sin0)/(injector.sign*sin(theta_new));
                // double z_new = injector.z - t_new*cos(theta_new);
                // double phi_new = phi;
                // if (rho_new < 0)
                // {
                //     rho_new *= -1.;
                //     phi_new += M_PI;
                //     direction = !direction;
                // }

                // double rInj = injector.Rinjection < 0 ? rArray.back() : injector.Rinjection;

                // if (rho_new >= std::min(rArray.back(), rInj)) // не расчитвать такой луч
                // {
                //     nFlyby++;
                //     extra_flyby++;
                //     continue;
                // }

                // //fout << rho_new << " " << theta_new << " " << z_new << " " << phi_new << " " << (direction ? -1 : 1) << "\n";
                // Line line(rho_new, theta_new, phi_new,  z_new, nz, nr, nphi, zArray, rArray,  phiArray, reader.t_epsilon, reader.t_epsilon_first, 
                //                 direction, injector.getRhoLarmor(reader.Bcenter), injector.Rinjection, reader.count_centers);
                // if (line.getLineWork())
                //     process(line, injector, distGamma(gen));
                // else {
                //     std::cerr << "не удалось построить линию!\n";
                //     break;
                // }
            }
        }

        nParticles += injector.particles;
        index++;
    }
    //fout.close();
    if (extra_flyby != 0)
        std::cerr << "extra-flyby: " << extra_flyby << "\n"; 
    countIonN(nCap, nF);
    countIonN(nCapCenter, nFCenter);
}

void Counter::printStartInfo() const 
{
    os << "# precision=" << reader.precision << "\n";
    os << "# normaN=" << reader.normaDensity << "\n";
    os << "#\n";
    os << "# mesh\n";
    os << "# \tz-axis\n# \t\tn " << nz << "\n";
    for (const double & z0 : zArray)
        os << "  \t\t\t" << z0 << "\n";
    os << "# \tr-axis\n# \t\tn " << nr << "\n";
    for (const double & r0 : rArray)
        os << "  \t\t\t" << r0 << "\n";
    os << "# \tphi-axis\n# \t\tn " << nphi << "\n";
    for (const double & phi0 : phiArray)
        os << "  \t\t\t" << phi0 << "\n";
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
        os << "# \t\tlines=" << injector.nLines  << "\n";
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
            {
                uint index = getIndex(iz, ir, iphi);
                os << getnCap(index, nCap) << " " << getnCap(index, nCapCenter) << " " << (intersectionCell[index] ? 1 : 0) << " " 
                << (intersectionCellCenter[index] ? 1 : 0) << " ";
            }
            os << "\n";
        }
    }

    os << "#\n";
    os << "#n from r:\n";
    double integral_1 = 0;
    double integral_2 = 0;
    for (uint ir = 0; ir < nr; ir++){
        integral_1 += M_PI*(rArray[ir+1]*rArray[ir+1]-rArray[ir]*rArray[ir])*nF[ir];
        integral_2 += M_PI*(rArray[ir+1]*rArray[ir+1]-rArray[ir]*rArray[ir])*nFCenter[ir];
        os << nF[ir] << " " << nFCenter[ir] << "\n";
    }
    os << "# integral nF 2pirdr = " << integral_1 << "\n";
    os << "# integral nFCenter 2pirdr = " << integral_2 << "\n";
    os << "\n";
}

Counter::~Counter()
{
    TimeProfiler::print(os);
}