#ifndef __INJECTOR_H__
#define __INJECTOR_H__

#include  "PhysicValues.h"
#include <cmath>

struct Injector
{
    double rho;
    double phi;
    double z;

    double sigma;
    double r0;
    double theta;
    unsigned particles;

    double E;
    uint M;
    uint Z;

    Injector(double rho, double phi, double z, double sigma, double r0, double theta, unsigned particles, double E=0., uint M=1, uint Z=1) : rho(rho), phi(phi), z(z), sigma(sigma),
                                                                                                r0(r0), theta(theta), particles(particles), E(E), M(M), Z(Z)
    {}


    double getOmege(double B) const {
        return Z*PhysicValues::ee * B / (M*PhysicValues::MP*PhysicValues::C);
    }

    double getV() const {
        return sqrt(2.*E*PhysicValues::EV_TO_ERG / (M*PhysicValues::MP));
    }

    double  getVperp() const {
        return getV()*sin(theta);
    }

    double getVx() const {
        return -getVperp() * sin(phi);
    }

    double getVy() const {
        return getVperp() * cos(phi);
    }

    double getVz() const {
        return getV() * cos(theta);
    }

    double getDeltaX(double B) const {
        return getVx()/getOmege(B);
    }

    double getDeltaY(double B) const {
        return getVy() /getOmege(B);
    };

};


#endif
