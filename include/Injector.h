#ifndef __INJECTOR_H__
#define __INJECTOR_H__

#include "Physics/PhysValues.h"
#include <cmath>

struct Injector
{
    double rho;
    double phi;
    double z;

    double sigma;
    double r0;
    double theta;
    double deltaTheta;
    double Rinjection;
    unsigned particles;

    double E;
    uint M;
    uint Z;
    bool plusDirection;
    int sign;


    uint nLines;


    Injector(double rho, double phi, double z, double sigma, double r0, double theta, double deltaTheta, double Rinjection, unsigned particles, double E=0., uint M=1, uint Z=1, bool plusDirection=true, uint nLines=1) : rho(rho), phi(phi), z(z), sigma(sigma),
                                                                                                r0(r0), theta(theta), deltaTheta(deltaTheta), Rinjection(Rinjection), 
                                                                                                particles(particles), E(E), M(M), Z(Z), plusDirection(plusDirection),
                                                                                                sign(plusDirection ? 1 : -1), nLines(nLines)
    {}


    double getOmege(double B) const {
        return Z*Physics::PhysicValues::e * B / (M*Physics::PhysicValues::MP*Physics::PhysicValues::C);
    }

    double getV() const {
        return sqrt(2.*E*Physics::PhysicValues::EV_TO_ERG / (M*Physics::PhysicValues::MP));
    }

    double  getVperp() const {
        return getV()*sin(theta);
    }

    double getVx() const {
        return -getVperp() * sin(phi) * sign;
    }

    double getVy() const {
        return getVperp() * cos(phi) * sign;
    }

    double getVz() const {
        return getV() * cos(theta);
    }

    double getRhoLarmor(double B) const {
        return getVperp() / getOmege(B);
    }

    double getDeltaX(double B) const {
        return getVy()/getOmege(B);
    }

    double getDeltaY(double B) const {
        return -getVx() /getOmege(B);
    };

};


#endif
