#ifndef __INJECTOR_H__
#define __INJECTOR_H__

struct Injector
{
    double rho;
    double phi;
    double z;

    double sigma;
    double r0;
    double theta;
    unsigned particles;

    Injector(double rho, double phi, double z, double sigma, double r0, double theta, unsigned particles) : rho(rho), phi(phi), z(z), sigma(sigma),
                                                                                                r0(r0), theta(theta), particles(particles)
    {}

};


#endif
