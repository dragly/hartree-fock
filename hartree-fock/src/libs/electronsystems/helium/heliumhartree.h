#ifndef HELIUMHARTREE_H
#define HELIUMHARTREE_H

#include <electronsystems/electronsystem.h>

#include <armadillo>

using namespace arma;

class Helium : public ElectronSystem
{
public:
    Helium();

    rowvec alpha;

    // Constants to be precalculated
    const double powPi5over2 = pow(M_PI, 5./2.);
    double banana() const;
    void setBanana(double banana);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);

    // ElectronSystem interface
public:
    virtual double uncoupledIntegral(int p, int q);

    virtual double coupledIntegral(int p, int r, int q, int s);
    virtual double overlapIntegral(int p, int q);
    virtual double additionalEnergyTerms();

    virtual uint nOrbitals();
    virtual uint nParticles();
};

inline uint Helium::nParticles() {
    return 2;
}

#endif // HELIUMHARTREE_H
