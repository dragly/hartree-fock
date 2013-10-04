#ifndef HELIUMHARTREE_H
#define HELIUMHARTREE_H

#include <src/electronsystems/electronsystem.h>

#include <armadillo>

using namespace arma;

class Helium : public ElectronSystem
{
public:
    Helium();

    virtual double electronInteractionIntegral(int p, int r, int q, int s);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);
    virtual double additionalEnergyTerms();

    virtual uint nOrbitals();
    virtual uint nParticles();

    rowvec alpha;

    // Constants to be precalculated
    const double powPi5over2 = pow(M_PI, 5./2.);
    double banana() const;
    void setBanana(double banana);
};

inline uint Helium::nParticles() {
    return 2;
}

#endif // HELIUMHARTREE_H
