#ifndef HYDROGENMOLECULE_H
#define HYDROGENMOLECULE_H

#include <src/basisfunctions/basisfunction.h>

#include <armadillo>

using namespace arma;

class HydrogenMolecule : public BasisFunction
{
public:
    HydrogenMolecule(double distance = 1.4);

    virtual double electronInteractionIntegral(int p, int r, int q, int s);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);

    virtual uint nOrbitals();
    virtual uint nParticles();

    mat R;

    vec alpha;

    // Constants to be precalculated
    const double powPi5over2 = pow(M_PI, 5./2.);
    double errorFunction(double arg);

    double nuclearRepulsion();
private:
    int nNuclei = 2;
    int nOrbitalsPerNuclei = 4;
};

inline uint HydrogenMolecule::nParticles() {
    return 2;
}

inline uint HydrogenMolecule::nOrbitals()
{
    return nNuclei * nOrbitalsPerNuclei;
}

#endif // HYDROGENMOLECULE_H
