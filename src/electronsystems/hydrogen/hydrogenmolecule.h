#ifndef HYDROGENMOLECULE_H
#define HYDROGENMOLECULE_H

#include <src/electronsystems/electronsystem.h>

#include <armadillo>
#include <fstream>

using namespace arma;

class HydrogenMolecule : public ElectronSystem
{
public:
    HydrogenMolecule(double distance = 1.4);
    virtual ~HydrogenMolecule();

    virtual double electronInteractionIntegral(int p, int r, int q, int s);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);

    virtual uint nOrbitals();
    virtual uint nParticles();

    mat R;

    rowvec alpha;

    // Constants to be precalculated
    const double powPi5over2 = pow(M_PI, 5./2.);
    double boysFunction(double arg);

    virtual double additionalEnergyTerms();
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
