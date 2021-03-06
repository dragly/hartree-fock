#ifndef HYDROGENMOLECULE_H
#define HYDROGENMOLECULE_H

#include <electronsystems/electronsystem.h>

#include <armadillo>
#include <fstream>

using namespace arma;

class HydrogenMolecule : public ElectronSystem
{
public:
    HydrogenMolecule(double distance = 1.4);
    virtual ~HydrogenMolecule();

    virtual double coupledIntegral(int p, int q, int r, int s);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);

    virtual uint nBasisFunctions();
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

    // ElectronSystem interface
public:
    double uncoupledIntegral(int p, int q);
};

inline uint HydrogenMolecule::nParticles() {
    return 2;
}

inline uint HydrogenMolecule::nBasisFunctions()
{
    return nNuclei * nOrbitalsPerNuclei;
}

#endif // HYDROGENMOLECULE_H
