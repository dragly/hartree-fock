#ifndef HYDROGENMOLECULE_H
#define HYDROGENMOLECULE_H

#include <src/basisfunctions/basisfunction.h>

#include <armadillo>

using namespace arma;

class HydrogenMolecule : public BasisFunction
{
public:
    HydrogenMolecule();

    virtual double electronInteractionIntegral(int p, int r, int q, int s);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);

    virtual uint nOrbitals();

    mat R;

    vec alpha;

    // Constants to be precalculated
    const double powPi5over2 = pow(M_PI, 5./2.);
    double errorFunction(double arg);
};

#endif // HYDROGENMOLECULE_H
