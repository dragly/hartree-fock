#ifndef MULTIHYDROGEN_H
#define MULTIHYDROGEN_H

#include <electronsystems/electronsystem.h>

#include <armadillo>
#include <fstream>

using namespace arma;

class MultiHydrogen : public ElectronSystem
{
public:
    MultiHydrogen(mat nucleiPositions);
    virtual ~MultiHydrogen();

    virtual double coupledIntegral(int p, int r, int q, int s);
    virtual double kineticIntegral(int p, int q);
    virtual double nuclearAttractionIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);

    virtual uint nBasisFunctions();
    virtual uint nParticles();

    mat R;

    vec alpha;

    // Constants to be precalculated
    const double powPi5over2 = pow(M_PI, 5./2.);
    double errorFunction(double arg);

    double additionalEnergyTerms();
private:
    uint m_nNuclei;
    uint m_nOrbitalsPerNuclei;

    // ElectronSystem interface
public:
    double uncoupledIntegral(int p, int q);
};

inline uint MultiHydrogen::nParticles() {
    return m_nNuclei;
}

inline uint MultiHydrogen::nBasisFunctions()
{
    return m_nNuclei * m_nOrbitalsPerNuclei;
}

#endif // MULTIHYDROGEN_H
