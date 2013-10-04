#ifndef MULTIHYDROGEN_H
#define MULTIHYDROGEN_H

#include <src/electronsystems/electronsystem.h>

#include <armadillo>
#include <fstream>

using namespace arma;

class MultiHydrogen : public ElectronSystem
{
public:
    MultiHydrogen(mat nucleiPositions);
    virtual ~MultiHydrogen();

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
    double errorFunction(double arg);

    double additionalEnergyTerms();
private:
    uint m_nNuclei;
    uint m_nOrbitalsPerNuclei;
};

inline uint MultiHydrogen::nParticles() {
    return m_nNuclei;
}

inline uint MultiHydrogen::nOrbitals()
{
    return m_nNuclei * m_nOrbitalsPerNuclei;
}

#endif // MULTIHYDROGEN_H
