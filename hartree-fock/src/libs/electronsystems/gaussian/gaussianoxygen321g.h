#ifndef GAUSSIANOXYGEN321G_H
#define GAUSSIANOXYGEN321G_H

#include <electronsystems/electronsystem.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

class GaussianOxygen321G : public ElectronSystem
{
public:
    GaussianOxygen321G();

    // ElectronSystem interface
public:
    virtual double coupledIntegral(int p, int q, int r, int s);
    virtual double uncoupledIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);
    virtual uint nBasisFunctions();
    virtual uint nParticles();
    virtual double additionalEnergyTerms();

private:
    uint m_nParticles;
    uint m_nBasisFunctions;
    int m_angularMomentumMax;
    double m_coreCharge;

    vector<GaussianContractedOrbital> m_basisFunctions;
    mat m_corePositions;
    GaussianElectronInteractionIntegral electronInteractionIntegral;
};

#endif // GAUSSIANOXYGEN321G_H
