#ifndef GAUSSIANSYSTEM_H
#define GAUSSIANSYSTEM_H

#include <electronsystems/electronsystem.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

class GaussianSystem : public ElectronSystem
{
public:
    GaussianSystem();

    // ElectronSystem interface
public:
    virtual double coupledIntegral(int p, int q, int r, int s);
    virtual double uncoupledIntegral(int p, int q);
    virtual double overlapIntegral(int p, int q);
    virtual uint nBasisFunctions();
    virtual uint nParticles();
    virtual double additionalEnergyTerms();

//    double particleDensity(double x, double y, double z);
    double particleDensity(const mat &C, double x, double y, double z);
protected:
    uint m_nParticles;
    uint m_nBasisFunctions;
    int m_angularMomentumMax;
    double m_coreCharge;

    vector<GaussianContractedOrbital> m_basisFunctions;
    mat m_corePositions;
    GaussianElectronInteractionIntegral electronInteractionIntegral;
};

#endif // GAUSSIANSYSTEM_H
