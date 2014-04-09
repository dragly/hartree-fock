#ifndef GAUSSIANSYSTEM_H
#define GAUSSIANSYSTEM_H

#include <electronsystems/electronsystem.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>

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
    double particleDensity(const mat &C, double x, double y, double z) const;
    void addCore(const GaussianCore& core);
protected:
    void setAngularMomentumMax(int angularMomentumMax);
private:
    uint m_nParticles;
    uint m_nBasisFunctions;
    int m_angularMomentumMax;
//    vector<double> m_coreCharges;

    vector<GaussianCore> m_cores;
    vector<GaussianContractedOrbital> m_basisFunctions;
//    vector<rowvec> m_corePositions;
    GaussianElectronInteractionIntegral electronInteractionIntegral;
    GaussianKineticIntegral kineticIntegral;
    GaussianColoumbAttractionIntegral coulombIntegral;
};

#endif // GAUSSIANSYSTEM_H
