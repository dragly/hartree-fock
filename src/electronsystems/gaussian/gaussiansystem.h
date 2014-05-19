#ifndef GAUSSIANSYSTEM_H
#define GAUSSIANSYSTEM_H

#include <electronsystems/electronsystem.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>

class Vector3;

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
    double corePotential(const Vector3 &position);
    double electronDensity(const mat &C, const Vector3 &position) const;
    void addCore(const GaussianCore& core);
    double orbitalDensity(uint orbital, const mat &C, double x, double y, double z) const;
    double electronPotential(const mat &C, const Vector3 position);
    double electronPotential(uint p, uint q, const Vector3 &position);
    double electrostaticPotential(const Vector3 &position);
    double electrostaticPotential(const mat &C, const Vector3 &position);
protected:
    void setAngularMomentumMax(int angularMomentumMax);
private:
    uint m_nParticles;
    uint m_nBasisFunctions;
    int m_angularMomentumMax;
//    vector<double> m_coreCharges;

    vector<GaussianCore> m_cores;
    vector<GaussianContractedOrbital> m_basisFunctions;
//    vector<Vector3> m_corePositions;
    GaussianElectronInteractionIntegral electronInteractionIntegral;
    GaussianKineticIntegral kineticIntegral;
    GaussianColoumbAttractionIntegral coulombIntegral;
};

#endif // GAUSSIANSYSTEM_H
