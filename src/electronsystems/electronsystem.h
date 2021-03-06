#ifndef ELECTRONSYSTEM_H
#define ELECTRONSYSTEM_H

#include <sys/types.h>
#include <stdexcept>
#include <iostream>

class ElectronSystem
{
public:
    ElectronSystem();

    virtual double coupledIntegral(int p, int r, int q, int s) = 0;
    virtual double uncoupledIntegral(int p, int q) = 0;
    virtual double overlapIntegral(int p, int q) = 0;

    virtual uint nBasisFunctions() = 0;
    virtual uint nParticles() = 0;
    virtual uint nParticlesUp();
    virtual uint nParticlesDown();

    virtual double additionalEnergyTerms() = 0;

    void setNParticlesDown(uint nParticlesDown);

private:
    bool m_nParticlesDownSet;
    uint m_nParticlesDown;
};

#endif // ELECTRONSYSTEM_H
