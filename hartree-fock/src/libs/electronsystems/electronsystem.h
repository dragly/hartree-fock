#ifndef ELECTRONSYSTEM_H
#define ELECTRONSYSTEM_H

#include <sys/types.h>

class ElectronSystem
{
public:
    ElectronSystem();

    virtual double electronInteractionIntegral(int p, int r, int q, int s) = 0;
    virtual double kineticIntegral(int p, int q) = 0;
    virtual double nuclearAttractionIntegral(int p, int q) = 0;
    virtual double overlapIntegral(int p, int q) = 0;

    virtual uint nOrbitals() = 0;
    virtual uint nParticles() = 0;
    virtual double additionalEnergyTerms() = 0;
};

#endif // ELECTRONSYSTEM_H
