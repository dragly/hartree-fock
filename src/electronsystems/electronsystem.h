#ifndef ELECTRONSYSTEM_H
#define ELECTRONSYSTEM_H

#include <sys/types.h>

class ElectronSystem
{
public:
    ElectronSystem();

    virtual double coupledIntegral(int p, int r, int q, int s) = 0;
    virtual double uncoupledIntegral(int p, int q) = 0;
    virtual double overlapIntegral(int p, int q) = 0;

    virtual uint nBasisFunctions() = 0;
    virtual uint nParticles() = 0;
    virtual uint nParticlesUp() {
        return nParticles() - nParticlesDown();
    }
    virtual uint nParticlesDown() {
        if(m_nParticlesDownSet) {
            return m_nParticlesDown;
        } else {
            return nParticles() / 2;
        }
    }

    virtual double additionalEnergyTerms() = 0;

    void setNParticlesDown(int nParticles);

private:
    bool m_nParticlesDownSet;
    int m_nParticlesDown;
};

#endif // ELECTRONSYSTEM_H
