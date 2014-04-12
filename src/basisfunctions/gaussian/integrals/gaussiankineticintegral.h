#ifndef GAUSSIANTYPEKINETICINTEGRAL_H
#define GAUSSIANTYPEKINETICINTEGRAL_H

#include <armadillo>
#include <math/hermiteexpansioncoefficient.h>

using namespace arma;

class HermiteExpansionCoefficient;
class GaussianPrimitiveOrbital;

class GaussianKineticIntegral
{
public:
    explicit GaussianKineticIntegral(int angularMomentumMax);
    void set(const rowvec &corePositionA, const rowvec &corePositionB,
             const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    double kineticIntegral(int dim, int iA, int iB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
protected:
    double m_exponentB;
    double m_exponentSum;

    HermiteExpansionCoefficient m_hermiteExpansionCoefficient;
};

#endif // GAUSSIANTYPEKINETICINTEGRAL_H
