#ifndef GAUSSIANTYPEKINETICINTEGRAL_H
#define GAUSSIANTYPEKINETICINTEGRAL_H

#include <armadillo>
//#include <math/hermiteexpansioncoefficient.h>
#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>

class Vector3;

using namespace arma;

class GaussianPrimitiveOrbital;

class GaussianKineticIntegral
{
public:
    explicit GaussianKineticIntegral(int angularMomentumMax);
    void set(const Vector3 &corePositionA, const Vector3 &corePositionB,
             const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    double kineticIntegral(int dim, int iA, int iB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double kineticIntegral(const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
protected:
    double m_exponentB;
    double m_exponentSum;

//    HermiteExpansionCoefficient m_hermiteExpansionCoefficient;
    GaussianOverlapIntegral m_overlapIntegral;
};

#endif // GAUSSIANTYPEKINETICINTEGRAL_H
