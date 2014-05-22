#ifndef GAUSSIANTYPEOVERLAPINTEGRAL_H
#define GAUSSIANTYPEOVERLAPINTEGRAL_H

#include <armadillo>
#include "math/hermiteexpansioncoefficient.h"

class GaussianPrimitiveOrbital;
class Vector3;

using namespace arma;

class GaussianOverlapIntegral
{
public:
    explicit GaussianOverlapIntegral(int angularMomentumMax);

    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double overlapIntegral(int dim, int i, int j);
    double overlapIntegral(const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    void set(Vector3 corePositionA, Vector3 corePositionB, const GaussianPrimitiveOrbital &primitiveA,
             const GaussianPrimitiveOrbital &primitiveB, bool expandForKinetic = false);
protected:
    double m_exponentSum;

    HermiteExpansionCoefficient m_hermiteExpansionCoefficients;
};

#endif // GAUSSIANTYPEOVERLAPINTEGRAL_H
