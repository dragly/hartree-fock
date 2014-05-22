#include "gaussianoverlapintegral.h"

#include "math/vector3.h"
#include <math/hermiteexpansioncoefficient.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

GaussianOverlapIntegral::GaussianOverlapIntegral(int angularMomentumMax) :
    m_hermiteExpansionCoefficients(angularMomentumMax + 1)
{
}

void GaussianOverlapIntegral::set(Vector3 corePositionA, Vector3 corePositionB,
                                  const GaussianPrimitiveOrbital& primitiveA,
                                  const GaussianPrimitiveOrbital& primitiveB,
                                  bool expandForKinetic)
{
    int expansion = 0;
    if(expandForKinetic) {
        expansion = 2;
    }
    m_exponentSum = primitiveA.exponent() + primitiveB.exponent();
    m_hermiteExpansionCoefficients.set(primitiveA.exponent(), primitiveB.exponent(),
                                       corePositionA, corePositionB,
                                       primitiveA.xExponent(), primitiveB.xExponent() + expansion,
                                       primitiveA.yExponent(), primitiveB.yExponent() + expansion,
                                       primitiveA.zExponent(), primitiveB.zExponent() + expansion);
    m_exponentSum = primitiveA.exponent() + primitiveB.exponent();
}

double GaussianOverlapIntegral::overlapIntegral(const GaussianPrimitiveOrbital& primitiveA,
                                                const GaussianPrimitiveOrbital& primitiveB)
{
    return overlapIntegral(primitiveA.xExponent(), primitiveA.yExponent(), primitiveA.zExponent(),
                           primitiveB.xExponent(), primitiveB.yExponent(), primitiveB.zExponent());
}

double GaussianOverlapIntegral::overlapIntegral(int dim, int iA, int iB)
{
    double p = m_exponentSum;
    const cube &E_dim = m_hermiteExpansionCoefficients[dim];
    return E_dim(iA,iB,0) * sqrt(M_PI / p);
}

double GaussianOverlapIntegral::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB)
{
    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}
