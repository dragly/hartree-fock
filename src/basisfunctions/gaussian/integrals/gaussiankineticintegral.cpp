#include "gaussiankineticintegral.h"

#include "math/vector3.h"
#include <math/hermiteexpansioncoefficient.h>
#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

GaussianKineticIntegral::GaussianKineticIntegral(int angularMomentumMax) :
    m_overlapIntegral(angularMomentumMax + 2)
{
}

void GaussianKineticIntegral::set(const Vector3& corePositionA, const Vector3& corePositionB,
                                  const GaussianPrimitiveOrbital& primitiveA, const GaussianPrimitiveOrbital& primitiveB)
{
    m_exponentB = primitiveB.exponent();
    m_exponentSum = primitiveA.exponent() + primitiveB.exponent();
    m_overlapIntegral.set(corePositionA, corePositionB, primitiveA, primitiveB, true);
}

double GaussianKineticIntegral::kineticIntegral(int dim, int iA, int iB) {
    double b = m_exponentB;
    double S_iA_iBnn = m_overlapIntegral.overlapIntegral(dim, iA, iB + 2);
    double S_iA_iB = m_overlapIntegral.overlapIntegral(dim, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= m_overlapIntegral.overlapIntegral(dim, iA, iB - 2);
    } else {
        S_iA_iBpp = 0;
    }
    return 4 * b * b * S_iA_iBnn - 2*b * (2*iB + 1) * S_iA_iB + iB * (iB - 1) * S_iA_iBpp;
}

double GaussianKineticIntegral::kineticIntegral(const GaussianPrimitiveOrbital& primitiveA,
                                                const GaussianPrimitiveOrbital& primitiveB) {
    return kineticIntegral(primitiveA.xExponent(), primitiveA.yExponent(), primitiveA.zExponent(),
                           primitiveB.xExponent(), primitiveB.yExponent(), primitiveB.zExponent());
}

double GaussianKineticIntegral::kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double T_iA_iB = kineticIntegral(0, iA, iB);
    double T_jA_jB = kineticIntegral(1, jA, jB);
    double T_kA_kB = kineticIntegral(2, kA, kB);

    double S_iA_iB = m_overlapIntegral.overlapIntegral(0, iA, iB);
    double S_jA_jB = m_overlapIntegral.overlapIntegral(1, jA, jB);
    double S_kA_kB = m_overlapIntegral.overlapIntegral(2, kA, kB);

    double result = T_iA_iB * S_jA_jB * S_kA_kB + S_iA_iB * T_jA_jB * S_kA_kB + S_iA_iB * S_jA_jB * T_kA_kB;
    result *= -0.5;
    return result; // TODO: Is there an error in the slides here?
}
