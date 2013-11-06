#include "gaussiankineticintegral.h"

#include <math/hermiteexpansioncoefficient.h>
#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>

GaussianTypeKineticIntegral::GaussianTypeKineticIntegral(rowvec corePositionA, rowvec corePositionB,
                                                         double exponentA, double exponentB,
                                                         int angularMomentumMax) :
    GaussianTypeKineticIntegral(exponentA, exponentB,
                                new HermiteExpansionCoefficient(exponentA, exponentB, corePositionA, corePositionB, angularMomentumMax))
{
    m_isResponsibleForDeletingHermiteExpansionObject = true;
}

GaussianTypeKineticIntegral::GaussianTypeKineticIntegral(double exponentA,  double exponentB,
                                                         HermiteExpansionCoefficient* hermiteExpansionCoefficient) :
    m_exponentB(exponentB),
    m_exponentSum(exponentA + exponentB),
    m_hermiteExpansionCoefficient(hermiteExpansionCoefficient),
    m_isResponsibleForDeletingHermiteExpansionObject(false)
{
}

GaussianTypeKineticIntegral::~GaussianTypeKineticIntegral()
{
    if(m_isResponsibleForDeletingHermiteExpansionObject) {
        delete m_hermiteExpansionCoefficient;
    }
}

double GaussianTypeKineticIntegral::kineticIntegral(int dim, int iA, int iB) {
    GaussianTypeOverlapIntegral overlapIntegral(m_exponentSum, m_hermiteExpansionCoefficient);
    double b = m_exponentB;
    double S_iA_iBnn = overlapIntegral.overlapIntegral(dim, iA, iB + 2);
    double S_iA_iB = overlapIntegral.overlapIntegral(dim, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= overlapIntegral.overlapIntegral(dim, iA, iB - 2);
    } else {
        S_iA_iBpp = 0;
    }
    return 4 * b * b * S_iA_iBnn - 2*b * (2*iB + 1) * S_iA_iB + iB * (iB + 1) * S_iA_iBpp;
}

double GaussianTypeKineticIntegral::kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    GaussianTypeOverlapIntegral overlapIntegral(m_exponentSum, m_hermiteExpansionCoefficient);
    double T_iA_iB = kineticIntegral(0, iA, iB);
    double T_jA_jB = kineticIntegral(1, jA, jB);
    double T_kA_kB = kineticIntegral(2, kA, kB);

    double S_iA_iB = overlapIntegral.overlapIntegral(0, iA, iB);
    double S_jA_jB = overlapIntegral.overlapIntegral(1, jA, jB);
    double S_kA_kB = overlapIntegral.overlapIntegral(2, kA, kB);

    double result = T_iA_iB * S_jA_jB * S_kA_kB + S_iA_iB * T_jA_jB * S_kA_kB + S_iA_iB * S_jA_jB * T_kA_kB;
    result *= -0.5;
    return result; // TODO: Is there an error in the slides here?
}
