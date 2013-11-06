#include "gaussianoverlapintegral.h"

#include <math/hermiteexpansioncoefficient.h>

GaussianTypeOverlapIntegral::GaussianTypeOverlapIntegral(rowvec corePositionA, rowvec corePositionB,
                                                         double exponentA, double exponentB,
                                                         int angularMomentumMax) :
    GaussianTypeOverlapIntegral(exponentA + exponentB,
                                new HermiteExpansionCoefficient(exponentA, exponentB, corePositionA, corePositionB, angularMomentumMax))
{
    m_isResponsibleForDeletingHermiteExpansionObject = true;
}

GaussianTypeOverlapIntegral::GaussianTypeOverlapIntegral(double exponentSum,
                                                         HermiteExpansionCoefficient* hermiteExpansionCoefficient) :
    m_exponentSum(exponentSum),
    m_hermiteExpansionCoefficient(hermiteExpansionCoefficient),
    m_isResponsibleForDeletingHermiteExpansionObject(false)
{
}

GaussianTypeOverlapIntegral::~GaussianTypeOverlapIntegral()
{
    if(m_isResponsibleForDeletingHermiteExpansionObject) {
        delete m_hermiteExpansionCoefficient;
    }
}

double GaussianTypeOverlapIntegral::overlapIntegral(int dim, int iA, int iB)
{
    double p = m_exponentSum;
    const cube &E_dim = (*m_hermiteExpansionCoefficient)[dim];
    return E_dim(iA,iB,0) * sqrt(M_PI / p);
}

double GaussianTypeOverlapIntegral::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}
