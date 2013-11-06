#include "gaussiancoloumbattractionintegral.h"

#include <hermiteintegral.h>
#include <math/hermiteexpansioncoefficient.h>

GaussianTypeColoumbAttractionIntegral::GaussianTypeColoumbAttractionIntegral(rowvec corePositionA, rowvec corePositionB, rowvec corePositionC,
                                                                             double exponentA, double exponentB,
                                                                             int angularMomentumMax) :
    GaussianTypeColoumbAttractionIntegral(exponentA + exponentB,0,0)
{
    double p = exponentA + exponentB;
    rowvec P = (exponentA * corePositionA + exponentB * corePositionB) / (exponentA + exponentB);
    rowvec PC = P - corePositionC;
    m_hermiteExpansionCoefficient = new HermiteExpansionCoefficient(exponentA, exponentB, corePositionA, corePositionB, angularMomentumMax);
    m_hermiteIntegral = new HermiteIntegral(p, PC, 2 * angularMomentumMax);
    m_isResponsibleForDeletingHermiteObjects = true;
}

GaussianTypeColoumbAttractionIntegral::GaussianTypeColoumbAttractionIntegral(double exponentSum,
                                                                             HermiteExpansionCoefficient* hermiteExpansionCoefficient,
                                                                             HermiteIntegral* hermiteIntegral) :
    m_hermiteExpansionCoefficient(hermiteExpansionCoefficient),
    m_hermiteIntegral(hermiteIntegral),
    m_exponentSum(exponentSum),
    m_isResponsibleForDeletingHermiteObjects(false)
{
}

GaussianTypeColoumbAttractionIntegral::~GaussianTypeColoumbAttractionIntegral()
{
    if(m_isResponsibleForDeletingHermiteObjects) {
        delete m_hermiteExpansionCoefficient;
        delete m_hermiteIntegral;
    }
}

double GaussianTypeColoumbAttractionIntegral::coloumbAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double result = 0;
//    const cube &E_x = (*m_hermiteExpansionCoefficient)[0];
//    const cube &E_y = (*m_hermiteExpansionCoefficient)[1];
//        const cube &E_z = (*m_hermiteExpansionCoefficient)[2];
    const HermiteExpansionCoefficient E = (*m_hermiteExpansionCoefficient);
    double p = m_exponentSum;
    const HermiteIntegral &R = (*m_hermiteIntegral);
    int tMax = iA + iB;
    int uMax = jA + jB;
    int vMax = kA + kB;
    for(int t = 0; t < tMax + 1; t++) {
        for(int u = 0; u < uMax + 1; u++) {
            for(int v = 0; v < vMax + 1; v++) {
                result += E(iA, jA, kA, iB, jB, kB, t, u, v) * R(0,t,u,v);
//                result += E_x(iA, iB, t) * E_y(jA, jB, u) * E_z(kA, kB, v) * R(0,t,u,v);
            }
        }
    }
    result *= 2*M_PI / p;
    return result;
}
