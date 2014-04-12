#include "gaussiancoloumbattractionintegral.h"

#include <hermiteintegral.h>
#include <math/hermiteexpansioncoefficient.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

GaussianColoumbAttractionIntegral::GaussianColoumbAttractionIntegral(int angularMomentumMax) :
    m_hermiteExpansionCoefficient(angularMomentumMax + 1),
    m_hermiteIntegral(2 * angularMomentumMax + 1)
{

}

//GaussianColoumbAttractionIntegral::GaussianColoumbAttractionIntegral(const rowvec& corePositionA, const rowvec& corePositionB,
//                                                                     const rowvec& corePositionC,
//                                                                     double exponentA, double exponentB,
//                                                                     int angularMomentumMax) :
//    GaussianColoumbAttractionIntegral(angularMomentumMax)
//{
//    set(corePositionA, corePositionB, corePositionC, exponentA, exponentB);
//}

void GaussianColoumbAttractionIntegral::set(const rowvec& corePositionA, const rowvec& corePositionB,
                                            const rowvec& corePositionC,
                                            const GaussianPrimitiveOrbital& primitiveA,
                                            const GaussianPrimitiveOrbital& primitiveB) {
    double exponentA = primitiveA.exponent();
    double exponentB = primitiveB.exponent();
    double p = exponentA + exponentB;
    rowvec P = (exponentA * corePositionA + exponentB * corePositionB) / (exponentA + exponentB);
    rowvec PC = P - corePositionC;
    m_exponentSum = exponentA + exponentB;
    int t = primitiveA.xExponent() + primitiveB.xExponent();
    int u = primitiveA.yExponent() + primitiveB.yExponent();
    int v = primitiveA.zExponent() + primitiveB.zExponent();
    m_hermiteIntegral.set(p, PC, t, u, v);
    m_hermiteExpansionCoefficient.set(exponentA, exponentB, corePositionA, corePositionB,
                                      primitiveA.xExponent(), primitiveB.xExponent(),
                                      primitiveA.yExponent(), primitiveB.yExponent(),
                                      primitiveA.zExponent(), primitiveB.zExponent());
}

double GaussianColoumbAttractionIntegral::coloumbAttractionIntegral(const GaussianPrimitiveOrbital& primitiveA,
                                                                    const GaussianPrimitiveOrbital& primitiveB) {
    return coloumbAttractionIntegral(primitiveA.xExponent(), primitiveA.yExponent(), primitiveA.zExponent(),
                                     primitiveB.xExponent(), primitiveB.yExponent(), primitiveB.zExponent());
}

double GaussianColoumbAttractionIntegral::coloumbAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double result = 0;
    //    const cube &E_x = (*m_hermiteExpansionCoefficient)[0];
    //    const cube &E_y = (*m_hermiteExpansionCoefficient)[1];
    //        const cube &E_z = (*m_hermiteExpansionCoefficient)[2];
    const HermiteExpansionCoefficient &E = m_hermiteExpansionCoefficient;
    double p = m_exponentSum;
    const HermiteIntegral &R = m_hermiteIntegral;
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
