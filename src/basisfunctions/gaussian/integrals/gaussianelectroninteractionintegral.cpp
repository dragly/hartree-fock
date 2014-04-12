#include "gaussianelectroninteractionintegral.h"

#include <math/hermiteexpansioncoefficient.h>
#include <hermiteintegral.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

GaussianElectronInteractionIntegral::GaussianElectronInteractionIntegral(int singleAngularMomentumMax) :
    m_hermiteIntegral(4 * singleAngularMomentumMax+1),
    m_hermiteExpansionCoefficientAB(singleAngularMomentumMax+1),
    m_hermiteExpansionCoefficientCD(singleAngularMomentumMax+1)
{
    reset(singleAngularMomentumMax);
}

void GaussianElectronInteractionIntegral::reset(int singleAngularMomentumMax)
{
    m_hermiteIntegral.reset(4*singleAngularMomentumMax+1);
    m_hermiteExpansionCoefficientAB.reset(singleAngularMomentumMax+1);
    m_hermiteExpansionCoefficientCD.reset(singleAngularMomentumMax+1);
}

void GaussianElectronInteractionIntegral::set(const rowvec &corePositionA, const rowvec &corePositionB, const rowvec &corePositionC, const rowvec &corePositionD, const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB, const GaussianPrimitiveOrbital &primitiveC, const GaussianPrimitiveOrbital &primitiveD)
{
    setAB(corePositionA, corePositionB, primitiveA, primitiveB);
    setCD(corePositionC, corePositionD, primitiveC, primitiveD);
}

void GaussianElectronInteractionIntegral::setAB(const rowvec &corePositionA, const rowvec &corePositionB,
                                                const GaussianPrimitiveOrbital& primitiveA, const GaussianPrimitiveOrbital& primitiveB)
{
    m_primitiveA = &primitiveA;
    m_primitiveB = &primitiveB;
    double a = primitiveA.exponent();
    double b = primitiveB.exponent();
    double p = a + b;
    m_exponentP = p;
    rowvec P = (a * corePositionA + b * corePositionB) / (a + b);
    m_centerOfMassP = P;
    m_hermiteExpansionCoefficientAB.set(a, b, corePositionA, corePositionB,
                                        primitiveA.xExponent(), primitiveB.xExponent(),
                                        primitiveA.yExponent(), primitiveB.yExponent(),
                                        primitiveA.zExponent(), primitiveB.zExponent());
}

void GaussianElectronInteractionIntegral::setCD(const rowvec &corePositionC, const rowvec &corePositionD,
                                                const GaussianPrimitiveOrbital& primitiveC, const GaussianPrimitiveOrbital& primitiveD)
{
    m_primitiveC = &primitiveC;
    m_primitiveD = &primitiveD;
    double c = primitiveC.exponent();
    double d = primitiveD.exponent();
    double p = m_exponentP;
    double q = c + d;
    m_exponentQ = q;
    rowvec P = m_centerOfMassP;
    rowvec Q = (c * corePositionC + d * corePositionD) / (c + d);
    m_centerOfMassQ = Q;
    rowvec PQ = P - Q;
    double alpha = p*q/(p+q);
    int tPlusTau = m_primitiveA->xExponent() + m_primitiveB->xExponent() + m_primitiveC->xExponent() + m_primitiveD->xExponent();
    int uPlusNu = m_primitiveA->yExponent() + m_primitiveB->yExponent() + m_primitiveC->yExponent() + m_primitiveD->yExponent();
    int vPlusPhi = m_primitiveA->zExponent() + m_primitiveB->zExponent() + m_primitiveC->zExponent() + m_primitiveD->zExponent();
    m_hermiteIntegral.set(alpha, PQ, tPlusTau, uPlusNu, vPlusPhi);
    m_hermiteExpansionCoefficientCD.set(c, d, corePositionC, corePositionD,
                                        primitiveC.xExponent(), primitiveD.xExponent(),
                                        primitiveC.yExponent(), primitiveD.yExponent(),
                                        primitiveC.zExponent(), primitiveD.zExponent());
}

GaussianElectronInteractionIntegral::~GaussianElectronInteractionIntegral()
{
}

double GaussianElectronInteractionIntegral::electronInteractionIntegral(const GaussianPrimitiveOrbital& primitiveA,
                                                                        const GaussianPrimitiveOrbital& primitiveB,
                                                                        const GaussianPrimitiveOrbital& primitiveC,
                                                                        const GaussianPrimitiveOrbital& primitiveD) {
    return electronInteractionIntegral(primitiveA.xExponent(), primitiveA.yExponent(), primitiveA.zExponent(),
                                       primitiveB.xExponent(), primitiveB.yExponent(), primitiveB.zExponent(),
                                       primitiveC.xExponent(), primitiveC.yExponent(), primitiveC.zExponent(),
                                       primitiveD.xExponent(), primitiveD.yExponent(), primitiveD.zExponent());
}

double GaussianElectronInteractionIntegral::electronInteractionIntegral(int iA, int jA, int kA,
                                                                        int iB, int jB, int kB,
                                                                        int iC, int jC, int kC,
                                                                        int iD, int jD, int kD) {
    double result = 0;
    double p = m_exponentP;
    double q = m_exponentQ;

    const HermiteIntegral &R = m_hermiteIntegral;
    const HermiteExpansionCoefficient &Eab = m_hermiteExpansionCoefficientAB;
    const HermiteExpansionCoefficient &Ecd = m_hermiteExpansionCoefficientCD;
    int tMax = iA + iB;
    int uMax = jA + jB;
    int vMax = kA + kB;
    int tauMax = iC + iD;
    int nuMax = jC + jD;
    int phiMax = kC + kD;
    for (int t = 0; t < tMax + 1; ++t) {
        for (int u = 0; u < uMax + 1; ++u) {
            for (int v = 0; v < vMax + 1; ++v) {
                for (int tau = 0; tau < tauMax + 1; ++tau) {
                    for (int nu = 0; nu < nuMax + 1; ++nu) {
                        for (int phi = 0; phi < phiMax + 1; ++phi) {
                            double product = 1;
                            product *= Eab(iA, jA, kA, iB, jB, kB, t, u, v);
                            product *= Ecd(iC, jC, kC, iD, jD, kD, tau, nu, phi);
//                            product *= pow(-1, tau + nu + phi);
                            product *= (1 - 2*((tau + nu + phi) % 2));
                            product *= R(0, t + tau, u + nu, v + phi);
                            result += product;
                        }
                    }
                }
            }
        }
    }

    result *= 2 * pow(M_PI, 5.0/2.0) / (p*q*sqrt(p + q));

    return result;
}
