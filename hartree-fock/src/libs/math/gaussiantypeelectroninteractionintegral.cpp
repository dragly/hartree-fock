#include "gaussiantypeelectroninteractionintegral.h"

#include <math/hermiteexpansioncoefficient.h>
#include <hermiteintegral.h>

GaussianTypeElectronInteractionIntegral::GaussianTypeElectronInteractionIntegral(const rowvec &corePositionA, const rowvec &corePositionB, const rowvec &corePositionC, const rowvec &corePositionD, double exponentA, double exponentB, double exponentC, double exponentD, int angularMomentumMax) :
    GaussianTypeElectronInteractionIntegral(exponentA, exponentB, exponentC, exponentD, angularMomentumMax, 0, 0, 0)
{
    double a = exponentA;
    double b = exponentB;
    double c = exponentC;
    double d = exponentD;
    double p = a + b;
    double q = c + d;
    rowvec P = (a * corePositionA + b * corePositionB) / (a + b);
    rowvec Q = (c * corePositionC + d * corePositionD) / (c + d);
    rowvec PQ = P - Q;
    double alpha = p*q/(p+q);
    m_hermiteIntegral = new HermiteIntegral(alpha, PQ, 4 * angularMomentumMax);
    m_hermiteExpansionCoefficientAB = new HermiteExpansionCoefficient(a, b, corePositionA, corePositionB, angularMomentumMax);
    m_hermiteExpansionCoefficientCD = new HermiteExpansionCoefficient(c, d, corePositionC, corePositionD, angularMomentumMax);
}

GaussianTypeElectronInteractionIntegral::GaussianTypeElectronInteractionIntegral(double exponentA, double exponentB,
                                                                                 double exponentC, double exponentD,
                                                                                 int angularMomentumMax,
                                                                                 HermiteExpansionCoefficient *hermiteExpansionCoefficientAB, HermiteExpansionCoefficient *hermiteExpansionCoefficientCD, HermiteIntegral *hermiteIntegral) :
    m_hermiteIntegral(hermiteIntegral),
    m_hermiteExpansionCoefficientAB(hermiteExpansionCoefficientAB),
    m_hermiteExpansionCoefficientCD(hermiteExpansionCoefficientCD),
    m_angularMomentumMax(angularMomentumMax)
{
    double a = exponentA;
    double b = exponentB;
    double c = exponentC;
    double d = exponentD;
    double p = a + b;
    double q = c + d;
    m_exponentP = p;
    m_exponentQ = q;
}

double GaussianTypeElectronInteractionIntegral::electronInteractionIntegral(int iA, int jA, int kA,
                                                                            int iB, int jB, int kB,
                                                                            int iC, int jC, int kC,
                                                                            int iD, int jD, int kD) {
    double result;
    double p = m_exponentP;
    double q = m_exponentQ;

    const HermiteIntegral &R = (*m_hermiteIntegral);
    const HermiteExpansionCoefficient &Eab = (*m_hermiteExpansionCoefficientAB);
    const HermiteExpansionCoefficient &Ecd = (*m_hermiteExpansionCoefficientCD);
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
                            result *= R(0, t + tau, u + nu, v + phi);
                            result *= pow(-1, tau + nu + phi);
                            result *= Ecd(iC, jC, kC, iD, jD, kD, tau, nu, phi);
                        }
                    }
                }
                result *= Eab(iA, jA, kA, iB, jB, kB, t, u, v);
            }
        }
    }

    result *= 2 * pow(M_PI, 5.0/2.0) / (p*q*sqrt(p + q));

    return result;
}
