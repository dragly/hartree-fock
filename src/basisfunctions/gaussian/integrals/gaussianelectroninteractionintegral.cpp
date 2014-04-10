#include "gaussianelectroninteractionintegral.h"

#include <math/hermiteexpansioncoefficient.h>
#include <hermiteintegral.h>

GaussianElectronInteractionIntegral::GaussianElectronInteractionIntegral(int angularMomentumMax) :
    m_hermiteIntegral(4*angularMomentumMax+1),
    m_hermiteExpansionCoefficientAB(angularMomentumMax+1),
    m_hermiteExpansionCoefficientCD(angularMomentumMax+1)
{
    reset(angularMomentumMax);
}

GaussianElectronInteractionIntegral::GaussianElectronInteractionIntegral(const rowvec &corePositionA, const rowvec &corePositionB,
                                                                         const rowvec &corePositionC, const rowvec &corePositionD,
                                                                         double exponentA, double exponentB, double exponentC, double exponentD, int angularMomentumMax) :
    GaussianElectronInteractionIntegral(angularMomentumMax)
{
    set(corePositionA, corePositionB, corePositionC, corePositionD, exponentA, exponentB, exponentC, exponentD);
}

//GaussianElectronInteractionIntegral::GaussianElectronInteractionIntegral(double exponentA, double exponentB,
//                                                                                 double exponentC, double exponentD,
//                                                                                 int angularMomentumMax,
//                                                                                 HermiteExpansionCoefficient *hermiteExpansionCoefficientAB, HermiteExpansionCoefficient *hermiteExpansionCoefficientCD, HermiteIntegral *hermiteIntegral) :
//    m_hermiteIntegral(hermiteIntegral),
//    m_hermiteExpansionCoefficientAB(hermiteExpansionCoefficientAB),
//    m_hermiteExpansionCoefficientCD(hermiteExpansionCoefficientCD),
//    m_angularMomentumMax(angularMomentumMax),
//    m_isResponsibleForFreeingIntegralsAndCoefficients(false)
//{
//    double a = exponentA;
//    double b = exponentB;
//    double c = exponentC;
//    double d = exponentD;
//    double p = a + b;
//    double q = c + d;
//    m_exponentP = p;
//    m_exponentQ = q;
//}

void GaussianElectronInteractionIntegral::reset(int angularMomentumMax)
{
    m_hermiteIntegral.reset(4*angularMomentumMax+1);
    m_hermiteExpansionCoefficientAB.reset(angularMomentumMax+1);
    m_hermiteExpansionCoefficientCD.reset(angularMomentumMax+1);
}

void GaussianElectronInteractionIntegral::set(const rowvec &corePositionA, const rowvec &corePositionB,
                                              const rowvec &corePositionC, const rowvec &corePositionD,
                                              double exponentA, double exponentB, double exponentC, double exponentD)
{
    setAB(corePositionA, corePositionB, exponentA, exponentB);
    setCD(corePositionC, corePositionD, exponentC, exponentD);
}

void GaussianElectronInteractionIntegral::setAB(const rowvec &corePositionA, const rowvec &corePositionB,
                                                double exponentA, double exponentB)
{
    double a = exponentA;
    double b = exponentB;
//    double c = exponentC;
//    double d = exponentD;
    double p = a + b;
    m_exponentP = p;
//    double q = c + d;
    rowvec P = (a * corePositionA + b * corePositionB) / (a + b);
    m_centerOfMassP = P;
//    rowvec Q = (c * corePositionC + d * corePositionD) / (c + d);
//    rowvec PQ = P - Q;
//    double alpha = p*q/(p+q);
//    m_hermiteIntegral.set(alpha, PQ);
    m_hermiteExpansionCoefficientAB.set(a, b, corePositionA, corePositionB);
//    m_hermiteExpansionCoefficientCD.set(c, d, corePositionC, corePositionD);
}

void GaussianElectronInteractionIntegral::setCD(const rowvec &corePositionC, const rowvec &corePositionD,
                                                double exponentC, double exponentD,
                                                int angularMomentum)
{
//    double a = exponentA;
//    double b = exponentB;
    double c = exponentC;
    double d = exponentD;
    double p = m_exponentP;
    double q = c + d;
    m_exponentQ = q;
    rowvec P = m_centerOfMassP;
    rowvec Q = (c * corePositionC + d * corePositionD) / (c + d);
    m_centerOfMassQ = Q;
    rowvec PQ = P - Q;
    double alpha = p*q/(p+q);
//    m_hermiteIntegral.set(alpha, PQ);
    int dimension = 0;
    if(angularMomentum != -1) {
        dimension = angularMomentum + 1;
    } else {
        dimension = -1;
    }
    m_hermiteIntegral.set(alpha, PQ, dimension);
//    m_hermiteExpansionCoefficientAB.set(a, b, corePositionA, corePositionB);
    m_hermiteExpansionCoefficientCD.set(c, d, corePositionC, corePositionD);
}

GaussianElectronInteractionIntegral::~GaussianElectronInteractionIntegral()
{
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
