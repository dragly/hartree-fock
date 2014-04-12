#ifndef GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
#define GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H

#include <armadillo>
#include <hermiteintegral.h>
#include <math/hermiteexpansioncoefficient.h>

using namespace arma;

class GaussianPrimitiveOrbital;

class GaussianElectronInteractionIntegral
{
public:
    GaussianElectronInteractionIntegral(int singleAngularMomentumMax);

    ~GaussianElectronInteractionIntegral();

    void reset(int singleAngularMomentumMax);

    void set(const rowvec &corePositionA, const rowvec &corePositionB,
             const rowvec &corePositionC, const rowvec &corePositionD,
             const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB,
             const GaussianPrimitiveOrbital &primitiveC, const GaussianPrimitiveOrbital &primitiveD);
    void setAB(const rowvec &corePositionA, const rowvec &corePositionB,
               const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    void setCD(const rowvec &corePositionC, const rowvec &corePositionD,
               const GaussianPrimitiveOrbital &primitiveC, const GaussianPrimitiveOrbital &primitiveD);
    double electronInteractionIntegral(int iA, int jA, int kA,
                                       int iB, int jB, int kB,
                                       int iC, int jC, int kC,
                                       int iD, int jD, int kD);
    double electronInteractionIntegral(const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB, const GaussianPrimitiveOrbital &primitiveC, const GaussianPrimitiveOrbital &primitiveD);
protected:
    HermiteIntegral m_hermiteIntegral;
    HermiteExpansionCoefficient m_hermiteExpansionCoefficientAB;
    HermiteExpansionCoefficient m_hermiteExpansionCoefficientCD;
    int m_angularMomentumMax;
    rowvec m_centerOfMassP;
    rowvec m_centerOfMassQ;
    double m_exponentP;
    double m_exponentQ;
};

#endif // GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
