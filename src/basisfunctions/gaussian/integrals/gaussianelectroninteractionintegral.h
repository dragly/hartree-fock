#ifndef GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
#define GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H

#include <armadillo>
#include <hermiteintegral.h>
#include <math/hermiteexpansioncoefficient.h>

class Vector3;

using namespace arma;

class GaussianPrimitiveOrbital;

class GaussianElectronInteractionIntegral
{
public:
    GaussianElectronInteractionIntegral(int singleAngularMomentumMax);

    ~GaussianElectronInteractionIntegral();

    void reset(int singleAngularMomentumMax);

    void set(const Vector3 &corePositionA, const Vector3 &corePositionB,
             const Vector3 &corePositionC, const Vector3 &corePositionD,
             const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB,
             const GaussianPrimitiveOrbital &primitiveC, const GaussianPrimitiveOrbital &primitiveD);
    void setAB(const Vector3 &corePositionA, const Vector3 &corePositionB,
               const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    void setCD(const Vector3 &corePositionC, const Vector3 &corePositionD,
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
    Vector3 m_centerOfMassP;
    Vector3 m_centerOfMassQ;
    double m_exponentP;
    double m_exponentQ;

    const GaussianPrimitiveOrbital *m_primitiveA;
    const GaussianPrimitiveOrbital *m_primitiveB;
    const GaussianPrimitiveOrbital *m_primitiveC;
    const GaussianPrimitiveOrbital *m_primitiveD;
};

#endif // GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
