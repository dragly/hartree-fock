#ifndef GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H
#define GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H

#include <armadillo>
#include <math/hermiteexpansioncoefficient.h>
#include <hermiteintegral.h>

class Vector3;

using namespace arma;

class GaussianPrimitiveOrbital;

class GaussianColoumbAttractionIntegral
{
public:
    GaussianColoumbAttractionIntegral(int angularMomentumMax);

    double coloumbAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    void set(const Vector3 &corePositionA, const Vector3 &corePositionB, const Vector3 &corePositionC,
             const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    double coloumbAttractionIntegral(const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
protected:
    HermiteExpansionCoefficient m_hermiteExpansionCoefficient;
    HermiteIntegral m_hermiteIntegral;
    double m_exponentSum;
};

#endif // GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H
