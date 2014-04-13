#ifndef GAUSSIANTYPEOVERLAPINTEGRAL_H
#define GAUSSIANTYPEOVERLAPINTEGRAL_H

#include <armadillo>

class HermiteExpansionCoefficient;
class GaussianPrimitiveOrbital;
class Vector3;

using namespace arma;

class GaussianOverlapIntegral
{
public:
    GaussianOverlapIntegral(Vector3 corePositionA, Vector3 corePositionB,
                            const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    GaussianOverlapIntegral(double exponentSum, HermiteExpansionCoefficient *hermiteExpansionCoefficient);
    virtual ~GaussianOverlapIntegral();


    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double overlapIntegral(int dim, int i, int j);
    double overlapIntegral(const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
protected:
    double m_exponentSum;

    HermiteExpansionCoefficient *m_hermiteExpansionCoefficient;

    bool m_isResponsibleForDeletingHermiteExpansionObject;
};

#endif // GAUSSIANTYPEOVERLAPINTEGRAL_H
