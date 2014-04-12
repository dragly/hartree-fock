#ifndef GAUSSIANTYPEOVERLAPINTEGRAL_H
#define GAUSSIANTYPEOVERLAPINTEGRAL_H

#include <armadillo>

using namespace arma;

class HermiteExpansionCoefficient;
class GaussianPrimitiveOrbital;

class GaussianOverlapIntegral
{
public:
    GaussianOverlapIntegral(rowvec corePositionA, rowvec corePositionB,
                            const GaussianPrimitiveOrbital &primitiveA, const GaussianPrimitiveOrbital &primitiveB);
    GaussianOverlapIntegral(double exponentSum, HermiteExpansionCoefficient *hermiteExpansionCoefficient);
    virtual ~GaussianOverlapIntegral();


    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double overlapIntegral(int dim, int i, int j);
protected:
    double m_exponentSum;

    HermiteExpansionCoefficient *m_hermiteExpansionCoefficient;

    bool m_isResponsibleForDeletingHermiteExpansionObject;
};

#endif // GAUSSIANTYPEOVERLAPINTEGRAL_H
