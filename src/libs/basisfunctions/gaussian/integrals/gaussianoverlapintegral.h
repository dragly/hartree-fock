#ifndef GAUSSIANTYPEOVERLAPINTEGRAL_H
#define GAUSSIANTYPEOVERLAPINTEGRAL_H

#include <armadillo>

using namespace arma;

class HermiteExpansionCoefficient;

class GaussianOverlapIntegral
{
public:
    GaussianOverlapIntegral(double exponentSum, HermiteExpansionCoefficient *hermiteExpansionCoefficient);
    GaussianOverlapIntegral(rowvec corePositionA, rowvec corePositionB,
                            double exponentA, double exponentB,
                            int angularMomentumMax);
    virtual ~GaussianOverlapIntegral();


    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double overlapIntegral(int dim, int i, int j);
protected:
    double m_exponentSum;

    HermiteExpansionCoefficient *m_hermiteExpansionCoefficient;

    bool m_isResponsibleForDeletingHermiteExpansionObject;
};

#endif // GAUSSIANTYPEOVERLAPINTEGRAL_H
