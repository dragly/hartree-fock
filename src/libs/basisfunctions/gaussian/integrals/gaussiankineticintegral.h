#ifndef GAUSSIANTYPEKINETICINTEGRAL_H
#define GAUSSIANTYPEKINETICINTEGRAL_H

#include <armadillo>

using namespace arma;

class HermiteExpansionCoefficient;

class GaussianKineticIntegral
{
public:
    explicit GaussianKineticIntegral(rowvec corePositionA, rowvec corePositionB, double exponentA, double exponentB, int angularMomentumMax);
    explicit GaussianKineticIntegral(double exponentA, double exponentB, HermiteExpansionCoefficient *hermiteExpansionCoefficient);
    virtual ~GaussianKineticIntegral();
    double kineticIntegral(int dim, int iA, int iB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
protected:
    double m_exponentB;
    double m_exponentSum;

    HermiteExpansionCoefficient *m_hermiteExpansionCoefficient;

    bool m_isResponsibleForDeletingHermiteExpansionObject;
};

#endif // GAUSSIANTYPEKINETICINTEGRAL_H
