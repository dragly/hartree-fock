#ifndef GAUSSIANTYPEKINETICINTEGRAL_H
#define GAUSSIANTYPEKINETICINTEGRAL_H

#include <armadillo>
#include <math/hermiteexpansioncoefficient.h>

using namespace arma;

class HermiteExpansionCoefficient;

class GaussianKineticIntegral
{
public:
    explicit GaussianKineticIntegral(int angularMomentumMax);
    explicit GaussianKineticIntegral(const rowvec& corePositionA, const rowvec& corePositionB, double exponentA, double exponentB, int angularMomentumMax);
    void set(const rowvec &corePositionA, const rowvec &corePositionB,
             double exponentA, double exponentB);
    double kineticIntegral(int dim, int iA, int iB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
protected:
    double m_exponentB;
    double m_exponentSum;

    HermiteExpansionCoefficient m_hermiteExpansionCoefficient;
};

#endif // GAUSSIANTYPEKINETICINTEGRAL_H
