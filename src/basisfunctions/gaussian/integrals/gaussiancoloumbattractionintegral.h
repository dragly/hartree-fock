#ifndef GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H
#define GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H

#include <armadillo>
#include <math/hermiteexpansioncoefficient.h>
#include <hermiteintegral.h>
using namespace arma;

class HermiteIntegral;
class HermiteExpansionCoefficient;

class GaussianColoumbAttractionIntegral
{
public:
    GaussianColoumbAttractionIntegral(int angularMomentumMax);
    GaussianColoumbAttractionIntegral(const rowvec &corePositionA, const rowvec &corePositionB, const rowvec &corePositionC,
                                      double exponentA, double exponentB,
                                      int angularMomentumMax);

    double coloumbAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    void set(const rowvec &corePositionA, const rowvec &corePositionB, const rowvec &corePositionC, double exponentA, double exponentB);
protected:
    HermiteExpansionCoefficient m_hermiteExpansionCoefficient;
    HermiteIntegral m_hermiteIntegral;
    double m_exponentSum;
};

#endif // GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H
