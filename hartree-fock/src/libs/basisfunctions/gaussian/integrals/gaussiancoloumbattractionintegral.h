#ifndef GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H
#define GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H

#include <armadillo>
using namespace arma;

class HermiteIntegral;
class HermiteExpansionCoefficient;

class GaussianColoumbAttractionIntegral
{
public:
    GaussianColoumbAttractionIntegral(rowvec corePositionA, rowvec corePositionB, rowvec corePositionC, double exponentA, double exponentB, int angularMomentumMax);
    GaussianColoumbAttractionIntegral(double exponentSum, HermiteExpansionCoefficient *hermiteExpansionCoefficient, HermiteIntegral *hermiteIntegral);
    ~GaussianColoumbAttractionIntegral();

    double coloumbAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB);

protected:
    HermiteExpansionCoefficient* m_hermiteExpansionCoefficient;
    HermiteIntegral *m_hermiteIntegral;
    double m_exponentSum;
    bool m_isResponsibleForDeletingHermiteObjects;
};

#endif // GAUSSIANTYPECOLOUMBATTRACTIONINTEGRAL_H