#ifndef GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
#define GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H

#include <armadillo>
using namespace arma;

class HermiteIntegral;
class HermiteExpansionCoefficient;

class GaussianTypeElectronInteractionIntegral
{
public:
    GaussianTypeElectronInteractionIntegral(const rowvec& corePositionA, const rowvec& corePositionB,
                                            const rowvec& corePositionC, const rowvec& corePositionD,
                                            double exponentA, double exponentB,
                                            double exponentC, double exponentD,
                                            int angularMomentumMax);

    GaussianTypeElectronInteractionIntegral(double exponentA, double exponentB, double exponentC, double exponentD,
                                            int angularMomentumMax,
                                            HermiteExpansionCoefficient* hermiteExpansionCoefficientAB,
                                            HermiteExpansionCoefficient* hermiteExpansionCoefficientCD,
                                            HermiteIntegral* hermiteIntegral);

    double electronInteractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB, int iC, int jC, int kC, int iD, int jD, int kD);
protected:
    HermiteIntegral* m_hermiteIntegral;
    HermiteExpansionCoefficient* m_hermiteExpansionCoefficientAB;
    HermiteExpansionCoefficient* m_hermiteExpansionCoefficientCD;
    int m_angularMomentumMax;
    double m_exponentP;
    double m_exponentQ;
};

#endif // GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
