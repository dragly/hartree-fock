#ifndef GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
#define GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H

#include <armadillo>
#include <hermiteintegral.h>
#include <math/hermiteexpansioncoefficient.h>

using namespace arma;


class GaussianElectronInteractionIntegral
{
public:
    GaussianElectronInteractionIntegral(int angularMomentumMax);
    GaussianElectronInteractionIntegral(const rowvec& corePositionA, const rowvec& corePositionB,
                                            const rowvec& corePositionC, const rowvec& corePositionD,
                                            double exponentA, double exponentB,
                                            double exponentC, double exponentD,
                                            int angularMomentumMax);

    ~GaussianElectronInteractionIntegral();

//    GaussianElectronInteractionIntegral(double exponentA, double exponentB, double exponentC, double exponentD,
//                                            int angularMomentumMax,
//                                            HermiteExpansionCoefficient* hermiteExpansionCoefficientAB,
//                                            HermiteExpansionCoefficient* hermiteExpansionCoefficientCD,
//                                            HermiteIntegral* hermiteIntegral);

    void reset(int angularMomentumMax);
    void set(const rowvec& corePositionA, const rowvec& corePositionB,
          const rowvec& corePositionC, const rowvec& corePositionD,
          double exponentA, double exponentB,
          double exponentC, double exponentD);
    void setAC(const rowvec &corePositionA, const rowvec &corePositionC, double exponentA, double exponentC);
    void setBD(const rowvec &corePositionB, const rowvec &corePositionD, double exponentB, double exponentD);

    double electronInteractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB, int iC, int jC, int kC, int iD, int jD, int kD);
protected:
    HermiteIntegral m_hermiteIntegral;
    HermiteExpansionCoefficient m_hermiteExpansionCoefficientAC;
    HermiteExpansionCoefficient m_hermiteExpansionCoefficientBD;
    int m_angularMomentumMax;
    rowvec m_centerOfMassP;
    rowvec m_centerOfMassQ;
    double m_exponentP;
    double m_exponentQ;
};

#endif // GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
