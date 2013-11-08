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
    void setAB(const rowvec &corePositionA, const rowvec &corePositionB, double exponentA, double exponentB);
    void setCD(const rowvec &corePositionC, const rowvec &corePositionD, double exponentC, double exponentD);

    double electronInteractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB, int iC, int jC, int kC, int iD, int jD, int kD);
protected:
    HermiteIntegral m_hermiteIntegral;
    HermiteExpansionCoefficient m_hermiteExpansionCoefficientAB;
    HermiteExpansionCoefficient m_hermiteExpansionCoefficientCD;
    int m_angularMomentumMax;
    rowvec m_centerOfMassP;
    rowvec m_centerOfMassQ;
    double m_exponentP;
    double m_exponentQ;
};

#endif // GAUSSIANTYPEELECTRONINTERACTIONINTEGRAL_H
