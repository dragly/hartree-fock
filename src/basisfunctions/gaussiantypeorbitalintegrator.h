#ifndef GAUSSIANTYPEORBITALINTEGRATOR_H
#define GAUSSIANTYPEORBITALINTEGRATOR_H

#include <armadillo>
#include <hermiteintegral.h>
#include <math/hermiteexpansioncoefficient.h>
#include <vector>

using namespace arma;
using namespace std;

class GaussianTypeOrbitalIntegrator
{
public:
    GaussianTypeOrbitalIntegrator();

    void setElectronPosition(const rowvec &electronPosition);
    void updateCoreElectronVector();

    rowvec corePositionB() const;
    void setCorePositionB(const rowvec &corePositionB);

    rowvec corePositionA() const;
    void setCorePositionA(const rowvec &corePositionA);

    rowvec overlapIntegrals(int maxAngularMomentum);

    uint maxAngularMomentumA() const;
    void setMaxAngularMomentumA(const uint &maxAngularMomentumA);

    uint maxAngularMomentumB() const;
    void setMaxAngularMomentumB(const uint &maxAngularMomentumB);

    vector<urowvec> combinationsA() const;
    const vector<urowvec> &combinationsB() const;


    double exponentA() const;
    void setExponentA(double exponentA);

    double exponentB() const;
    void setExponentB(double exponentB);

    void reset();
    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double overlapIntegral(int dim, int i, int j);

    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double kineticIntegral(int dim, int iA, int iB);
//    bool checkIndexCombinationForE(int iA, int iB, int t);
    //    double boysFunction(double arg);
    double coreColoumbIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec corePositionC() const;
    void setCorePositionC(const rowvec &corePositionC);
    void setupR();

protected:
    void regenerateCombinationsA();
    void regenerateCombinationsB();
    void regenerateCombinations(bool isA);

    double m_exponentA;
    double m_exponentB;

    double m_totalExponent;

    double m_reducedExponent;

    double m_weightA;
    double m_weightB;

    rowvec m_corePositionA;
    rowvec m_corePositionB;
    rowvec m_corePositionC;

    rowvec m_centerOfMass;
    rowvec m_centerOfMassDiffA;
    rowvec m_centerOfMassDiffB;

    rowvec m_relativeSeparation;

    uint m_angularMomentumAMax;
    uint m_angularMomentumBMax;

    bool m_isDirty;

    HermiteExpansionCoefficient m_E; // t, i, j
    HermiteIntegral m_R;

    vector<urowvec> m_combinationsA;
    vector<urowvec> m_combinationsB;
    void setupE();
//    void setupR();
};

#endif // GAUSSIANTYPEORBITALINTEGRATOR_H
