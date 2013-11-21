#ifndef GAUSSIANTYPEORBITALINTEGRATOR_H
#define GAUSSIANTYPEORBITALINTEGRATOR_H

#include <armadillo>
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
private:
    void regenerateCombinationsA();
    void regenerateCombinationsB();
    void regenerateCombinations(bool isA);

    double m_exponentA;
    double m_exponentB;

    double m_weightA;
    double m_weightB;

    rowvec m_corePositionA;
    rowvec m_corePositionB;

    uint m_maxAngularMomentumA;
    uint m_maxAngularMomentumB;

    bool m_isDirty;

    cube m_E[3]; // t, i, j

    vector<urowvec> m_combinationsA;
    vector<urowvec> m_combinationsB;
    void setupE();
};

#endif // GAUSSIANTYPEORBITALINTEGRATOR_H
