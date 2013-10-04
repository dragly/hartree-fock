#ifndef GAUSSIANTYPEORBITAL_H
#define GAUSSIANTYPEORBITAL_H

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
    vector<urowvec> combinationsB() const;


    double orbitalExponentA() const;
    void setOrbitalExponentA(double orbitalExponentA);

    double orbitalExponentB() const;
    void setOrbitalExponentB(double orbitalExponentB);

    void reset();
private:
    void regenerateCombinationsA();
    void regenerateCombinationsB();
    void regenerateCombinations(bool isA);

    rowvec m_corePositionA;
    rowvec m_corePositionB;
    urowvec m_xyzPowers;
    uint m_maxAngularMomentumA;
    uint m_maxAngularMomentumB;

    double m_orbitalExponentA;
    double m_orbitalExponentB;

    bool m_isDirty;

    cube m_Ex; // t, i, j

    vector<urowvec> m_combinationsA;
    vector<urowvec> m_combinationsB;
    void setupE();
};

#endif // GAUSSIANTYPEORBITAL_H
