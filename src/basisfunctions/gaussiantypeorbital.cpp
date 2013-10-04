#include "gaussiantypeorbital.h"

using namespace std;
using namespace arma;

GaussianTypeOrbitalIntegrator::GaussianTypeOrbitalIntegrator() :
    m_orbitalExponentA(0),
    m_orbitalExponentB(0),
    m_isDirty(true),
    m_corePositionA(zeros<rowvec>(3)),
    m_corePositionB(zeros<rowvec>(3))
{
    setMaxAngularMomentumA(0);
    setMaxAngularMomentumB(0);
}

void GaussianTypeOrbitalIntegrator::reset() {
    setupE();
}

void GaussianTypeOrbitalIntegrator::setupE() {
    uint maxL = max(maxAngularMomentumA(), maxAngularMomentumB());
    m_Ex = zeros(maxL + 1, maxL + 1, maxL + 1);
    double a = m_orbitalExponentA;
    double b = m_orbitalExponentB;
    double mu = a * b / (a + b);
    double xDiff = m_corePositionA(0) - m_corePositionB(0);
    m_Ex(0,0,0) = exp(-mu * xDiff);

    cout << m_Ex;
}

rowvec GaussianTypeOrbitalIntegrator::corePositionB() const
{
    return m_corePositionB;
}

void GaussianTypeOrbitalIntegrator::setCorePositionB(const rowvec &corePositionB)
{
    m_corePositionB = corePositionB;
}
rowvec GaussianTypeOrbitalIntegrator::corePositionA() const
{
    return m_corePositionA;
}

void GaussianTypeOrbitalIntegrator::setCorePositionA(const rowvec &corePositionA)
{
    m_corePositionA = corePositionA;
}

rowvec GaussianTypeOrbitalIntegrator::overlapIntegrals(int maxAngularMomentum) {
    int l = maxAngularMomentum;
    l = 2;
    return zeros(0);
}

uint GaussianTypeOrbitalIntegrator::maxAngularMomentumA() const
{
    return m_maxAngularMomentumA;
}

void GaussianTypeOrbitalIntegrator::setMaxAngularMomentumA(const uint &maxAngularMomentum)
{
    m_maxAngularMomentumA = maxAngularMomentum;
    if(maxAngularMomentum > 4) {
        cout << "WARNING: Creating GTO orbitals for angular momentums over 4. Subshells s,p,d,f,g,?" << endl;
    }
    regenerateCombinationsA();
}

void GaussianTypeOrbitalIntegrator::regenerateCombinationsA() {
    regenerateCombinations(true);
}

void GaussianTypeOrbitalIntegrator::regenerateCombinationsB()
{
    regenerateCombinations(false);
}
double GaussianTypeOrbitalIntegrator::orbitalExponentB() const
{
    return m_orbitalExponentB;
}

void GaussianTypeOrbitalIntegrator::setOrbitalExponentB(double orbitalExponentB)
{
    m_orbitalExponentB = orbitalExponentB;
}

double GaussianTypeOrbitalIntegrator::orbitalExponentA() const
{
    return m_orbitalExponentA;
}

void GaussianTypeOrbitalIntegrator::setOrbitalExponentA(double orbitalExponentA)
{
    m_orbitalExponentA = orbitalExponentA;
}


vector<urowvec> GaussianTypeOrbitalIntegrator::combinationsB() const
{
    return m_combinationsB;
}

void GaussianTypeOrbitalIntegrator::regenerateCombinations(bool isA = true) {
    m_isDirty = true;
    vector<urowvec> combinations;
    uint l = m_maxAngularMomentumA;
    for(uint i = 0; i <= l; i++) {
        for(uint j = 0; j <= l; j++) {
            for(uint k = 0; k <= l; k++) {
                if(i + j + k > l) {
                    continue;
                }
                urowvec combination = {i, j, k};
                combinations.push_back(combination);
            }
        }
    }
    if(isA) {
        m_combinationsA = combinations;
    } else {
        m_combinationsB = combinations;
    }
}

uint GaussianTypeOrbitalIntegrator::maxAngularMomentumB() const
{
    return m_maxAngularMomentumB;
}

void GaussianTypeOrbitalIntegrator::setMaxAngularMomentumB(const uint &maxAngularMomentumB)
{
    m_maxAngularMomentumB = maxAngularMomentumB;
}

vector<urowvec> GaussianTypeOrbitalIntegrator::combinationsA() const
{
    return m_combinationsA;
}
