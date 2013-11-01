#include "gaussiantypeorbitalintegrator.h"

#include <math/boysfunction.h>
#include <math/boysfunctionintermediate.h>
#include <set>

using namespace std;
using namespace arma;

GaussianTypeOrbitalIntegrator::GaussianTypeOrbitalIntegrator() :
    m_exponentA(0),
    m_exponentB(0),
    m_totalExponent(0),
    m_reducedExponent(0),
    m_weightA(0),
    m_weightB(0),
    m_corePositionA(zeros<rowvec>(3)),
    m_corePositionB(zeros<rowvec>(3)),
    m_corePositionC(zeros<rowvec>(3)),
    m_isDirty(true)
{
    setMaxAngularMomentumA(0);
    setMaxAngularMomentumB(0);
}

void GaussianTypeOrbitalIntegrator::reset() {
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    double mu = a * b / (a + b);
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;

    rowvec AB = A - B;
    rowvec P = (a * A + b * B) / p;
    rowvec PA = P - A;
    rowvec PB = P - B;

    m_totalExponent = p;
    m_reducedExponent = mu;
    m_relativeSeparation = AB;
    m_centerOfMass = P;
    m_centerOfMassDiffA = PA;
    m_centerOfMassDiffB = PB;

    setupE();
    setupR();
}

void GaussianTypeOrbitalIntegrator::setupE() {
    m_E = HermiteExpansionCoefficient(m_exponentA, m_exponentB, m_corePositionA, m_corePositionB, m_angularMomentumAMax);
}

double GaussianTypeOrbitalIntegrator::overlapIntegral(int dim, int iA, int iB)
{
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    const cube &E_dim = m_E[dim];
    return E_dim(iA,iB,0) * sqrt(M_PI / p);
}


double GaussianTypeOrbitalIntegrator::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}

double GaussianTypeOrbitalIntegrator::coreColoumbIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double result = 0;
    const cube &E_x = m_E[0];
    const cube &E_y = m_E[1];
    const cube &E_z = m_E[2];
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    const HermiteIntegral &R = m_R;
    int tMax = iA + iB;
    int uMax = jA + jB;
    int vMax = kA + kB;
    for(int t = 0; t < tMax + 1; t++) {
        for(int u = 0; u < uMax + 1; u++) {
            for(int v = 0; v < vMax + 1; v++) {
                result += E_x(iA, iB, t) * E_y(jA, jB, u) * E_z(kA, kB, v) * R(0,t,u,v);
            }
        }
    }
    result *= 2*M_PI / p;
    return result;
}

double GaussianTypeOrbitalIntegrator::kineticIntegral(int dim, int iA, int iB) {
    double b = m_exponentB;
    double S_iA_iBnn = overlapIntegral(dim, iA, iB + 2);
    double S_iA_iB = overlapIntegral(dim, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= overlapIntegral(dim, iA, iB - 2);
    } else {
        S_iA_iBpp = 0;
    }
    return 4 * b * b * S_iA_iBnn - 2*b * (2*iB + 1) * S_iA_iB + iB * (iB + 1) * S_iA_iBpp;
}

double GaussianTypeOrbitalIntegrator::kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double T_iA_iB = kineticIntegral(0, iA, iB);
    double T_jA_jB = kineticIntegral(1, jA, jB);
    double T_kA_kB = kineticIntegral(2, kA, kB);

    double S_iA_iB = overlapIntegral(0, iA, iB);
    double S_jA_jB = overlapIntegral(1, jA, jB);
    double S_kA_kB = overlapIntegral(2, kA, kB);

    double result = T_iA_iB * S_jA_jB * S_kA_kB + S_iA_iB * T_jA_jB * S_kA_kB + S_iA_iB * S_jA_jB * T_kA_kB;
    result *= -0.5;
    return result; // TODO: Is there an error in the slides here?
}

bool hasIndices(const vector<int>& neededIndices, const vector<vector<int>>& calculatedIndices) {
    bool found = false;
    for(vector<int> indices : calculatedIndices) {
        if(indices[0] == neededIndices[0]
                && indices[1] == neededIndices[1]
                && indices[2] == neededIndices[2]
                && indices[3] == neededIndices[3]) {
            found = true;
        }
    }
    if(!found) {
        cout << "ERROR! Needed element not found" << endl;
        return false;
    }
    return true;
}

void GaussianTypeOrbitalIntegrator::setupR() {
    double p = m_totalExponent;
    const rowvec &C = m_corePositionC;
    const rowvec &P = m_centerOfMass;
    const rowvec &PC = P - C;
    m_R = HermiteIntegral(p, PC, m_angularMomentumAMax + m_angularMomentumBMax);
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
    (void)l;
    return zeros(0);
}

uint GaussianTypeOrbitalIntegrator::maxAngularMomentumA() const
{
    return m_angularMomentumAMax;
}

void GaussianTypeOrbitalIntegrator::setMaxAngularMomentumA(const uint &maxAngularMomentum)
{
    m_angularMomentumAMax = maxAngularMomentum;
    if(maxAngularMomentum > 4) {
        cout << "WARNING: Creating gaussian type orbitals for angular momentums over 4. Subshells s,p,d,f,g,?" << endl;
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
double GaussianTypeOrbitalIntegrator::exponentB() const
{
    return m_exponentB;
}

void GaussianTypeOrbitalIntegrator::setExponentB(double exponentB)
{
    m_exponentB = exponentB;
}

double GaussianTypeOrbitalIntegrator::exponentA() const
{
    return m_exponentA;
}

void GaussianTypeOrbitalIntegrator::setExponentA(double exponentA)
{
    m_exponentA = exponentA;
}

const vector<urowvec> &GaussianTypeOrbitalIntegrator::combinationsB() const
{
    return m_combinationsB;
}

void GaussianTypeOrbitalIntegrator::regenerateCombinations(bool isA = true) {
    m_isDirty = true;
    vector<urowvec> combinations;
    uint l = m_angularMomentumAMax;
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
rowvec GaussianTypeOrbitalIntegrator::corePositionC() const
{
    return m_corePositionC;
}

void GaussianTypeOrbitalIntegrator::setCorePositionC(const rowvec &corePositionC)
{
    m_corePositionC = corePositionC;
}


uint GaussianTypeOrbitalIntegrator::maxAngularMomentumB() const
{
    return m_angularMomentumBMax;
}

void GaussianTypeOrbitalIntegrator::setMaxAngularMomentumB(const uint &maxAngularMomentumB)
{
    m_angularMomentumBMax = maxAngularMomentumB;
}

vector<urowvec> GaussianTypeOrbitalIntegrator::combinationsA() const
{
    return m_combinationsA;
}
