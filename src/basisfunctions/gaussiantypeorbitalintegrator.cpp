#include "gaussiantypeorbitalintegrator.h"

#include <src/math/boysfunction.h>

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

    cout << "p=" << p << endl;
    cout << "mu=" << mu << endl;
    cout << "AB=" << AB << endl;
    cout << "P=" << P << endl;

    cout << "PA = " << PA << endl;
    cout << "PB = " << PB << endl;

    m_totalExponent = p;
    m_reducedExponent = mu;
    m_relativeSeparation = AB;
    m_centerOfMass = P;
    m_centerOfMassDiffA = PA;
    m_centerOfMassDiffB = PB;

    setupE();
    rowvec testPosition = {0,0,0};
    setupR(testPosition);
}

bool GaussianTypeOrbitalIntegrator::checkIndexCombinationForE(int iA, int iB, int t) {
    if(t < 0 || t > (iA + iB) || iA < 0 || iB < 0) {
        return false;
    } else {
        return true;
    }
}

void GaussianTypeOrbitalIntegrator::setupE() {
    int maxL = max(maxAngularMomentumA(), maxAngularMomentumB());
    for(int dim = 0; dim < 3; dim++) {
        int maxiA = maxL + 3;
        int maxiB = maxL + 3;
        // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
        // needed values available for the algorithm.
        int maxt = maxL + 3;
        m_E[dim] = zeros(maxiA, maxiB, maxt);

        double p = m_totalExponent;
        double mu = m_reducedExponent;
        double dim_AB = m_relativeSeparation(dim);
        double dim_PA = m_centerOfMassDiffA(dim);
        double dim_PB = m_centerOfMassDiffB(dim);

        //    vector<string> performed;

        // First row is special
        m_E[dim](0,0,0) = exp(-mu * dim_AB * dim_AB);
        cout << "KAB=" << m_E[dim](0,0,0) << endl;
        //    performed.push_back("000");
        for(int iB = 0; iB < maxiB; iB++) {
            for(int t = 0; t < maxt - 1; t++) {
                int iA = 0;
                if(iA == 0 && iB == 0 && t == 0) {
                    continue;
                }
                // p = previous, n = next
                // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
                int iBp = iB - 1;
                int tp = t - 1;
                int tn = t + 1;
                double E_iA_iBp_tp = 0;
                if(checkIndexCombinationForE(iA, iBp, tp)) {
                    E_iA_iBp_tp = m_E[dim](iA, iBp, tp);
                }
                double E_iA_iBp_t = 0;
                if(checkIndexCombinationForE(iA, iBp, t)) {
                    E_iA_iBp_t = m_E[dim](iA, iBp, t);
                }
                double E_iA_iBp_tn = 0;
                if(checkIndexCombinationForE(iA, iBp, tn)) {
                    E_iA_iBp_tn = m_E[dim](iA, iBp, tn);
                }
                m_E[dim](iA,iB,t) = 1 / (2*p) * E_iA_iBp_tp + dim_PB * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
            }
        }
        for(int iA = 1; iA < maxiA; iA++) {
            for(int iB = 0; iB < maxiB; iB++) {
                for(int t = 0; t < maxt - 1; t++) {
                    // p = previous, n = next
                    // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
                    int iAp = iA - 1;
                    int tp = t - 1;
                    int tn = t + 1;
                    double E_iAp_iB_tp = 0;
                    if(checkIndexCombinationForE(iAp, iB, tp)) {
                        E_iAp_iB_tp = m_E[dim](iAp, iB, tp);
                    }
                    double E_iAp_iB_t = 0;
                    if(checkIndexCombinationForE(iAp, iB, t)) {
                        E_iAp_iB_t = m_E[dim](iAp, iB, t);
                    }
                    double E_iAp_iB_tn = 0;
                    if(checkIndexCombinationForE(iAp, iB, tn)) {
                        E_iAp_iB_tn = m_E[dim](iAp, iB, tn);
                    }
                    m_E[dim](iA,iB,t) = 1 / (2*p) * E_iAp_iB_tp + dim_PA * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
                }
            }
        }
//        cout << m_E[dim];
    }
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
//    cube* E = m_E;
//    double a = m_exponentA;
//    double b = m_exponentB;
//    double p = a + b;
    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
    //    return E[0](iA,iB,0) * E[1](jA,jB,0) * E[2](kA,kB,0) * pow((M_PI / p), 3./2.);
}

double GaussianTypeOrbitalIntegrator::kineticIntegral(int dim, int iA, int iB) {
//    double a = m_exponentA;
    double b = m_exponentB;
//    double p = a + b;
//    const cube &E_dim = m_E[dim];
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

//double GaussianTypeOrbitalIntegrator::boysFunction(double arg) {
//    if (arg < 1.0E-6){
//        return 1.0;
//    } else {
//        arg = sqrt(arg);
//        double f = 1.0/arg * erf(arg) *sqrt(M_PI)/2.0;
//        return f;
//    }
//}

void GaussianTypeOrbitalIntegrator::setupR(const rowvec &corePositionC) {
//    for(int dim = 0; dim < 3; dim++) {
//        double p = m_totalExponent;
//        double mu = m_reducedExponent;
//        double dim_AB = m_relativeSeparation(dim);
//        double dim_P = m_centerOfMass(dim);
//        double dim_PA = m_centerOfMassDiffA(dim);
//        double dim_PB = m_centerOfMassDiffB(dim);
//        double R0000 = boysFunction();
//    }
//    for(int dim = 0; dim < 3; dim++) {
//    double a = m_exponentA;
    double p = m_totalExponent;
//    double mu = m_reducedExponent;
//    const rowvec &AB = m_relativeSeparation;
    const rowvec &C = corePositionC;
    const rowvec &P = m_centerOfMass;
    const rowvec &PC = P - C;
//    const rowvec &PA = m_centerOfMassDiffA;
//    const rowvec &PB = m_centerOfMassDiffB;
    cout << "dot = " << dot(PC, PC) << endl;
    double boysArg = p * dot(PC, PC);
    BoysFunction boysFunction(boysArg, 20);
    double F0 = boysFunction.result(0);
    double F1 = boysFunction.result(1);
    double F2 = boysFunction.result(2);
    double R0000 = F0;
    double R1000 = -2 * p * F1;
    double R2000 = (-2 * p)*(-2*p) * F2;
    cout << "R0000 = " << R0000 << endl;
    cout << "R1000 = " << R1000 << endl;
    cout << "R2000 = " << R2000 << endl;
    cube *E = m_E;
    cout << "Result = " << 2 * M_PI / p * R0000 * E[0](0,0,0) * E[1](0,0,0) * E[2](0,0,0) << endl;
}

//double GaussianTypeOrbitalIntegrator::coreIntegral() {

//}

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
    return m_maxAngularMomentumA;
}

void GaussianTypeOrbitalIntegrator::setMaxAngularMomentumA(const uint &maxAngularMomentum)
{
    m_maxAngularMomentumA = maxAngularMomentum;
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
