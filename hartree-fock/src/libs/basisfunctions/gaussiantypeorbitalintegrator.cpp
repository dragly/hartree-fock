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

//    cout << "p=" << p << endl;
//    cout << "mu=" << mu << endl;
//    cout << "AB=" << AB << endl;
//    cout << "P=" << P << endl;

//    cout << "PA = " << PA << endl;
//    cout << "PB = " << PB << endl;

    m_totalExponent = p;
    m_reducedExponent = mu;
    m_relativeSeparation = AB;
    m_centerOfMass = P;
    m_centerOfMassDiffA = PA;
    m_centerOfMassDiffB = PB;

    setupE();
    setupR();
    cout << m_R[0] << endl;
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
        int maxt = 2 * maxL;
        m_E[dim] = zeros(maxiA, maxiB, maxt);

        double p = m_totalExponent;
        double mu = m_reducedExponent;
        double dim_AB = m_relativeSeparation(dim);
        double dim_PA = m_centerOfMassDiffA(dim);
        double dim_PB = m_centerOfMassDiffB(dim);

        //    vector<string> performed;

        // First row is special
        m_E[dim](0,0,0) = exp(-mu * dim_AB * dim_AB);
//        cout << "KAB=" << m_E[dim](0,0,0) << endl;
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

double GaussianTypeOrbitalIntegrator::coreColoumbIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double result = 0;
    const cube &E_x = m_E[0];
    const cube &E_y = m_E[1];
    const cube &E_z = m_E[2];
    const vector<cube> &R = m_R;
    int tMax = iA + iB;
    int uMax = jA + jB;
    int vMax = kA + kB;
    for(int t = 0; t < tMax + 1; t++) {
        for(int u = 0; u < uMax + 1; u++) {
            for(int v = 0; v < vMax + 1; v++) {
                result += E_x(iA, iB, t) * E_y(jA, jB, u) * E_z(kA, kB, v) * R[0](t,u,v);
            }
        }
    }
    return result;
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
    const rowvec &C = m_corePositionC;
    const rowvec &P = m_centerOfMass;
    const rowvec &PC = P - C;
    cout << "PC " << PC << endl;
//    const rowvec &PA = m_centerOfMassDiffA;
//    const rowvec &PB = m_centerOfMassDiffB;
    double boysArg = p * dot(PC, PC);
//    int lMax = 2;
    int tMax = m_angularMomentumAMax + m_angularMomentumBMax;
    int tuvSumMax = 3 * tMax;
    int nMax = 3 * tMax;
    BoysFunctionIntermediate boysFunctionIntermediate(nMax, 1000, 0, 50, 1e5);
    BoysFunction boysFunction(boysArg, nMax + 1, &boysFunctionIntermediate);
//    vector<vector<int>> calculatedIndices;
    m_R.clear();
    // Initialize and allocate the cubes for R
    for(int n = 0; n < nMax + 1; n++) {
        m_R.push_back(zeros(nMax, nMax, nMax));
    }
    // Calculate R0_tuv
    for(int n = 0; n < nMax + 1; n++) {
        m_R[n](0,0,0) = pow(-2*p, n) * boysFunction.result(n);
//        cout << m_R[n](0,0,0) << endl;
//        vector<int> indices = {n,0,0,0};
//        calculatedIndices.push_back(indices);
    }
    for(int tuvSum = 1; tuvSum < tuvSumMax + 1; tuvSum++) {
        cout << "All elements for which t + u + v = " << tuvSum << endl;
        for(int n = 0; n < nMax - tuvSum; n++) {
            for(int t = 0; t < tMax + 1; t++) {
                for(int u = 0; u < tMax + 1; u++) {
                    for(int v = 0; v < tMax + 1; v++) {
                        if(t + u + v != tuvSum || t + u + v == 0) {
                            continue;
                        }
//                        if(n == 0 && t == 0 && u == 1 && v == 0) {
//                            cout << "woot" << endl;
//                        }
                        int largestElement = max(t,max(u,v));
                        // Needed indices
                        int t2 = t;
                        int u2 = u;
                        int v2 = v;
                        int t3 = t;
                        int u3 = u;
                        int v3 = v;
                        // Resulting factors
                        int factor2 = t;
                        double factor3 = PC(0);
                        // Choose indices based on method (slide 34)
                        if(largestElement == t) {
//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "t biggest" << endl;
//                            }
                            t2 = t - 2;
                            t3 = t - 1;
                            factor2 = t - 1;
                            factor3 = PC(0);
                        } else if(largestElement == u) {
//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "u biggest" << endl;
//                            }
                            u2 = u - 2;
                            u3 = u - 1;
                            factor2 = u - 1;
                            factor3 = PC(1);
                        } else {
                            v2 = v - 2;
                            v3 = v - 1;
                            factor2 = v - 1;
                            factor3 = PC(2);
                        }
                        double currentValue = 0;
//                        vector<int> currentIndices = {n, tRec1, uRec1, vRec1};
                        if(t2 >= 0 && u2 >= 0 && v2 >= 0) {

//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "Adding factor2 " << endl;
//                            }
                            currentValue += factor2 * m_R[n+1](t2, u2, v2);
//                            vector<int> neededIndices = {n+1,tRec2,uRec2,vRec2};
//                            if(!hasIndices(neededIndices, calculatedIndices)) {
//                                cout << "Target: " <<  currentIndices[0] << " " << currentIndices[1] << " " << currentIndices[2] << " " << currentIndices[3] << endl;
//                                cout << "Needs: " <<  neededIndices[0] << " " << neededIndices[1] << " " << neededIndices[2] << " " << neededIndices[3] << endl;
//                                break;
//                            }
                        }
                        if(t3 >= 0 && u3 >= 0 && v3 >= 0) {
//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "Adding factor3 " << endl;
//                                cout << n + 1 << " " << t3 << " " << u3 << " " << v3 << endl;
//                                cout << factor3 << endl;
//                                cout << PC(1) << endl;
//                            }
                            currentValue += factor3 * m_R[n+1](t3, u3, v3);
//                            vector<int> neededIndices2 = {n+1,tRec3,uRec3,vRec3};
//                            if(!hasIndices(neededIndices2, calculatedIndices)) {
//                                cout << "Target: " <<  currentIndices[0] << " " << currentIndices[1] << " " << currentIndices[2] << " " << currentIndices[3] << endl;
//                                cout << "Needs: " <<  neededIndices2[0] << " " << neededIndices2[1] << " " << neededIndices2[2] << " " << neededIndices2[3] << endl;
//                                break;
//                            }
                        }

//                        if(n == 0 && t == 1 && u == 0 && v == 0) {
//                            cout << "currentValue " << currentValue << endl;
//                        }
                        m_R[n](t, u, v) = currentValue;
//                        cout << "Calculated: " <<  currentIndices[0] << " " << currentIndices[1] << " " << currentIndices[2] << " " << currentIndices[3] << endl;
//                        calculatedIndices.push_back(currentIndices);
                    }
                }
            }
        }
    }
//    cout << banana.find(elementa) << endl;

//    double F0 = boysFunction.result(0);
//    double F1 = boysFunction.result(1);
//    double F2 = boysFunction.result(2);
//    double R0000 = F0;
//    double R1000 = -2 * p * F1;
//    double R2000 = (-2 * p)*(-2*p) * F2;
//    cout << "R0000 = " << R0000 << endl;
//    cout << "R1000 = " << R1000 << endl;
//    cout << "R2000 = " << R2000 << endl;
//    cube *E = m_E;
//    cout << "Result = " << 2 * M_PI / p * R0000 * E[0](0,0,0) * E[1](0,0,0) * E[2](0,0,0) << endl;
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
