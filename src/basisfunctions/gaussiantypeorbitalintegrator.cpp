#include "gaussiantypeorbitalintegrator.h"

using namespace std;
using namespace arma;

GaussianTypeOrbitalIntegrator::GaussianTypeOrbitalIntegrator() :
    m_exponentA(0),
    m_exponentB(0),
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
    setupE();
}

void GaussianTypeOrbitalIntegrator::setupE() {
    int maxL = max(maxAngularMomentumA(), maxAngularMomentumB());
    for(int dim = 0; dim < 3; dim++) {
        int maxi = maxL + 3;
        int maxj = maxL + 3;
        // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
        // needed values available for the algorithm.
        int maxt = maxL + 3;
        m_E[dim] = zeros(maxi, maxj, maxt);

        const rowvec &A = m_corePositionA;
        const rowvec &B = m_corePositionB;

        double a = m_exponentA;
        double b = m_exponentB;
        double p = a + b;
        double mu = a * b / (a + b);
        double dim_AB = A(dim) - B(dim);
        double P_dim = (a * A(dim) + b * B(dim)) / p;
        double dim_PA = P_dim - A(dim);
        double dim_PB = P_dim - B(dim);

        cout << "dim_PA = " << dim_PA << endl;
        cout << "dim_PB = " << dim_PB << endl;

        //    vector<string> performed;

        // First row is special
        m_E[dim](0,0,0) = exp(-mu * dim_AB * dim_AB);
        //    performed.push_back("000");
        for(int j = 0; j < maxj; j++) {
            for(int t = 0; t < maxt - 1; t++) {
                int i = 0;
                if(i == 0 && j == 0 && t == 0) {
                    continue;
                }
                // p = previous, n = next
                // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
                int jp = j - 1;
                int tp = t - 1;
                int tn = t + 1;
                double E_i_jp_tp = 0;
                if(!(tp < 0 || tp > (i + jp) || jp < 0)) {
                    E_i_jp_tp = m_E[dim](i, jp, tp);
                }
                double E_i_jp_t = 0;
                if(!(t > (i + jp) || jp < 0)) {
                    E_i_jp_t = m_E[dim](i, jp, t);
                }
                double E_i_jp_tn = 0;
                if(!(tn > (i + jp) || jp < 0)) {
                    E_i_jp_tn = m_E[dim](i, jp, tn);
                }
                m_E[dim](i,j,t) = 1 / (2*p) * E_i_jp_tp + dim_PB * E_i_jp_t +  (t + 1)*E_i_jp_tn;
            }
        }
        for(int i = 1; i < maxi; i++) {
            for(int j = 0; j < maxj; j++) {
                for(int t = 0; t < maxt - 1; t++) {
                    // p = previous, n = next
                    // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
                    int ip = i - 1;
                    int tp = t - 1;
                    int tn = t + 1;
                    double E_ip_j_tp = 0;
                    if(!(tp < 0 || tp > (ip + j) || ip < 0)) {
                        E_ip_j_tp = m_E[dim](ip, j, tp);
                    }
                    double E_ip_j_t = 0;
                    if(!(t > (ip + j) || ip < 0)) {
                        E_ip_j_t = m_E[dim](ip, j, t);
                    }
                    double E_ip_j_tn = 0;
                    if(!(tn > (ip + j) || ip < 0)) {
                        E_ip_j_tn = m_E[dim](ip, j, tn);
                    }
                    m_E[dim](i,j,t) = 1 / (2*p) * E_ip_j_tp + dim_PA * E_ip_j_t +  (t + 1)*E_ip_j_tn;
                }
            }
        }
        cout << m_E[dim];
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
    cube* E = m_E;
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
//    return E[0](iA,iB,0) * E[1](jA,jB,0) * E[2](kA,kB,0) * pow((M_PI / p), 3./2.);
}

double GaussianTypeOrbitalIntegrator::kineticIntegral(int dim, int iA, int iB) {
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    const cube &E_dim = m_E[dim];
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
