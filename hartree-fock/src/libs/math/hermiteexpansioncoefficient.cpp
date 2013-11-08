#include "hermiteexpansioncoefficient.h"

HermiteExpansionCoefficient::HermiteExpansionCoefficient(int dimension) :
    m_a(0),
    m_b(0),
    m_A(zeros<rowvec>(3)),
    m_B(zeros<rowvec>(3))
{
    reset(dimension);
}

HermiteExpansionCoefficient::HermiteExpansionCoefficient(double a, double b, rowvec A, rowvec B, int angularMomentumMax, bool setupImmediately) :
    HermiteExpansionCoefficient(angularMomentumMax)
{
    set(a,b,A,B,setupImmediately);
}

void HermiteExpansionCoefficient::reset(int dimension)
{
    m_dimension = dimension;
    int maxL = m_dimension;
    int maxiA = m_dimension;
    int maxiB = m_dimension;
    // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
    // needed values available for the algorithm.
    int maxt = 2 * maxL;
    for(int dim = 0; dim < 3; dim++) {
        cube& dim_E = m_E[dim];
        dim_E = zeros(maxiA, maxiB, maxt + 1);
    }
}

void HermiteExpansionCoefficient::set(double a, double b, rowvec A, rowvec B, bool setupImmediately)
{
    m_a = a;
    m_b = b;
    m_A = A;
    m_B = B;
    if(setupImmediately) {
        setupE();
    }
}

bool HermiteExpansionCoefficient::checkIndexCombinationForE(int iA, int iB, int t) {
    if(t < 0 || t > (iA + iB) || iA < 0 || iB < 0) {
        return false;
    } else {
        return true;
    }
}

void HermiteExpansionCoefficient::setupE() {
    int maxL = m_dimension;
    int maxiA = m_dimension;
    int maxiB = m_dimension;
    // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
    // needed values available for the algorithm.
    int maxt = 2 * maxL;
    double a = m_a;
    double b = m_b;
    double p = a + b;
    double mu = a * b / (a + b);
    m_AB = m_A - m_B;
    m_P = (a * m_A + b * m_B) / p;
    m_PA = m_P - m_A;
    m_PB = m_P - m_B;
    for(int dim = 0; dim < 3; dim++) {
        cube& dim_E = m_E[dim];
//        dim_E = zeros(maxiA, maxiB, maxt + 1);
        double dim_AB = m_AB(dim);
        double dim_PA = m_PA(dim);
        double dim_PB = m_PB(dim);

        //    vector<string> performed;

        // First row is special
        dim_E(0,0,0) = exp(-mu * dim_AB * dim_AB);
//        cout << "KAB=" << dim_E(0,0,0) << endl;
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
                    E_iA_iBp_tp = dim_E(iA, iBp, tp);
                }
                double E_iA_iBp_t = 0;
                if(checkIndexCombinationForE(iA, iBp, t)) {
                    E_iA_iBp_t = dim_E(iA, iBp, t);
                }
                double E_iA_iBp_tn = 0;
                if(checkIndexCombinationForE(iA, iBp, tn)) {
                    E_iA_iBp_tn = dim_E(iA, iBp, tn);
                }
                dim_E(iA,iB,t) = 1 / (2*p) * E_iA_iBp_tp + dim_PB * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
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
                        E_iAp_iB_tp = dim_E(iAp, iB, tp);
                    }
                    double E_iAp_iB_t = 0;
                    if(checkIndexCombinationForE(iAp, iB, t)) {
                        E_iAp_iB_t = dim_E(iAp, iB, t);
                    }
                    double E_iAp_iB_tn = 0;
                    if(checkIndexCombinationForE(iAp, iB, tn)) {
                        E_iAp_iB_tn = dim_E(iAp, iB, tn);
                    }
                    dim_E(iA,iB,t) = 1 / (2*p) * E_iAp_iB_tp + dim_PA * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
                }
            }
        }
//        cout << dim_E;
    }
}
