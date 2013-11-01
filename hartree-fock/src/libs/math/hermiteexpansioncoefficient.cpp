#include "hermiteexpansioncoefficient.h"

HermiteExpansionCoefficient::HermiteExpansionCoefficient() :
    HermiteExpansionCoefficient(0,rowvec(3),0,rowvec(3),0,false)
{
}

HermiteExpansionCoefficient::HermiteExpansionCoefficient(double a, rowvec A, double b, rowvec B, int angularMomentumMax, bool setupImmediately) :
    m_a(a),
    m_b(b),
    m_A(A),
    m_B(B),
    m_angularMomentumMax(angularMomentumMax)
{
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
    int maxL = m_angularMomentumMax;
    for(int dim = 0; dim < 3; dim++) {
        int maxiA = maxL + 3;
        int maxiB = maxL + 3;
        // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
        // needed values available for the algorithm.
        int maxt = 2 * maxL;
        m_E[dim] = zeros(maxiA, maxiB, maxt);

        double a = m_a;
        double b = m_b;
        double p = a + b;
        double mu = a * b / (a + b);
        rowvec AB = m_A - m_B;
        rowvec P = (a * m_A + b * m_B) / p;
        rowvec PA = P - m_A;
        rowvec PB = P - m_B;
        double dim_AB = AB(dim);
        double dim_PA = PA(dim);
        double dim_PB = PB(dim);

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
