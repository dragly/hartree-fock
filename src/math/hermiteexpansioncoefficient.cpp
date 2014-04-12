#include "hermiteexpansioncoefficient.h"

HermiteExpansionCoefficient::HermiteExpansionCoefficient(int dimensionMax) :
    m_a(0),
    m_b(0),
    m_A(zeros<rowvec>(3)),
    m_B(zeros<rowvec>(3))
{
    reset(dimensionMax);
}

void HermiteExpansionCoefficient::reset(int dimensionMax)
{
    m_dimensionMax = dimensionMax;
    int maxL = m_dimensionMax;
    int maxiA = m_dimensionMax;
    int maxiB = m_dimensionMax;
    // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
    // needed values available for the algorithm.
    int maxt = 2 * maxL;
    for(int dim = 0; dim < 3; dim++) {
        cube& dim_E = m_E[dim];
        dim_E = zeros(maxiA, maxiB, maxt + 1);
    }
}

void HermiteExpansionCoefficient::set(double a, double b, rowvec A, rowvec B,
                                      int iA, int iB, int jA, int jB, int kA, int kB)
{
    m_a = a;
    m_b = b;
    m_A = A;
    m_B = B;
    setupE(iA, iB, jA, jB, kA, kB);
}

bool HermiteExpansionCoefficient::checkIndexCombinationForE(int iA, int iB, int t) {
    if(t < 0 || t > (iA + iB) || iA < 0 || iB < 0) {
        return false;
    } else {
        return true;
    }
}

void HermiteExpansionCoefficient::setupE(int iA, int iB, int jA, int jB, int kA, int kB) {
    int maxiAs[3];
    maxiAs[0] = iA + 1;
    maxiAs[1] = jA + 1;
    maxiAs[2] = kA + 1;

    int maxiBs[3];
    maxiBs[0] = iB + 1;
    maxiBs[1] = jB + 1;
    maxiBs[2] = kB + 1;

    // Since we are only going to need E_i_j_0 for the integrals, we cap t to i and j to only have the
    // needed values available for the algorithm.
    int maxts[3];
    maxts[0] = iA + iB + 1;
    maxts[1] = jA + jB + 1;
    maxts[2] = kA + kB + 1;

    double a = m_a;
    double b = m_b;
    double p = a + b;
    double mu = a * b / (a + b);
    m_AB = m_A - m_B;
    m_P = (a * m_A + b * m_B) / p;
    m_PA = m_P - m_A;
    m_PB = m_P - m_B;
    for(int dim = 0; dim < 3; dim++) {
        cube& E = m_E[dim];
        double AB = m_AB(dim);
        double PA = m_PA(dim);
        double PB = m_PB(dim);

        int maxiA = maxiAs[dim];
        int maxiB = maxiBs[dim];
        int maxt = maxts[dim];

        // First row is special
        E(0,0,0) = exp(-mu * AB * AB);

        // Build the first row, iA = 0
        for(int iB = 0; iB < maxiB; iB++) {
            for(int t = 0; t < maxt; t++) {
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
                    E_iA_iBp_tp = E(iA, iBp, tp);
                }
                double E_iA_iBp_t = 0;
                if(checkIndexCombinationForE(iA, iBp, t)) {
                    E_iA_iBp_t = E(iA, iBp, t);
                }
                double E_iA_iBp_tn = 0;
                if(checkIndexCombinationForE(iA, iBp, tn)) {
                    E_iA_iBp_tn = E(iA, iBp, tn);
                }
                E(iA,iB,t) = 1 / (2*p) * E_iA_iBp_tp + PB * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
            }
        }

        // Iterate over all rows for iA >= 1
        for(int iA = 1; iA < maxiA; iA++) {
            for(int iB = 0; iB < maxiB; iB++) {
                for(int t = 0; t < maxt; t++) {
                    // p = previous, n = next
                    // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
                    int iAp = iA - 1;
                    int tp = t - 1;
                    int tn = t + 1;
                    double E_iAp_iB_tp = 0;
                    if(checkIndexCombinationForE(iAp, iB, tp)) {
                        E_iAp_iB_tp = E(iAp, iB, tp);
                    }
                    double E_iAp_iB_t = 0;
                    if(checkIndexCombinationForE(iAp, iB, t)) {
                        E_iAp_iB_t = E(iAp, iB, t);
                    }
                    double E_iAp_iB_tn = 0;
                    if(checkIndexCombinationForE(iAp, iB, tn)) {
                        E_iAp_iB_tn = E(iAp, iB, tn);
                    }
                    E(iA,iB,t) = 1 / (2*p) * E_iAp_iB_tp + PA * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
                }
            }
        }
    }
}
