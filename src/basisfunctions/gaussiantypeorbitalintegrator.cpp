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
    m_Ex = zeros(maxL + 1, maxL + 1, maxL + 1);

    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;

    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    double mu = a * b / (a + b);
    double X_AB = A(0) - B(0);
    double P_x = (a * A(0) + b * B(0)) / p;
    double X_PA = P_x - A(0);
    double X_PB = P_x - B(0);

    cout << "X_PA = " << X_PA << endl;
    cout << "X_PB = " << X_PB << endl;

    vector<string> performed;

    // First row is special
    m_Ex(0,0,0) = exp(-mu * X_AB * X_AB);
    performed.push_back("000");
    for(int j = 0; j < maxL + 1; j++) {
        for(int t = 0; t < maxL + 1; t++) {
            int i = 0;
            if(i == 0 && j == 0 && t == 0) {
                continue;
            }
            // p = previous, n = next
            // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
            int jp = j - 1;
            int tp = t - 1;
            int tn = t + 1;
            double E_i_jp_tp;
            if(tp < 0 || tp > (i + jp) || jp < 0) {
                E_i_jp_tp = 0;
            } else {
                stringstream idToFind;
                idToFind << i << jp << tp;
                if(find(performed.begin(), performed.end(), idToFind.str()) == performed.end()) {
                    cout << "Could not find E_i_jp_tp " << idToFind.str() << endl;
                }
                E_i_jp_tp = m_Ex(i, jp, tp);
            }
            double E_i_jp_t;
            if(t > (i + jp) || jp < 0) {
                E_i_jp_t = 0;
            } else {
                stringstream idToFind;
                idToFind << i << jp << t;
                if(find(performed.begin(), performed.end(), idToFind.str()) == performed.end()) {
                    cout << "Could not find E_i_jp_t" << idToFind.str() << endl;
                }
                E_i_jp_t = m_Ex(i, jp, t);
            }
            double E_i_jp_tn;
            if(tn > maxL || tn > (i + jp) || jp < 0) {
                E_i_jp_tn = 0;
            } else {
                stringstream idToFind;
                idToFind << i << jp << tn;
                if(find(performed.begin(), performed.end(), idToFind.str()) == performed.end()) {
                    cout << "Could not find E_i_jp_tn" << idToFind.str() << endl;
                }
                E_i_jp_tn = m_Ex(i, jp, tn);
            }
            m_Ex(i,j,t) = 1 / (2*p) * E_i_jp_tp + X_PB * E_i_jp_t +  (t + 1)*E_i_jp_tn;

            cout << i << j << t << " = " << m_Ex(i,j,t) << endl;

            stringstream id;
            id << i << j << t;
            performed.push_back(id.str());
        }
    }
    cout << "Done first row" << endl;
    cout << m_Ex;
    for(int i = 1; i < maxL + 1; i++) {
        for(int j = 0; j < maxL + 1; j++) {
            for(int t = 0; t < maxL + 1; t++) {
                // p = previous, n = next
                // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
                int ip = i - 1;
                int tp = t - 1;
                int tn = t + 1;
                double E_ip_j_tp;
                if(tp < 0 || tp > (ip + j) || ip < 0) {
                    E_ip_j_tp = 0;
                } else {
                    stringstream idToFind;
                    idToFind << ip << j << tp;
                    if(find(performed.begin(), performed.end(), idToFind.str()) == performed.end()) {
                        cout << "Could not find E_ip_j_tp" << idToFind.str() << endl;
                    }
                    E_ip_j_tp = m_Ex(ip, j, tp);
                }
                double E_ip_j_t;
                if(t > (ip + j) || ip < 0) {
                    E_ip_j_t = 0;
                } else {
                    stringstream idToFind;
                    idToFind << ip << j << t;
                    if(find(performed.begin(), performed.end(), idToFind.str()) == performed.end()) {
                        cout << "Could not find E_ip_j_t" << idToFind.str() << endl;
                    }
                    E_ip_j_t = m_Ex(ip, j, t);
                }
                double E_ip_j_tn;
                if(tn > maxL || tn > (ip + j) || ip < 0) {
                    E_ip_j_tn = 0;
                } else {
                    stringstream idToFind;
                    idToFind << ip << j << tn;
                    if(find(performed.begin(), performed.end(), idToFind.str()) == performed.end()) {
                        cout << "Could not find E_ip_j_tn" << idToFind.str() << endl;
                    }
                    E_ip_j_tn = m_Ex(ip, j, tn);
                }
                m_Ex(i,j,t) = 1 / (2*p) * E_ip_j_tp + X_PA * E_ip_j_t +  (t + 1)*E_ip_j_tn;
                cout << i << j << t << " = " << m_Ex(i,j,t) << endl;
                stringstream id;
                id << i << j << t;
                performed.push_back(id.str());
            }
        }
    }

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
