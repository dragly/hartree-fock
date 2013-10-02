#include "hartreefocksolver.h"

#include <src/basisfunctions/basisfunction.h>

#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;

HartreeFockSolver::HartreeFockSolver(BasisFunction *basisFunction) :
    m_basisFunction(basisFunction)
{
    cout << setprecision(20);
    allocateQMemory();
    reset();
}

HartreeFockSolver::~HartreeFockSolver()
{
    cleanUpQMemory();
}

void HartreeFockSolver::reset() {
    setuph();
    setupS();
    setupQ();
    resetC();
}

void HartreeFockSolver::setuph() {
    BasisFunction* f = m_basisFunction;
    uint nOrbitals = f->nOrbitals();
    h.reset();
    h = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            h(p,q) = f->kineticIntegral(p,q) + f->nuclearAttractionIntegral(p,q);
        }
    }
}

void HartreeFockSolver::setupS() {
    BasisFunction* f = m_basisFunction;
    uint nOrbitals = f->nOrbitals();
    S.reset();
    S = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            S(p,q) = f->overlapIntegral(p, q);
        }
    }
}

void HartreeFockSolver::allocateQMemory() {
    if(!isQAllocated) {
        BasisFunction* f = m_basisFunction;
        uint n = f->nOrbitals();
        QData = new double[n*n*n*n];
        Q = new double***[n];
        for(uint p = 0; p < n; p++) {
            Q[p] = new double**[n];
            for(uint r = 0; r < n; r++) {
                Q[p][r] = new double *[n];
                for(uint q = 0; q < n; q++) {
                    Q[p][r][q] = &QData[n*n*n*p + n*n*r + n*q];
                    for(uint s = 0; s < n; s++) {
                        Q[p][r][q][s] = 0;
                    }
                }
            }
        }
        isQAllocated = true;
    }
}

void HartreeFockSolver::cleanUpQMemory() {
    if(isQAllocated) {
        BasisFunction* f = m_basisFunction;
        uint n = f->nOrbitals();
        for (uint i = 0; i < n; ++i) {
            for (uint j = 0; j < n; ++j){
                delete [] Q[i][j];
            }
            delete [] Q[i];
        }
        delete [] Q;
        delete []QData;
    }
}

void HartreeFockSolver::setupQ() {
    BasisFunction* f = m_basisFunction;
    uint n = f->nOrbitals();
    for(uint p = 0; p < n; p++) {
        for(uint r = 0; r < n; r++) {
            for(uint q = 0; q < n; q++) {
                for(uint s = 0; s < n; s++) {
                    Q[p][r][q][s] = f->electronInteractionIntegral(p, r, q, s);
                }
            }
        }
    }
}

void HartreeFockSolver::resetC() {
    BasisFunction* f = m_basisFunction;
    C.reset();
    C = ones(f->nOrbitals(), f->nParticles() / 2);
}

void HartreeFockSolver::advance() {
    BasisFunction* f = m_basisFunction;
    uint no = f->nOrbitals();
    uint nk = f->nParticles() / 2;

    normalizeCwithRegardsToS();
    setupP();
    setupF();

    vec s;
    mat U;
    eig_sym(s, U, S);

    mat V = U*diagmat(1.0/sqrt(s));

    F = V.t() * F * V;

    vec eps;
    mat Cmat;
    eig_sym(eps, Cmat, F);


    C = V*Cmat.submat(0, 0, no - 1, nk - 1);
    normalizeCwithRegardsToS();

    setupP();

    double energy = 0;

    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            energy += P(p,q) * h(p,q);
        }
    }

    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            for(uint r = 0; r < no; r++) {
                for(uint s = 0; s < no; s++) {
                    energy += 0.25 * Qtilde(p, q, r, s) * P(p,q) * P(r,s);
                }
            }
        }
    }
    m_energy = energy;
}

void HartreeFockSolver::normalizeCwithRegardsToS(){

    BasisFunction* f = m_basisFunction;
    uint no = f->nOrbitals();
    uint nk = f->nParticles() / 2;

    for(uint k = 0; k < nk; k++) {
        double factor = 0.0;
        for(uint p = 0; p < no; p++){
            for(uint q = 0; q < no; q++){
                factor += C(p,k) * S(p,q) * C(q,k);
            }
        }
        C.col(k) = C.col(k) / sqrt(factor);
    }
}

double HartreeFockSolver::Qtilde(int p, int q, int r, int s) {
    return 2 * Q[p][r][q][s] - Q[p][r][s][q];
}

void HartreeFockSolver::setupF() {
    BasisFunction* f = m_basisFunction;
    uint n = f->nOrbitals();
    F = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            F(p,q) = h(p,q);
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    F(p,q) += 0.5 * Qtilde(p, q, r, s) * P(r,s);
                }
            }
        }
    }
}

void HartreeFockSolver::setupP() {
    BasisFunction* f = m_basisFunction;
    uint no = f->nOrbitals();
    uint nk = f->nParticles() / 2;
    P = zeros(no,no);
    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            for(uint k = 0; k < nk; k++) {
                P(p,q) += 2 * C(p,k) * C(q,k);
            }
        }
    }
}
