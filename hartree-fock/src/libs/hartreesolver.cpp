#include "hartreesolver.h"

#include <electronsystems/electronsystem.h>

#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;

HartreeSolver::HartreeSolver(ElectronSystem *basisFunction) :
    m_basisFunction(basisFunction)
{
    cout << setprecision(20);
    allocateQMemory();
    reset();
}

HartreeSolver::~HartreeSolver()
{
    cleanUpQMemory();
}

void HartreeSolver::reset() {
    setuph();
    setupS();
    setupQ();
    resetC();
}

void HartreeSolver::setuph() {
    ElectronSystem* f = m_basisFunction;
    uint n = f->nOrbitals();
    h.reset();
    h = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            h(p,q) = f->kineticIntegral(p,q) + f->nuclearAttractionIntegral(p,q);
        }
    }
}

void HartreeSolver::setupS() {
    ElectronSystem* f = m_basisFunction;
    uint n = f->nOrbitals();
    S.reset();
    S = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            S(p,q) = f->overlapIntegral(p, q);
        }
    }
}

void HartreeSolver::allocateQMemory() {
    if(!isQAllocated) {
        ElectronSystem* f = m_basisFunction;
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

void HartreeSolver::cleanUpQMemory() {
    if(isQAllocated) {
        ElectronSystem* f = m_basisFunction;
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

void HartreeSolver::setupQ() {
    ElectronSystem* f = m_basisFunction;
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

void HartreeSolver::resetC() {
    ElectronSystem* f = m_basisFunction;
    uint n = f->nOrbitals();
    C.reset();
    C = ones(n);
}

void HartreeSolver::advance() {
    normalizeCwithRegardsToS();
    setupF();

    vec s;
    mat U;
    eig_sym(s, U, S);

    mat V = U*diagmat(1.0/sqrt(s));

    F = V.t() * F * V;

    vec eps;
    mat Cmat;
    eig_sym(eps, Cmat, F);

    C = V*Cmat.col(0);
    normalizeCwithRegardsToS();

    double energy = 0;

    ElectronSystem* f = m_basisFunction;
    uint n = f->nOrbitals();

    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            energy += 2 * C(p) * C(q) * h(p,q);
        }
    }

    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    energy += Q[p][r][q][s] * C(p) * C(q) * C(r) * C(s);
                }
            }
        }
    }
    m_energy = energy;
}

void HartreeSolver::normalizeCwithRegardsToS(){
    double factor = 0.0;

    for(uint i= 0; i < C.n_elem; i++){
        for(uint j= 0; j < C.n_elem; j++){
            factor += C(i)*S(i,j)*C(j);
        }
    }

    C = C/sqrt(factor);
}

void HartreeSolver::setupF() {
    ElectronSystem* f = m_basisFunction;
    uint n = f->nOrbitals();
    F = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    F(p,q) += Q[p][r][q][s] * C(r) * C(s);
                }
            }
        }
    }
    F = F + h;
}
