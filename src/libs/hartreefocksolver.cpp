#include "hartreefocksolver.h"

#include <electronsystems/electronsystem.h>

#include <armadillo>
#include <iomanip>
#include <cmath>

using namespace arma;
using namespace std;

HartreeFockSolver::HartreeFockSolver(ElectronSystem *basisFunction) :
    m_electronSystem(basisFunction)
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
    //    normalizeCwithRegardsToS();
    setupP();
}

void HartreeFockSolver::setuph() {
    ElectronSystem* f = m_electronSystem;
    uint nOrbitals = f->nBasisFunctions();
    h.reset();
    h = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            h(p,q) = f->uncoupledIntegral(p,q);
        }
    }
}

void HartreeFockSolver::setupS() {
    ElectronSystem* f = m_electronSystem;
    uint nOrbitals = f->nBasisFunctions();
    S.reset();
    S = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            S(p,q) = f->overlapIntegral(p, q);
        }
    }
}

void HartreeFockSolver::allocateQMemory() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    Q.set_size(n,n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            Q(i,j) = zeros(n,n);
        }
    }
}

void HartreeFockSolver::cleanUpQMemory() {
    Q.reset();
}

void HartreeFockSolver::setupQ() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    for(uint p = 0; p < n; p++) {
        for(uint r = 0; r < n; r++) {
            for(uint q = 0; q < n; q++) {
                for(uint s = 0; s < n; s++) {
                    Q(p,r)(q,s) = f->coupledIntegral(p, q, r, s); // NOTE: Indexes changed on purpose in element and call
                }
            }
        }
    }
}

void HartreeFockSolver::resetC() {
    ElectronSystem* f = m_electronSystem;
    C.reset();
    C = zeros(f->nBasisFunctions(), f->nParticles() / 2);
}

void HartreeFockSolver::advance() {
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;
    nk = max(int(nk), 2); // TODO Fix this hack for the Hydrogen atom
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
                    energy += 0.25 * Qtilde(p, q, r, s) * P(p,q) * P(s,r);
                }
            }
        }
    }
    energy += m_electronSystem->additionalEnergyTerms();
    m_energy = energy;
}

void HartreeFockSolver::normalizeCwithRegardsToS(){
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
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
    return 2 * Q(p,r)(q,s) - Q(p,r)(s,q);
}

void HartreeFockSolver::setupF() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    F = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            F(p,q) = h(p,q);
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    F(p,q) += 0.5 * Qtilde(p, q, r, s) * P(s,r);
                }
            }
        }
    }
}

void HartreeFockSolver::setupP() {
    P = 2 * C * C.t();
}
