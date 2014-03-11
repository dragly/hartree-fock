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
    m_uncoupledMatrix.reset();
    m_uncoupledMatrix = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            m_uncoupledMatrix(p,q) = f->uncoupledIntegral(p,q);
        }
    }
}

void HartreeFockSolver::setupS() {
    ElectronSystem* f = m_electronSystem;
    uint nOrbitals = f->nBasisFunctions();
    m_overlapMatrix.reset();
    m_overlapMatrix = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            m_overlapMatrix(p,q) = f->overlapIntegral(p, q);
        }
    }
}

void HartreeFockSolver::allocateQMemory() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    m_coupledMatrix.set_size(n,n);
    for(int i = 0; i < int(n); i++) {
        for(int j = 0; j < int(n); j++) {
            m_coupledMatrix(i,j) = zeros(n,n);
        }
    }
}

void HartreeFockSolver::cleanUpQMemory() {
    m_coupledMatrix.reset();
}

void HartreeFockSolver::setupQ() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    for(uint p = 0; p < n; p++) {
        for(uint r = 0; r < n; r++) {
            for(uint q = 0; q < n; q++) {
                for(uint s = 0; s < n; s++) {
                    m_coupledMatrix(p,r)(q,s) = f->coupledIntegral(p, q, r, s); // NOTE: Indexes changed on purpose in element and call
                }
            }
        }
    }
}

void HartreeFockSolver::resetC() {
    ElectronSystem* f = m_electronSystem;
    m_coefficientMatrix.reset();
    m_coefficientMatrix = zeros(f->nBasisFunctions(), f->nParticles() / 2);
}

void HartreeFockSolver::advance() {
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;
    nk = max(int(nk), 2); // TODO Fix this hack for the Hydrogen atom
    setupF();

    vec s;
    mat U;
    eig_sym(s, U, m_overlapMatrix);

    mat V = U*diagmat(1.0/sqrt(s));

    mat &F = m_fockMatrix;
    F = V.t() * F * V;

    vec eps;
    mat Cmat;
    eig_sym(eps, Cmat, m_fockMatrix);


    mat &C = m_coefficientMatrix;
    C = V*Cmat.submat(0, 0, no - 1, nk - 1);
    normalizeCwithRegardsToS();

    setupP();

    double energy = 0;

    mat& P = m_densityMatrix;
    mat& h = m_uncoupledMatrix;
    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            energy += P(p,q) * h(p,q);
        }
    }

    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            for(uint r = 0; r < no; r++) {
                for(uint s = 0; s < no; s++) {
                    double Qtilde = coupledMatrixTilde(p, q, r, s);
                    energy += 0.25 * Qtilde * P(p,q) * P(s,r);
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

    const mat& S = m_overlapMatrix;
    mat& C = m_coefficientMatrix;

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

double HartreeFockSolver::coupledMatrixTilde(int p, int q, int r, int s) { // TODO: Find a better name for tilde
    field<mat>& Q = m_coupledMatrix;
    return 2 * Q(p,r)(q,s) - Q(p,r)(s,q);
}

void HartreeFockSolver::setupF() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    mat &F = m_fockMatrix;
    mat &P = m_densityMatrix;
    mat &h = m_uncoupledMatrix;
    F = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            F(p,q) = h(p,q);
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    double Qtilde = coupledMatrixTilde(p, q, r, s);
                    F(p,q) += 0.5 * Qtilde * P(s,r);
                }
            }
        }
    }
}

void HartreeFockSolver::setupP() {
    mat &P = m_densityMatrix;
    mat &C = m_coefficientMatrix;
    P = 2 * C * C.t();
}
