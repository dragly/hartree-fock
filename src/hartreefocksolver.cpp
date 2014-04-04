#include "hartreefocksolver.h"

#include <electronsystems/electronsystem.h>

#include <armadillo>
#include <iomanip>
#include <cmath>

using namespace arma;
using namespace std;

HartreeFockSolver::HartreeFockSolver(ElectronSystem *basisFunction) :
    m_electronSystem(basisFunction),
    m_convergenceTreshold(1e-8),
    m_nIterationsMax(1e3)
{
    cout << setprecision(20);
    allocateCoupledMatrix();
    reset();
}

HartreeFockSolver::~HartreeFockSolver()
{
    cleanUpCoupledMatrix();
}

void HartreeFockSolver::reset() {
    setupUncoupledMatrix();
    setupOverlapMatrix();
    setupCoupledMatrix();
    resetCoefficientMatrix();
//    normalizeCoefficientMatrix();
    setupDensityMatrix();
}

void HartreeFockSolver::resetCoefficientMatrix() {
    ElectronSystem* f = m_electronSystem;
    m_coefficientMatrix.reset();
    m_coefficientMatrix = zeros(f->nBasisFunctions(), f->nParticles() / 2);
}

void HartreeFockSolver::allocateCoupledMatrix() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    m_coupledMatrix.set_size(n,n);
    for(int i = 0; i < int(n); i++) {
        for(int j = 0; j < int(n); j++) {
            m_coupledMatrix(i,j) = zeros(n,n);
        }
    }
}

void HartreeFockSolver::cleanUpCoupledMatrix() {
    m_coupledMatrix.reset();
}

void HartreeFockSolver::setupOverlapMatrix() {
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

void HartreeFockSolver::setupUncoupledMatrix() {
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

void HartreeFockSolver::setupCoupledMatrix() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    for(uint p = 0; p < n; p++) {
        for(uint r = 0; r < n; r++) {
            for(uint q = p; q < n; q++) {
                for(uint s = r; s < n; s++) {
                    // NOTE: Indexes changed on purpose in element and call due to
                    // notation differences between Thijssen and Helgaker
                    m_coupledMatrix(p,r)(q,s) = f->coupledIntegral(p, q, r, s);
                }
            }
        }
    }
    for(uint p = 0; p < n; p++) {
        for(uint r = 0; r < n; r++) {
            for(uint q = p; q < n; q++) {
                for(uint s = r; s < n; s++) {
                    double originalValue = m_coupledMatrix(p,r)(q,s);
                    m_coupledMatrix(q,s)(p,r) = originalValue;
                    m_coupledMatrix(q,r)(p,s) = originalValue;
                    m_coupledMatrix(p,s)(q,r) = originalValue;
                    m_coupledMatrix(r,p)(s,q) = originalValue;
                    m_coupledMatrix(s,p)(r,q) = originalValue;
                    m_coupledMatrix(r,q)(s,p) = originalValue;
                    m_coupledMatrix(s,q)(r,p) = originalValue;
                }
            }
        }
    }
}

void HartreeFockSolver::setupFockMatrix() {
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

double HartreeFockSolver::coupledMatrixTilde(int p, int q, int r, int s) { // TODO: Find a better name for tilde
    field<mat>& Q = m_coupledMatrix;
    return 2 * Q(p,r)(q,s) - Q(p,r)(s,q);
}


void HartreeFockSolver::normalizeCoefficientMatrix(uint nParticles, mat &coefficientMatrix){
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
    uint nk = nParticles;

    const mat& S = m_overlapMatrix;
    mat& C = coefficientMatrix;

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

void HartreeFockSolver::advance() {
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;
    setupFockMatrix();

    vec s;
    mat U;
    eig_sym(s, U, m_overlapMatrix);

    mat V = U*diagmat(1.0/sqrt(s));

    mat &F = m_fockMatrix;
    F = V.t() * F * V;

    mat Cprime;
    eig_sym(m_fockEnergies, Cprime, m_fockMatrix);


    mat &C = m_coefficientMatrix;
    C = V*Cprime.submat(0, 0, no - 1, nk - 1);
    normalizeCoefficientMatrix(nk, C);

    setupDensityMatrix();

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

int HartreeFockSolver::iterationsUsed() const
{
    return m_iterationsUsed;
}

void HartreeFockSolver::setupDensityMatrix() {
    mat &P = m_densityMatrix;
    mat &C = m_coefficientMatrix;
    mat tempP = 2 * C * C.t();
//    P = tempP;
    double mixFactor = 0.9;
    if(P.n_elem > 0) {
        P = mixFactor * P + (1 - mixFactor) * tempP; // smoothing
    } else {
        P = tempP;
    }
}

void HartreeFockSolver::solve() {
    for(int i = 0; i < m_nIterationsMax; i++) {
        vec previousFockEnergies = m_fockEnergies;
        advance();
        if(i > 0) {
            if(sum(abs(m_fockEnergies - previousFockEnergies)) / m_fockEnergies.n_elem < m_convergenceTreshold) {
                m_iterationsUsed = i;
                break;
            }
        }
    }
}
double HartreeFockSolver::convergenceTreshold() const
{
    return m_convergenceTreshold;
}

void HartreeFockSolver::setConvergenceTreshold(double convergenceTreshold)
{
    m_convergenceTreshold = convergenceTreshold;
}

void HartreeFockSolver::setElectronSystem(ElectronSystem *basisFunction) {
    m_electronSystem = basisFunction;
}

ElectronSystem *HartreeFockSolver::electronSystem() {
    return m_electronSystem;
}

double HartreeFockSolver::energy()
{
    return m_energy;
}

const mat &HartreeFockSolver::coefficientMatrix() const
{
    return m_coefficientMatrix;
}

const mat &HartreeFockSolver::overlapMatrix() const {
    return m_overlapMatrix;
}

const mat &HartreeFockSolver::densityMatrix() const
{
    return m_densityMatrix;
}

const field<mat> &HartreeFockSolver::coupledMatrix() const
{
    return m_coupledMatrix;
}

int HartreeFockSolver::nIterationsMax() const
{
    return m_nIterationsMax;
}

void HartreeFockSolver::setNIterationsMax(int nIterationsMax)
{
    m_nIterationsMax = nIterationsMax;
}

const mat &HartreeFockSolver::uncoupledMatrix() const
{
    return m_uncoupledMatrix;
}
