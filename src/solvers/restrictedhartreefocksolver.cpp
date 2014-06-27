#include "restrictedhartreefocksolver.h"

#include <electronsystems/electronsystem.h>

/*!
 * \class RestrictedHartreeFockSolver
 * \brief Solver for the Roothan equation
 */

RestrictedHartreeFockSolver::RestrictedHartreeFockSolver(ElectronSystem *electronSystem) :
    HartreeFockSolver(electronSystem),
    m_initialCoefficientMatrixSetManually(false)
{
}

RestrictedHartreeFockSolver::~RestrictedHartreeFockSolver()
{

}

void RestrictedHartreeFockSolver::resetFockMatrix() {
    uint n = electronSystem()->nBasisFunctions();
    m_fockMatrix = zeros(n,n);
}

void RestrictedHartreeFockSolver::resetCoefficientMatrix() {
    ElectronSystem* f = electronSystem();
    m_coefficientMatrix.reset();
    if(!m_initialCoefficientMatrixSetManually
            || m_initialCoefficientMatrix.n_rows != f->nBasisFunctions()
            || m_initialCoefficientMatrix.n_cols != f->nParticlesUp()) {
        m_initialCoefficientMatrix = zeros(f->nBasisFunctions(), f->nParticles() / 2);
    }
    m_coefficientMatrix = m_initialCoefficientMatrix;
}

void RestrictedHartreeFockSolver::setupFockMatrix() {
    ElectronSystem* f = electronSystem();
    uint n = f->nBasisFunctions();
    mat &F = m_fockMatrix;
    const mat &P = m_densityMatrix;
    const mat &h = uncoupledMatrix();
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

void RestrictedHartreeFockSolver::setInitialCoefficientMatrix(const mat &coefficients)
{
    m_initialCoefficientMatrixSetManually = true;
    m_initialCoefficientMatrix = coefficients;
}

void RestrictedHartreeFockSolver::setup() {
    HartreeFockSolver::setup();
    resetCoefficientMatrix();
    resetFockMatrix();
    setupDensityMatrix();
}

void RestrictedHartreeFockSolver::advance() {
    HartreeFockSolver::advance();
    ElectronSystem* f = electronSystem();
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;
    setupFockMatrix();

    const mat &V = transformationMatrix();

    const mat &F = m_fockMatrix;
    mat Fprime = V.t() * F * V;

    mat Cprime;
    eig_sym(m_fockEnergies, Cprime, Fprime);


    mat &C = m_coefficientMatrix;
    C = V*Cprime.submat(0, 0, no - 1, nk - 1);
    normalizeCoefficientMatrix(nk, C);

    setupDensityMatrix();

    double energy = 0;

    const mat& P = m_densityMatrix;
    const mat& h = uncoupledMatrix();
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
    energy += electronSystem()->additionalEnergyTerms();
    m_energy = energy;
}

void RestrictedHartreeFockSolver::solve() {
    HartreeFockSolver::solve();
    for(int i = 0; i < nIterationsMax(); i++) {
        vec previousFockEnergies = m_fockEnergies;
        advance();
        if(i > 0) {
            double averageEnergyChange = sum(abs(m_fockEnergies - previousFockEnergies))
                    / m_fockEnergies.n_elem;
            if(averageEnergyChange < convergenceTreshold()) {
                m_iterationsUsed = i;
                break;
            }
        }
    }
}

double RestrictedHartreeFockSolver::energy()
{
    return m_energy;
}

const mat &RestrictedHartreeFockSolver::densityMatrix() const
{
    return m_densityMatrix;
}

const mat &RestrictedHartreeFockSolver::coefficientMatrix() const
{
    return m_coefficientMatrix;
}

double RestrictedHartreeFockSolver::coupledMatrixTilde(int p, int q, int r, int s)
{
    const field<mat>& Q = coupledMatrix();
    return 2 * Q(p,r)(q,s) - Q(p,r)(s,q);
}

void RestrictedHartreeFockSolver::setupDensityMatrix() {
    mat &P = m_densityMatrix;
    mat &C = m_coefficientMatrix;
    mat tempP = 2 * C * C.t();
//    P = tempP;

    if(P.n_elem > 0) {
        P = densityMixFactor() * P + (1 - densityMixFactor()) * tempP; // smoothing
    } else {
        P = tempP;
    }
}

