#include "unrestrictedhartreefocksolver.h"

#include <electronsystems/electronsystem.h>

UnrestrictedHartreeFockSolver::UnrestrictedHartreeFockSolver(ElectronSystem *system) :
    HartreeFockSolver(system)
{
}

void UnrestrictedHartreeFockSolver::reset() {
    HartreeFockSolver::reset();
    resetCoefficientMatrices();
}

void UnrestrictedHartreeFockSolver::resetCoefficientMatrices() {
    ElectronSystem* f = electronSystem();
    m_coefficientMatrixUp.reset();
    m_coefficientMatrixUp = zeros(f->nBasisFunctions(), f->nParticlesUp());
    m_coefficientMatrixDown.reset();
    m_coefficientMatrixDown = zeros(f->nBasisFunctions(), f->nParticlesDown());
}

void UnrestrictedHartreeFockSolver::setupFockMatrices() {
    ElectronSystem* f = electronSystem();
    uint n = f->nBasisFunctions();

    mat &F = m_fockMatrixUp;
    for (int i = 0; i < 2; ++i) {
        if(i == 1) {
            F = m_fockMatrixDown;
        }
        const mat &P = densityMatrix();
        const mat &h = uncoupledMatrix();
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
}

void UnrestrictedHartreeFockSolver::normalizeCoefficientMatrices(){
    ElectronSystem* f = electronSystem();
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;

    const mat& S = overlapMatrix();
    mat& C = m_coefficientMatrixUp;
    for(int i = 0; i < 2; i++) {
        if(i == 1) {
            C = m_coefficientMatrixDown;
        }
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
}

void UnrestrictedHartreeFockSolver::setupDensityMatrices() {
    mat &P = m_densityMatrixUp;
    mat &C = m_coefficientMatrixUp;
    for (int i = 0; i < 2; ++i) {
        P = m_densityMatrixDown;
        C = m_coefficientMatrixDown;
        mat tempP = 2 * C * C.t();
    //    P = tempP;
        double mixFactor = 0.9;
        if(P.n_elem > 0) {
            P = mixFactor * P + (1 - mixFactor) * tempP; // smoothing
        } else {
            P = tempP;
        }
    }
}

void UnrestrictedHartreeFockSolver::advance() {
    ElectronSystem* f = electronSystem();
    uint no = f->nBasisFunctions();
    uint nkUp = f->nParticlesUp();
    uint nkDown = f->nParticlesDown();
    setupFockMatrices();

    vec s;
    mat U;
    eig_sym(s, U, overlapMatrix());
    mat V = U*diagmat(1.0/sqrt(s));

    mat &F = m_fockMatrixUp;
    mat &C = m_coefficientMatrixUp;
    int nk = nkUp;
    vec &fockEnergies = m_fockEnergiesUp;
    for (int i = 0; i < 2; ++i) {
        if(i == 1) {
            F = m_fockMatrixDown;
            C = m_coefficientMatrixDown;
            nk = nkDown;
            fockEnergies = m_fockEnergiesDown;
        }
        F = V.t() * F * V;

        mat Cprime;
        eig_sym(fockEnergies, Cprime, F);

        C = V*Cprime.submat(0, 0, no - 1, nk - 1);
    }
    normalizeCoefficientMatrices();
    setupDensityMatrices();

    double energy = 0;
    mat& P = m_densityMatrixUp;
    for (int i = 0; i < 2; ++i) {
        if(i == 1) {
            P = m_densityMatrixDown;
        }

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
    }
    energy += electronSystem()->additionalEnergyTerms();
    m_energyUHF = energy;
}

double UnrestrictedHartreeFockSolver::energy()
{
    return m_energyUHF;
}
