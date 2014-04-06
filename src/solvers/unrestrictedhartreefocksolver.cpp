#include "unrestrictedhartreefocksolver.h"

#include <electronsystems/electronsystem.h>
#include <iomanip>

using std::setprecision;
using std::fixed;

UnrestrictedHartreeFockSolver::UnrestrictedHartreeFockSolver(ElectronSystem *system) :
    HartreeFockSolver(system)
{
    resetCoefficientMatrices();
    setupDensityMatrices();
    setupFockMatrices();
}

void UnrestrictedHartreeFockSolver::reset() {
    HartreeFockSolver::reset();
    resetCoefficientMatrices();
    setupDensityMatrices();
    setupFockMatrices();
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

    mat &Fu = m_fockMatrixUp;
    mat &Fd = m_fockMatrixDown;
    const mat &Pu = m_densityMatrixUp;
    const mat &Pd = m_densityMatrixDown;
    const mat &h = uncoupledMatrix();
    const field<mat> &Q = coupledMatrix();
    Fu = zeros(n,n);
    Fd = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            Fu(p,q) = h(p,q);
            Fd(p,q) = h(p,q);
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    Fu(p,q) += Pu(s,r) * (Q(p,r)(q,s) - Q(p,r)(s,q)) + Pd(s,r) * Q(p,r)(q,s);
                    Fd(p,q) += Pd(s,r) * (Q(p,r)(q,s) - Q(p,r)(s,q)) + Pu(s,r) * Q(p,r)(q,s);
                }
            }
        }
    }
}

void UnrestrictedHartreeFockSolver::setupDensityMatrices() {
    mat &Pu = m_densityMatrixUp;
    mat &Pd = m_densityMatrixDown;
    mat &Cu = m_coefficientMatrixUp;
    mat &Cd = m_coefficientMatrixDown;
    mat tempPu = Cu * Cu.t();
    mat tempPd = Cd * Cd.t();
    double mixFactor = 0.5;
    if(Pu.n_elem > 0 && Pd.n_elem > 0) {
        Pu = mixFactor * Pu + (1 - mixFactor) * tempPu; // smoothing
        Pd = mixFactor * Pd + (1 - mixFactor) * tempPd; // smoothing
    } else {
        Pu = tempPu;
        Pd = tempPd;
    }
}

void UnrestrictedHartreeFockSolver::advance() {
    ElectronSystem* f = electronSystem();
    uint no = f->nBasisFunctions();
    uint nkUp = f->nParticlesUp();
    uint nkDn = f->nParticlesDown();

    vec s;
    mat U;
    eig_sym(s, U, overlapMatrix());
    mat V = U*diagmat(1.0/sqrt(s));

    mat &Fu = m_fockMatrixUp;
    mat &Fd = m_fockMatrixDown;
    mat &Cu = m_coefficientMatrixUp;
    mat &Cd = m_coefficientMatrixDown;
    vec &fockEnergiesUp = m_fockEnergiesUp;
    vec &fockEnergiesDn = m_fockEnergiesDown;

    mat FprimeUp = V.t() * Fu * V;
    mat FprimeDn = V.t() * Fd * V;

    mat CprimeUp;
    mat CprimeDn;

    eig_sym(fockEnergiesUp, CprimeUp, FprimeUp);
    eig_sym(fockEnergiesDn, CprimeDn, FprimeDn);

    Cu = V*CprimeUp.submat(0, 0, no - 1, nkUp - 1);
    Cd = V*CprimeDn.submat(0, 0, no - 1, nkDn - 1);

    normalizeCoefficientMatrix(f->nParticlesUp(), m_coefficientMatrixUp);
    normalizeCoefficientMatrix(f->nParticlesDown(), m_coefficientMatrixDown);

    setupDensityMatrices();
    setupFockMatrices();

    double energy = 0;
    mat &Pu = m_densityMatrixUp;
    mat &Pd = m_densityMatrixDown;

    const mat& h = uncoupledMatrix();

    energy += 0.5 * accu( (Pu + Pd) % h + Fu % Pu + Fd % Pd);

    energy += electronSystem()->additionalEnergyTerms();
    m_energyUHF = energy;
}

void UnrestrictedHartreeFockSolver::solve() {
    setupFockMatrices();
    for(int i = 0; i < nIterationsMax(); i++) {
        vec previousFockEnergies = m_fockEnergiesUp;
        advance();
        if(i > 0) {
            if(sum(abs(m_fockEnergiesUp - previousFockEnergies)) / m_fockEnergiesUp.n_elem < convergenceTreshold()) {
                m_iterationsUsed = i;
                break;
            }
        }
    }
}

double UnrestrictedHartreeFockSolver::energy()
{
    return m_energyUHF;
}
