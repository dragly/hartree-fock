#include "unrestrictedhartreefocksolver.h"
#include "electronsystems/electronsystem.h"

#include <iomanip>

using std::setprecision;
using std::fixed;

UnrestrictedHartreeFockSolver::UnrestrictedHartreeFockSolver(ElectronSystem *system) :
    HartreeFockSolver(system),
    m_initialCoefficientMatricesSetManually(false)
{
}

void UnrestrictedHartreeFockSolver::setup()
{
    HartreeFockSolver::setup();
    resetCoefficientMatrices();
    resetFockMatrices();
    setupDensityMatrices();
    if(m_densityMatrixUp.n_rows > 0 && m_densityMatrixUp.n_cols > 1) {
        m_densityMatrixUp(0,1) = 0.1; // Added asymmetry between the spin up and spin down orbitals
    }
    setupFockMatrices();
}

void UnrestrictedHartreeFockSolver::resetFockMatrices() {
    uint n = electronSystem()->nBasisFunctions();
    m_fockMatrixUp = zeros(n,n);
    m_fockMatrixDown = zeros(n,n);
}

void UnrestrictedHartreeFockSolver::resetCoefficientMatrices() {
    ElectronSystem* f = electronSystem();

    if(!m_initialCoefficientMatricesSetManually
            || m_initialCoefficientMatrixUp.n_rows != f->nBasisFunctions()
            || m_initialCoefficientMatrixUp.n_cols != f->nParticlesUp()
            || m_initialCoefficientMatrixDown.n_rows != f->nBasisFunctions()
            || m_initialCoefficientMatrixDown.n_cols != f->nParticlesDown()) {
        m_initialCoefficientMatrixUp = zeros(f->nBasisFunctions(), f->nParticlesUp());
        m_initialCoefficientMatrixDown = zeros(f->nBasisFunctions(), f->nParticlesDown());
    }

    m_coefficientMatrixUp = m_initialCoefficientMatrixUp;
    m_coefficientMatrixDown = m_initialCoefficientMatrixDown;
}

void UnrestrictedHartreeFockSolver::randomizeCoefficientMatrices() {
    ElectronSystem* f = electronSystem();
    m_coefficientMatrixUp = randn(f->nBasisFunctions(), f->nParticlesUp());
    m_coefficientMatrixDown = randn(f->nBasisFunctions(), f->nParticlesDown());
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
    for(uint p = 0; p < n; p++) {
        for(uint q = p; q < n; q++) {
            double hpq = h(p,q);
            double Fupq = hpq;
            double Fdpq = hpq;
            // TODO: Apply the symmetry in Q for increased performance
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    double Qprqs = Q(p,r)(q,s);
                    double Qprsq = Q(p,r)(s,q);
                    double QprqsMinusQprsq = Qprqs - Qprsq;
                    double Pusr = Pu(s,r);
                    double Pdsr = Pd(s,r);
                    Fupq += Pusr * QprqsMinusQprsq + Pdsr * Qprqs;
                    Fdpq += Pdsr * QprqsMinusQprsq + Pusr * Qprqs;
                }
            }
            // Set elements and apply symmetry
            Fu(p,q) = Fupq;
            Fu(q,p) = Fupq;
            Fd(p,q) = Fdpq;
            Fd(q,p) = Fdpq;
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
    HartreeFockSolver::advance();
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

    if(nkUp > 0) {
        Cu = V*CprimeUp.submat(0, 0, no - 1, nkUp - 1);
    }
    if(nkDn > 0) {
        Cd = V*CprimeDn.submat(0, 0, no - 1, nkDn - 1);
    }

    normalizeCoefficientMatrix(f->nParticlesUp(), m_coefficientMatrixUp);
    normalizeCoefficientMatrix(f->nParticlesDown(), m_coefficientMatrixDown);

    setupDensityMatrices();
    setupFockMatrices();
}

void UnrestrictedHartreeFockSolver::solve() {
    HartreeFockSolver::solve();
    m_noConvergenceCounter = 0;
    for(int i = 0; i < nIterationsMax(); i++) {
        m_iterationsUsed = i + 1;
        vec previousFockEnergies = m_fockEnergiesUp;
        advance();
        if(i > 0) {
            double stdDeviation = sum(abs(m_fockEnergiesUp - previousFockEnergies)) / m_fockEnergiesUp.n_elem;
            if(stdDeviation < convergenceTreshold()) {
                break;
            }
        }
    }
    calculateEnergy();
}

void UnrestrictedHartreeFockSolver::calculateEnergy()
{
    double energy = 0;
    mat &Pu = m_densityMatrixUp;
    mat &Pd = m_densityMatrixDown;
    mat &Fu = m_fockMatrixUp;
    mat &Fd = m_fockMatrixDown;

    const mat& h = uncoupledMatrix();
    energy += 0.5 * accu( (Pu + Pd) % h + Fu % Pu + Fd % Pd);
    energy += electronSystem()->additionalEnergyTerms();
    m_energyUHF = energy;
}


const mat &UnrestrictedHartreeFockSolver::coeffcientMatrixUp()
{
    return m_coefficientMatrixUp;
}

const mat &UnrestrictedHartreeFockSolver::coeffcientMatrixDown()
{
    return m_coefficientMatrixDown;
}

void UnrestrictedHartreeFockSolver::setInitialCoefficientMatrices(const mat &up, const mat &down)
{
    m_initialCoefficientMatricesSetManually = true;
    m_initialCoefficientMatrixUp = up;
    m_initialCoefficientMatrixDown = down;
}

double UnrestrictedHartreeFockSolver::energy()
{
    return m_energyUHF;
}
