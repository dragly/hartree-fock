#include "unrestrictedhartreefocksolver.h"
#include "electronsystems/electronsystem.h"

#include <iomanip>

using std::setprecision;
using std::fixed;

UnrestrictedHartreeFockSolver::UnrestrictedHartreeFockSolver(ElectronSystem *system) :
    HartreeFockSolver(system),
    m_initialCoefficientMatricesSetManually(false),
    m_diisSampleCount(10),
    m_diisStartingIteration(20),
    m_isDiisEnabled(true)
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
    for(int i = 0; i < nIterationsMax(); i++) {
        m_iterationsUsed = i + 1;
        vec previousFockEnergies = m_fockEnergiesUp;
        if(m_isDiisEnabled && i > m_diisStartingIteration) {
            performDIIS();
        }
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


const mat &UnrestrictedHartreeFockSolver::coeffcientMatrixUp() const
{
    return m_coefficientMatrixUp;
}

const mat &UnrestrictedHartreeFockSolver::coeffcientMatrixDown() const
{
    return m_coefficientMatrixDown;
}

const mat &UnrestrictedHartreeFockSolver::densityMatrixUp() const
{
    return m_densityMatrixUp;
}

const mat &UnrestrictedHartreeFockSolver::densityMatrixDown() const
{
    return m_densityMatrixDown;
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

/*!
 * \brief UnrestrictedHartreeFockSolver::performDIIS helps convergence by the DIIS procedure.
 *
 * This implementation is a forked version of one made by Milad Hobbi Mobarhan
 * (source: https://github.com/miladh/HF)
 */
void UnrestrictedHartreeFockSolver::performDIIS()
{
    mat &Fu = m_fockMatrixUp;
    mat &Fd = m_fockMatrixDown;
    mat &Pu = m_densityMatrixUp;
    mat &Pd = m_densityMatrixDown;
    const mat &S = overlapMatrix();

    m_errorsU.push_back(Fu*Pu*S - S*Pu*Fu);
    m_fockMatricesU.push_back(m_fockMatrixUp);

    m_errorsD.push_back(Fd*Pd*S - S*Pd*Fd);
    m_fockMatricesD.push_back(m_fockMatrixDown);

    if(signed(m_errorsU.size()) > m_diisSampleCount){
        m_errorsU.erase(m_errorsU.begin());
        m_fockMatricesU.erase(m_fockMatricesU.begin());

        m_errorsD.erase(m_errorsD.begin());
        m_fockMatricesD.erase(m_fockMatricesD.begin());

        mat Au = zeros(m_errorsU.size()+1, m_errorsU.size()+1);
        mat Ad = zeros(m_errorsD.size()+1, m_errorsD.size()+1);

        for(uint i = 0; i < m_errorsU.size(); i++){
            Au(i, m_errorsU.size()) = -1;
            Au(m_errorsU.size(), i) = -1;

            Ad(i, m_errorsD.size()) = -1;
            Ad(m_errorsD.size(), i) = -1;

            for(uint j = 0; j < m_errorsU.size(); j++){
                Au(i,j) = trace(m_errorsU.at(i) * m_errorsU.at(j));
                Ad(i,j) = trace(m_errorsD.at(i) * m_errorsD.at(j));
            }
        }

        vec bu = zeros(Au.n_rows);
        vec bd = zeros(Ad.n_rows);

        bu(bu.n_elem -1) = -1.0;
        bd(bd.n_elem -1) = -1.0;
        bu = inv(Au) * bu;
        bd = inv(Ad) * bd;

        for(uint i = 0; i < bu.n_elem - 1; i++){
            m_fockMatrixUp += bu(i) * m_fockMatricesU.at(i);
            m_fockMatrixDown += bd(i) * m_fockMatricesD.at(i);
        }
    }
}
int UnrestrictedHartreeFockSolver::diisStartingIteration() const
{
    return m_diisStartingIteration;
}

void UnrestrictedHartreeFockSolver::setDiisStartingIteration(int diisStartingIteration)
{
    m_diisStartingIteration = diisStartingIteration;
}

int UnrestrictedHartreeFockSolver::diisSampleCount() const
{
    return m_diisSampleCount;
}

void UnrestrictedHartreeFockSolver::setDiisSampleCount(int diisIterationCount)
{
    m_diisSampleCount = diisIterationCount;
}

bool UnrestrictedHartreeFockSolver::isDiisEnabled() const
{
    return m_isDiisEnabled;
}

void UnrestrictedHartreeFockSolver::setDiisEnabled(bool useDIIS)
{
    m_isDiisEnabled = useDIIS;
}

