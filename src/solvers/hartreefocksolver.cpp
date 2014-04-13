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
    m_nIterationsMax(1e3),
    m_densityMixFactor(0.5),
    m_hasBeenSetup(false)
{
}

HartreeFockSolver::~HartreeFockSolver()
{
    cleanUpCoupledMatrix();
}

void HartreeFockSolver::setup()
{
    allocateCoupledMatrix();
    setupIntegralMatrices();
    m_hasBeenSetup = true;
}

void HartreeFockSolver::advance()
{
    if(!m_hasBeenSetup) {
        setup();
    }
}

void HartreeFockSolver::setupIntegralMatrices()
{
    setupUncoupledMatrix();
    setupOverlapMatrix();
    setupCoupledMatrix();
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

int HartreeFockSolver::iterationsUsed() const
{
    return m_iterationsUsed;
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

const mat &HartreeFockSolver::overlapMatrix() const {
    return m_overlapMatrix;
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

double HartreeFockSolver::densityMixFactor() const
{
    return m_densityMixFactor;
}

void HartreeFockSolver::setDensityMixFactor(double densityMixFactor)
{
    m_densityMixFactor = densityMixFactor;
}
