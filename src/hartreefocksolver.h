#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <armadillo>

using namespace arma;

// Forward declarations
class ElectronSystem;

class HartreeFockSolver
{
public:
    explicit HartreeFockSolver(ElectronSystem *electronSystem);

    virtual ~HartreeFockSolver();

    double alpha[4];

    virtual void reset();
    virtual void advance();

    void setElectronSystem(ElectronSystem *electronSystem);
    ElectronSystem *electronSystem();

    virtual double energy();
    const mat& coefficientMatrix() const;
    double convergenceTreshold() const;
    void setConvergenceTreshold(double convergenceTreshold);

    int iterationsUsed() const;

    void solve();
    int nIterationsMax() const;
    void setNIterationsMax(int nIterationsMax);

    const mat &uncoupledMatrix() const;
    const mat &overlapMatrix() const;
    const mat &densityMatrix() const;
    const field<mat> &coupledMatrix() const;

    double densityMixFactor() const;
    void setDensityMixFactor(double densityMixFactor);

protected:
    void normalizeCoefficientMatrix(uint nParticles, mat &coefficientMatrix);
    double coupledMatrixTilde(int p, int q, int r, int s);
private:
    mat m_uncoupledMatrix;
    mat m_overlapMatrix;
    mat m_fockMatrix;
    mat m_densityMatrix;
    mat m_coefficientMatrix;
    field<mat> m_coupledMatrix;
    vec m_fockEnergies;

    ElectronSystem *m_electronSystem;

    void resetCoefficientMatrix();
    void setupFockMatrix();
    void setupDensityMatrix();
    void setupCoupledMatrix();
    void setupOverlapMatrix();
    void setupUncoupledMatrix();
    void cleanUpCoupledMatrix();
    void allocateCoupledMatrix();

    int m_iterationsUsed;

    double *m_QData;
    double m_energy = 0;
    double m_convergenceTreshold;
    double m_previousFockEnergyRMS;
    int m_nIterationsMax;
    double m_densityMixFactor;
};

#endif // HARTREEFOCKSOLVER_H
