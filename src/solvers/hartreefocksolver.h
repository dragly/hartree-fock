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

    virtual void setup();
    virtual void advance();
    virtual void solve();

    void setElectronSystem(ElectronSystem *electronSystem);
    ElectronSystem *electronSystem();

    virtual double energy() = 0;
    double convergenceTreshold() const;
    void setConvergenceTreshold(double convergenceTreshold);

    int iterationsUsed() const;

    int nIterationsMax() const;
    void setNIterationsMax(int nIterationsMax);

    const mat &uncoupledMatrix() const;
    const mat &overlapMatrix() const;
    const mat &transformationMatrix() const;
    const field<mat> &coupledMatrix() const;

    double densityMixFactor() const;
    void setDensityMixFactor(double densityMixFactor);

protected:
    void setupIntegralMatrices();
    void normalizeCoefficientMatrix(uint nParticles, mat &coefficientMatrix);
    int m_iterationsUsed;
private:
    void setupCoupledMatrix();
    void setupOverlapMatrix();
    void setupUncoupledMatrix();
    void cleanUpCoupledMatrix();
    void allocateCoupledMatrix();
    void setupTransformationMatrix();

    ElectronSystem *m_electronSystem;

    // Matrices
    mat m_uncoupledMatrix;
    mat m_overlapMatrix;
    mat m_transformationMatrix;
    field<mat> m_coupledMatrix;

    // Variables
    double *m_QData;
    double m_convergenceTreshold;
    double m_previousFockEnergyRMS;
    int m_nIterationsMax;
    double m_densityMixFactor;
    bool m_hasBeenSetup;
};

#endif // HARTREEFOCKSOLVER_H
