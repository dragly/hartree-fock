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

    void reset();
    void advance();

    inline void setElectronSystem(ElectronSystem *electronSystem);
    ElectronSystem *electronSystem();

    inline double energy();
    const mat &coefficientMatrix() const;
    const mat &overlapMatrix();
    double convergenceTreshold() const;
    void setConvergenceTreshold(double convergenceTreshold);

    int iterationsUsed() const;

    void solve();
private:
    mat m_uncoupledMatrix;
    mat m_overlapMatrix;
    mat m_fockMatrix;
    mat m_densityMatrix;
    mat m_coefficientMatrix;
    field<mat> m_coupledMatrix;
    vec m_fockEnergies;

    ElectronSystem *m_electronSystem;

    void resetC();
    void setupF();
    void setupP();
    void setupQ();
    void setupS();
    void setuph();
    void normalizeCwithRegardsToS();
    void cleanUpQMemory();
    void allocateQMemory();

    int m_iterationsUsed;

    double coupledMatrixTilde(int p, int q, int r, int s);
    double *m_QData;
    double m_energy = 0;
    double m_convergenceTreshold;
    double m_previousFockEnergyRMS;
};

inline void HartreeFockSolver::setElectronSystem(ElectronSystem *basisFunction) {
    m_electronSystem = basisFunction;
}

inline ElectronSystem *HartreeFockSolver::electronSystem() {
    return m_electronSystem;
}

inline double HartreeFockSolver::energy()
{
    return m_energy;
}

inline const mat &HartreeFockSolver::coefficientMatrix() const
{
    return m_coefficientMatrix;
}

inline const mat &HartreeFockSolver::overlapMatrix() {
    return m_overlapMatrix;
}

#endif // HARTREEFOCKSOLVER_H
