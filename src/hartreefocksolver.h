#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <armadillo>

using namespace arma;

// Forward declarations
class ElectronSystem;

class HartreeFockSolver
{
public:
    explicit HartreeFockSolver(ElectronSystem *basisFunction);

    virtual ~HartreeFockSolver();

    double alpha[4];

    void reset();
    void advance();

    inline void setBasisFunction(ElectronSystem *basisFunction);
    ElectronSystem *basisFunction();

    inline double energy();
private:
    mat h;
    mat S;
    mat F;
    mat P;
    mat C;

    double ****Q;

    ElectronSystem *m_basisFunction;

    void resetC();
    void setupF();
    void setupP();
    void setupQ();
    void setupS();
    void setuph();
    void normalizeCwithRegardsToS();
    double Qtilde(int p, int q, int r, int s);

    void cleanUpQMemory();
    void allocateQMemory();

    double *QData;

    bool isQAllocated = false;

    double m_energy = 0;
};

inline void HartreeFockSolver::setBasisFunction(ElectronSystem *basisFunction) {
    m_basisFunction = basisFunction;
}

inline ElectronSystem *HartreeFockSolver::basisFunction() {
    return m_basisFunction;
}

inline double HartreeFockSolver::energy()
{
    return m_energy;
}

#endif // HARTREEFOCKSOLVER_H
