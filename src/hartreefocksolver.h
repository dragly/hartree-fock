#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <armadillo>

using namespace arma;

// Forward declarations
class BasisFunction;

class HartreeFockSolver
{
public:
    explicit HartreeFockSolver(BasisFunction *basisFunction);

    virtual ~HartreeFockSolver();

    double alpha[4];

    void reset();
    void advance();

    inline void setBasisFunction(BasisFunction *basisFunction);

    inline double energy();
private:
    mat h;
    mat S;
    mat F;
    mat P;
    mat C;

    double ****Q;

    BasisFunction *m_basisFunction;

    void resetC();
    void setupF();
    void setupP();
    void setupQ();
    void setupS();
    void setuph();
    void setupAlpha();
    void normalizeCwithRegardsToS();
    double Qtilde(int p, int q, int r, int s);

    void cleanUpQMemory();
    void allocateQMemory();

    double *QData;

    bool isQAllocated = false;

    double m_energy = 0;
};

inline void HartreeFockSolver::setBasisFunction(BasisFunction *basisFunction) {
    m_basisFunction = basisFunction;
}

inline double HartreeFockSolver::energy()
{
    return m_energy;
}

#endif // HARTREEFOCKSOLVER_H
