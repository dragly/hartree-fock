#ifndef HARTREESOLVER_H
#define HARTREESOLVER_H

#include <armadillo>

using namespace arma;

// Forward declarations
class ElectronSystem;

class HartreeSolver
{
public:
    explicit HartreeSolver(ElectronSystem *basisFunction);

    virtual ~HartreeSolver();

    double alpha[4];

    void reset();
    void advance();

    inline void setBasisFunction(ElectronSystem *basisFunction);

    inline double energy();
private:
    mat h;
    mat S;
    mat F;
    vec C;

    double ****Q;

    ElectronSystem *m_basisFunction;

    void resetC();
    void setupF();
    void setupQ();
    void setupS();
    void setuph();
    void setupAlpha();
    void normalizeCwithRegardsToS();

    void cleanUpQMemory();
    void allocateQMemory();

    double *QData;

    bool isQAllocated = false;

    double m_energy = 0;
};

inline void HartreeSolver::setBasisFunction(ElectronSystem *basisFunction) {
    m_basisFunction = basisFunction;
}

inline double HartreeSolver::energy()
{
    return m_energy;
}

#endif // HARTREESOLVER_H
