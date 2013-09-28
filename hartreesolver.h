#ifndef HARTREESOLVER_H
#define HARTREESOLVER_H

#include <armadillo>

using namespace arma;

// Forward declarations
class BasisFunction;

class HartreeSolver
{
public:
    explicit HartreeSolver(BasisFunction *basisFunction);

    virtual ~HartreeSolver();

    double alpha[4];

    void reset();
    void advance();

    inline void setBasisFunction(BasisFunction *basisFunction);
private:
    mat h;
    mat S;
    mat F;
    vec C;

    double ****Q;

    BasisFunction *m_basisFunction;

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
};

inline void HartreeSolver::setBasisFunction(BasisFunction *basisFunction) {
    m_basisFunction = basisFunction;
}

#endif // HARTREESOLVER_H
