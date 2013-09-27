#ifndef HARTREESOLVER_H
#define HARTREESOLVER_H

#include <armadillo>

using namespace arma;

class HartreeSolver
{
public:
    explicit HartreeSolver();

    virtual ~HartreeSolver();

    double alpha[4];
    double basisFunction(int index, vec position);
    double matrixElement(int p, int r, int q, int s);

    const double powPi5over2 = pow(M_PI, 5./2.);
    double kineticIntegral(int p, int q);
    double nuclearAttractionIntegral(int p, int q);
    double overlapIntegral(int p, int q);

    void reset();
    void advance();
private:
    mat h;
    mat S;
    mat F;
    vec C;

    double ****Q;

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

#endif // HARTREESOLVER_H
