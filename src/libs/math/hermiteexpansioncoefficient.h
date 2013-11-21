#ifndef HERMITEEXPANSIONCOEFFICIENT_H
#define HERMITEEXPANSIONCOEFFICIENT_H

#include <armadillo>

using namespace arma;
using namespace std;

class HermiteExpansionCoefficient
{
public:
    HermiteExpansionCoefficient(int dimension);
    explicit HermiteExpansionCoefficient(double a, double b, rowvec A, rowvec B, int angularMomentumMax, bool setupImmediately = true);

    void reset(int dimension);
    void set(double a, double b, rowvec A, rowvec B, bool setupImmediately = true);

    const cube &operator [](const uword row) const;
    cube &operator[](const uword row);
    void setupE();
    bool checkIndexCombinationForE(int iA, int iB, int t);
    double operator ()(int iA, int jA, int kA, int iB, int jB, int kB, const uword t, const uword u, const uword v) const;
protected:
    double m_a;
    double m_b;
    rowvec m_A;
    rowvec m_B;
    int m_dimension;
    cube m_E[3]; // t, i, j
    rowvec m_AB;
    rowvec m_P;
    rowvec m_PA;
    rowvec m_PB;
};

inline const cube &HermiteExpansionCoefficient::operator[](const uword dimension) const {
    return m_E[dimension];
}

inline cube &HermiteExpansionCoefficient::operator[](const uword row)
{
    return const_cast<cube&>(static_cast<const HermiteExpansionCoefficient&>(*this)[row]);
}

inline double HermiteExpansionCoefficient::operator()(int iA, int jA, int kA, int iB, int jB, int kB, const uword t, const uword u, const uword v) const {
    return m_E[0](iA, iB, t) * m_E[1](jA, jB, u) * m_E[2](kA, kB, v);
}

#endif // HERMITEEXPANSIONCOEFFICIENT_H
