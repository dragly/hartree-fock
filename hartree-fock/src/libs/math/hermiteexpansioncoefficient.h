#ifndef HERMITEEXPANSIONCOEFFICIENT_H
#define HERMITEEXPANSIONCOEFFICIENT_H

#include <armadillo>

using namespace arma;
using namespace std;

class HermiteExpansionCoefficient
{
public:
    HermiteExpansionCoefficient();
    explicit HermiteExpansionCoefficient(double a, double b, rowvec A, rowvec B, int angularMomentumMax, bool setupImmediately = true);

    const cube &operator [](const uword row) const;
    cube &operator[](const uword row);
    void setupE();
    bool checkIndexCombinationForE(int iA, int iB, int t);
protected:
    double m_a;
    double m_b;
    rowvec m_A;
    rowvec m_B;
    int m_angularMomentumMax;
    cube m_E[3]; // t, i, j
};

inline const cube &HermiteExpansionCoefficient::operator[](const uword dimension) const {
    return m_E[dimension];
}

inline cube &HermiteExpansionCoefficient::operator[](const uword row)
{
    return const_cast<cube&>(static_cast<const HermiteExpansionCoefficient&>(*this)[row]);
}

#endif // HERMITEEXPANSIONCOEFFICIENT_H
