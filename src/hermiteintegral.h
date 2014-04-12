#ifndef HERMITEINTEGRAL_H
#define HERMITEINTEGRAL_H

#include <armadillo>

using namespace arma;
using namespace std;

class HermiteIntegral
{
public:
    HermiteIntegral(int dimensionMax);
    HermiteIntegral(double alpha, const rowvec &A, int dimension, bool setupImmediately = true);

    void reset(int dimension);
    void set(double alpha, const rowvec &A, bool setupImmediately = true);
    void set(double alpha, const rowvec &A, int dimension);
    void setupR();

    const cube &operator [](const uword row) const;
    cube &operator[](const uword row);
    const cube &operator()(const uword row) const;
    cube &operator()(const uword row);
    double operator ()(const uword n, const uword t, const uword u, const uword v) const;
protected:
    double m_alpha;
    rowvec m_A;
    field<cube> m_R;
    int m_dimensionMax;
};

inline const cube &HermiteIntegral::operator()(const uword row) const {
    return m_R(row);
}

inline double HermiteIntegral::operator()(const uword n, const uword t, const uword u, const uword v) const {
    return m_R(n)(t,u,v);
}

inline const cube &HermiteIntegral::operator[](const uword row) const {
    return m_R[row];
}

// Non-const version based on Scott Meyer's Item 3 in Effective C++ Third Edition

inline cube &HermiteIntegral::operator()(const uword row)
{
    return const_cast<cube&>(static_cast<const HermiteIntegral&>(*this)(row));
}

inline cube &HermiteIntegral::operator[](const uword row)
{
    return const_cast<cube&>(static_cast<const HermiteIntegral&>(*this)[row]);
}

#endif // HERMITEINTEGRAL_H
