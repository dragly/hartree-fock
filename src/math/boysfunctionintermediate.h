#ifndef BOYSFUNCTIONINTERMEDIATE_H
#define BOYSFUNCTIONINTERMEDIATE_H

#include <armadillo>

using namespace arma;

class BoysFunctionIntermediate
{
public:
    BoysFunctionIntermediate(int levelMax, int nValues = 10000, double limitMin = 0.09, double limitMax = 30.0, int nIntegralValues=1e4);

    double result(double arg, int n) const;
    void updateResults();
protected:
    const int m_taylorExpansionOrder = 6;
    double m_dx;

    uint m_nValues;
    int m_levelMax;
    double m_limitMin;
    double m_limitMax;
    uint m_nIntegralValues;

    mat m_results;
    vec m_args;
    double integrand(double x, double t, double n) const;
    double directIntegral(double x, double n) const;
    double fastPow(double x, int a) const;
};

#endif // BOYSFUNCTIONINTERMEDIATE_H
