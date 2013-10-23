#ifndef BOYSFUNCTIONINTERMEDIATE_H
#define BOYSFUNCTIONINTERMEDIATE_H

#include <armadillo>

using namespace arma;

class BoysFunctionIntermediate
{
public:
    BoysFunctionIntermediate(int levelMax, int nValues = 10000, double limitMin = 0.09, double limitMax = 30.0);

    double result(double arg, int n) const;
    void updateResults();
protected:
    double m_limitMin;
    double m_limitMax;
    const int m_taylorExpansionOrder = 6;
    double m_dx;

    uint m_nValues;
    int m_levelMax;

    mat m_results;
    vec m_args;
};

#endif // BOYSFUNCTIONINTERMEDIATE_H
