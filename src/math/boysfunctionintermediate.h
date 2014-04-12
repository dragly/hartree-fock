#ifndef BOYSFUNCTIONINTERMEDIATE_H
#define BOYSFUNCTIONINTERMEDIATE_H

#include <armadillo>

using namespace arma;

class BoysFunctionIntermediate
{
public:
    BoysFunctionIntermediate(int levelMax = 20, int nValues = 1000, double limitMin = 0.0, double limitMax = 50.0, int nIntegralValues=1e7);

    double result(double arg, int n) const;
    void updateResults();

    // TODO: Ensure we are using the correct levelMax
    static BoysFunctionIntermediate& getInstance() {
        static BoysFunctionIntermediate instance;
        return instance;
    }

protected:
    const int m_taylorExpansionOrder = 6;
    double m_dx;

    uint m_nValues;
    int m_levelMax;
    int m_neededLevelMax;
    double m_limitMin;
    double m_limitMax;
    uint m_nIntegralValues;

    mat m_results;
    vec m_args;
    double integrand(double x, double t, double n) const;
    double directIntegral(double x, double n) const;
    double fastPow(double x, int a) const;

private:
    // Added to protect from copying
    BoysFunctionIntermediate(const BoysFunctionIntermediate&);
    void operator=(const BoysFunctionIntermediate&);

    double factorialInverseTable[40];
};

#endif // BOYSFUNCTIONINTERMEDIATE_H
