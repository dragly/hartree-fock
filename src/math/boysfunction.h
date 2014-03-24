#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <armadillo>

class BoysFunctionIntermediate;

using namespace arma;

class BoysFunction
{
public:
    BoysFunction(double arg, int levelMax = 0, BoysFunctionIntermediate* intermediate = 0);

    double result(int level = 0) const;

    // Static helper functions
    static double calculateAsymptopticForm(double arg, int level);
    static double calculateTaylorExpansion(double arg, int level);
    static double calculateZeroLevel(double arg);
    static double doubleFactorial(double n);
    static double factorial(double n);
protected:
    vec m_results;
    BoysFunctionIntermediate* m_intermediate;
};

#endif // BOYSFUNCTION_H
