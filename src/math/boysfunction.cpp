#include "boysfunction.h"
#include <iomanip>
#include <math/boysfunctionintermediate.h>

using namespace std;

/*!
 * \class BoysFunction
 * \brief Calculates the Boys function efficiently for any argument and level
 *
 * The argument and the highest level is specified upon construction,
 * while the other levels are fetched on demand.
 */

BoysFunction::BoysFunction()
{

}

BoysFunction::BoysFunction(double arg, int levelMax, BoysFunctionIntermediate *intermediate)
{
    set(arg, levelMax, intermediate);
}

void BoysFunction::set(double arg, int levelMax, BoysFunctionIntermediate *intermediate)
{
    double limitMin = 0;
    double limitMax = 50;
    if(intermediate == 0) {
        m_intermediate = &(BoysFunctionIntermediate::getInstance());
    } else {
        m_intermediate = intermediate;
    }
    double x = arg;
    m_results = zeros(levelMax + 1);
    double expmx = exp(-x);
    if(arg < limitMin || arg > limitMax) {
        if(arg < limitMin) {
            m_results(levelMax) = calculateTaylorExpansion(arg, levelMax);
        } else if(arg > limitMax) {
            m_results(levelMax) = calculateAsymptopticForm(arg, levelMax);
        }
        // Iterate down
        for(int n = levelMax; n > 0; n--) {
            double Fn = m_results(n);
            m_results(n - 1) = (2*x*Fn + expmx) / (2*n - 1);
        }
    } else {
        // Get results from intermediate area
        for(int n = 0; n < levelMax + 1; n++) {
            m_results(n) = m_intermediate->result(arg, n);
        }
    }
}

double BoysFunction::calculateAsymptopticForm(double arg, int level) {
    double n = level;
    double dfac = doubleFactorial(2 * n - 1);
    double result = dfac / pow(2, n + 1) * sqrt(M_PI / pow(arg, 2*n + 1));
    return result;
}

double BoysFunction::doubleFactorial(double n) {
    if(n == -1) {
        return 1;
    }
    double result = 1;
    for(double i = n; i >= 1; i -= 2) {
        result *= i;
    }
    return result;
}

double BoysFunction::factorial(double n) {
    if(n == 0) {
        return 1;
    }
    double result = 1;
    double startPoint = 2;
    for(double i = startPoint; i <= n; i++) {
        result *= i;
    }
    return result;
}

double BoysFunction::calculateTaylorExpansion(double arg, int level) {
    double n = level;
    double result = 0;
    for(int k = 0; k < 100; k++) {
        result += pow(-arg, k) / (factorial(k) * (2 * n + 2 * k + 1));
    }
    return result;
}

double BoysFunction::calculateZeroLevel(double arg) {
    if (arg < 1.0E-6){
        return 1.0;
    } else {
        arg = sqrt(arg);
        double f = 1.0/arg * erf(arg) *sqrt(M_PI)/2.0;
        return f;
    }
}

double BoysFunction::result(int level) const
{
    return m_results(level);
}
