#include "boysfunctionintermediate.h"

#include <math/boysfunction.h>
#include <string>

using namespace std;

BoysFunctionIntermediate::BoysFunctionIntermediate(int levelMax, int nValues, double limitMin, double limitMax, int nIntegralValues) :
    m_nValues(nValues),
    m_levelMax(levelMax),
    m_limitMin(limitMin),
    m_limitMax(limitMax),
    m_nIntegralValues(nIntegralValues)
{
    m_dx = (limitMax - limitMin) / (nValues - 1);
    updateResults();
}

void BoysFunctionIntermediate::updateResults() {
//    cout << "BoysFunctionIntermediate::updateResults(): starting" << endl;
    //    int nValues = m_nValues;
    //    vec previousResults = zeros(m_levelMax + 1 + (m_taylorExpansionOrder + 1) * nValues);
    //    vec nextResults = zeros(m_levelMax + 1 + (m_taylorExpansionOrder + 1) * nValues);
    //    // Calculate for all n in first x (i = 0)
    //    for(int n = 0; n < m_levelMax + 1 + m_taylorExpansionOrder * nValues + 1; n++) {
    //        uint i = 0;
    //        double x = m_args[i];
    //        previousResults(n) = BoysFunction::calculateTaylorExpansion(x, n);
    //    }

    //    for(int i = 0; i < m_nValues - 1; i++) {
    //        double dx = m_dx;
    //        for(int n = 0; n < m_levelMax + 1 + m_taylorExpansionOrder * (nValues - i); n++) {
    //            double sumResult = 0.0;
    //            for(int k = 0; k < m_taylorExpansionOrder; k++) {
    //                double F_n_k_x = previousResults(n + k);
    //                double kFac = BoysFunction::factorial(k);

    //                sumResult += F_n_k_x * pow(-dx, k) / kFac;
    //            }
    //            nextResults(n) = sumResult;
    //        }
    //        for(int n = 0; n < m_levelMax + 1; n++) {
    //            m_results(i+1, n) = nextResults(n);
    //        }
    //        previousResults = nextResults;
    //    }
    // Try to load results from file
    stringstream fileNameStream;
    fileNameStream << "boys_function_data"
//                   << "_lmax_" << m_levelMax
                   << "_nx" << m_nValues
                   << "_limmin_" << m_limitMin
                   << "_limmax_" << m_limitMax
                   << "_nt" << m_nIntegralValues
                   << ".arma";
    bool allGood = m_results.load(fileNameStream.str());
    if(allGood) {
        if(m_results.n_cols < uint(m_levelMax)) {
            allGood = false;
        }
    }
    if(!allGood) {
        cout << "BoysFunctionIntermediate::updateResults(): Boys function data file not found. Generating now." << endl;
        cout << "BoysFunctionIntermediate::updateResults(): Filename: " << fileNameStream.str() << endl;
        m_results = zeros(m_nValues, m_levelMax + 1 + m_taylorExpansionOrder + 1); // + 6 for the Taylor expansions
        // Could not load results from file. Generate results and write file instead.
        for(uint i = 0; i < m_nValues; i++) {
            double x = i * m_dx;
            cout << "x = " << x << endl;
            for(uint n = 0; n < uint(m_levelMax + 1) + m_taylorExpansionOrder + 1; n++) {
                m_results(i,n) = directIntegral(x, n);
            }
        }
        cout << "BoysFunctionIntermediate::updateResults(): Generating results done. Writing to file " << fileNameStream.str() << endl;
        m_results.save(fileNameStream.str());
    }
}

double BoysFunctionIntermediate::result(double arg, int n) const {
    // Linear interpolation
    double dx = m_dx;
    int closestIndex = (arg - m_limitMin + dx / 2.0) / dx; // + dx / 2.0 always gives the closest index
//    double closestArg = m_args[closestIndex];
    double closestArg = closestIndex * dx;
    double difference = arg - closestArg;
    double sumResult = 0.0;
    for(int k = 0; k < m_taylorExpansionOrder + 1; k++) {
        double F_n_k_x = m_results(closestIndex, n + k);
        double kFac = BoysFunction::factorial(k);
        sumResult += F_n_k_x * fastPow(-difference, k) / kFac;
    }
    return sumResult;
}

double BoysFunctionIntermediate::directIntegral(double x, double n) const {
    double t0 = 0.0;
    double t1 = 1.0;
    double nt = m_nIntegralValues;
    double dt = (t1 - t0) / nt;
    double integralResult = 0;
    for(int i = 0; i < nt; i++) {
        double t = t0 + dt * i;
        integralResult += 2 * integrand(x, t + dt, n);
    }
    integralResult -= integrand(x, 0, n);
    integralResult -= integrand(x, (nt - 1)*dt, n);
    integralResult *= dt / 2.0;
    return integralResult;
}

double BoysFunctionIntermediate::fastPow(double x, int a) const {
    double powResult = 1;
    for(int i = 0; i < a; i++) {
        powResult *= x;
    }
    return powResult;
}

double BoysFunctionIntermediate::integrand(double x, double t, double n) const {
    return exp(-x*t*t) * fastPow(t, 2*n);
}
