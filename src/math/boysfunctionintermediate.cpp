#include "boysfunctionintermediate.h"

#include <math/boysfunction.h>
#include <string>

using namespace std;

/*!
 * \class BoysFunctionIntermediate
 * \brief Helper class for fast evaluations in BoysFunction
 */

BoysFunctionIntermediate::BoysFunctionIntermediate(int levelMax, int nValues, double limitMin, double limitMax, int nIntegralValues) :
    m_nValues(nValues),
    m_levelMax(levelMax),
    m_neededLevelMax(m_levelMax + m_taylorExpansionOrder + 1),
    m_limitMin(limitMin),
    m_limitMax(limitMax),
    m_nIntegralValues(nIntegralValues)
{
    factorialInverseTable[0] = 1.0;
    factorialInverseTable[1] = 1.0;
    factorialInverseTable[2] = 0.5;
    factorialInverseTable[3] = 0.166666666667;
    factorialInverseTable[4] = 0.0416666666667;
    factorialInverseTable[5] = 0.00833333333333;
    factorialInverseTable[6] = 0.00138888888889;
    factorialInverseTable[7] = 0.000198412698413;
    factorialInverseTable[8] = 2.48015873016e-05;
    factorialInverseTable[9] = 2.7557319224e-06;
    factorialInverseTable[10] = 2.7557319224e-07;
    factorialInverseTable[11] = 2.50521083854e-08;
    factorialInverseTable[12] = 2.08767569879e-09;
    factorialInverseTable[13] = 1.60590438368e-10;
    factorialInverseTable[14] = 1.14707455977e-11;
    factorialInverseTable[15] = 7.64716373182e-13;
    factorialInverseTable[16] = 4.77947733239e-14;
    factorialInverseTable[17] = 2.81145725435e-15;
    factorialInverseTable[18] = 1.56192069686e-16;
    factorialInverseTable[19] = 8.22063524662e-18;
    factorialInverseTable[20] = 4.11031762331e-19;
    factorialInverseTable[21] = 1.95729410634e-20;
    factorialInverseTable[22] = 8.89679139245e-22;
    factorialInverseTable[23] = 3.86817017063e-23;
    factorialInverseTable[24] = 1.6117375711e-24;
    factorialInverseTable[25] = 6.44695028438e-26;
    factorialInverseTable[26] = 2.47959626322e-27;
    factorialInverseTable[27] = 9.1836898638e-29;
    factorialInverseTable[28] = 3.27988923707e-30;
    factorialInverseTable[29] = 1.13099628864e-31;
    factorialInverseTable[30] = 3.76998762882e-33;
    factorialInverseTable[31] = 1.21612504155e-34;
    factorialInverseTable[32] = 3.80039075485e-36;
    factorialInverseTable[33] = 1.15163356208e-37;
    factorialInverseTable[34] = 3.38715753552e-39;
    factorialInverseTable[35] = 9.67759295863e-41;
    factorialInverseTable[36] = 2.68822026629e-42;
    factorialInverseTable[37] = 7.26546017915e-44;
    factorialInverseTable[38] = 1.91196320504e-45;
    factorialInverseTable[39] = 4.90246975651e-47;

    m_dx = (limitMax - limitMin) / (nValues - 1);
    m_dxInverse = 1.0 / m_dx;
    updateResults();
}

void BoysFunctionIntermediate::updateResults() {
    // Try to load results from file
    stringstream fileNameStream;
    fileNameStream << "boys_function_data"
//                   << "_lmax_" << m_levelMax
                   << "_nx" << m_nValues
                   << "_limmin_" << m_limitMin
                   << "_limmax_" << m_limitMax
                   << "_nt" << m_nIntegralValues
                   << ".arma";
    bool allGood = m_results.load("boys_tabulated.dat");
    if(!allGood) {
        cout << "BoysFunctionIntermediate::updateResults(): Boys function data file not found. Generating now." << endl;
    }
    if(allGood) {
        if(m_results.n_cols < uint(m_neededLevelMax)) {
            cout << "BoysFunctionIntermediate::updateResults(): Boys function file has " << m_results.n_cols << " cols, "
                 << "but we need " << m_neededLevelMax << endl;
            allGood = false;
        }
    }
    if(!allGood) {
        cout << "BoysFunctionIntermediate::updateResults(): Filename: " << fileNameStream.str() << endl;
        cout << "BoysFunctionIntermediate::updateResults(): Max level: " << m_levelMax << endl;
        m_results = zeros(m_nValues, m_neededLevelMax); // + 6 for the Taylor expansions
        // Could not load results from file. Generate results and write file instead.
        for(uint i = 0; i < m_nValues; i++) {
            double x = i * m_dx;
            cout << "x = " << x << endl;
            for(uint n = 0; n < uint(m_neededLevelMax); n++) {
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
    double dxI = m_dxInverse;
    int closestIndex = (arg - m_limitMin + dx * 0.5) * dxI; // + dx / 2.0 always gives the closest index
    double closestArg = closestIndex * dx;
    double difference = arg - closestArg;
    double sumResult = m_results(closestIndex, n);
    double powResult = 1;
    for(int k = 1; k < m_taylorExpansionOrder + 1; k++) {
        double F_n_k_x = m_results(closestIndex, n + k);
        powResult *= -difference;
        sumResult += F_n_k_x * powResult * factorialInverseTable[k];
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
