#include "boysfunctionintermediate.h"

#include <src/math/boysfunction.h>

BoysFunctionIntermediate::BoysFunctionIntermediate(int levelMax, int nValues, double limitMin, double limitMax) :
    m_nValues(nValues),
    m_levelMax(levelMax),
    m_limitMin(limitMin),
    m_limitMax(limitMax)
{
    m_results = zeros(nValues, levelMax + 1); // + 6 for the Taylor expansion
    m_args = linspace(m_limitMin, m_limitMax, nValues);
    m_dx = m_args(1) - m_args(0);
    updateResults();
}

void BoysFunctionIntermediate::updateResults() {
    cout << "BoysFunctionIntermediate::updateResults(): starting" << endl;
    int nValues = m_nValues;
    vec previousResults = zeros(m_levelMax + 1 + (m_taylorExpansionOrder + 1) * nValues);
    vec nextResults = zeros(m_levelMax + 1 + (m_taylorExpansionOrder + 1) * nValues);
    // Calculate for all n in first x (i = 0)
    for(int n = 0; n < m_levelMax + 1 + m_taylorExpansionOrder * nValues + 1; n++) {
        uint i = 0;
        double x = m_args[i];
        previousResults(n) = BoysFunction::calculateTaylorExpansion(x, n);
    }

    for(int i = 0; i < m_nValues - 1; i++) {
        double dx = m_dx;
        for(int n = 0; n < m_levelMax + 1 + m_taylorExpansionOrder * (nValues - i); n++) {
            double sumResult = 0.0;
            for(int k = 0; k < m_taylorExpansionOrder; k++) {
                double F_n_k_x = previousResults(n + k);
                double kFac = BoysFunction::factorial(k);

                sumResult += F_n_k_x * pow(-dx, k) / kFac;
            }
            nextResults(n) = sumResult;
        }
        for(int n = 0; n < m_levelMax + 1; n++) {
            m_results(i+1, n) = nextResults(n);
        }
        previousResults = nextResults;
    }
    cout << "BoysFunctionIntermediate::updateResults(): done" << endl;
}

double BoysFunctionIntermediate::result(double arg, int n) const {
    // Linear interpolation
    double dx = m_dx;
    int index0 = (arg - m_limitMin) / dx;
    int index1 = index0 + 1;
    double arg0 = m_args[index0];
    double result0 = m_results(index0, n);
    double result1 = m_results(index1, n);
    double weight = (arg - arg0) / dx;
    double returnResult = (1-weight) * result0 + weight * result1;
    return returnResult;
}
