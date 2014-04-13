#include "hermiteintegral.h"

#include "math/vector3.h"
#include <math/boysfunction.h>

HermiteIntegral::HermiteIntegral(int dimensionMax) :
    m_alpha(0),
    m_A(Vector3::createZeros()),
    m_dimensionMax(dimensionMax)
{
    reset(dimensionMax);
}

void HermiteIntegral::reset(int dimension)
{
    m_dimensionMax = dimension;
    int tMax = m_dimensionMax;
    int nMax = m_dimensionMax;
    m_R.reset();
    // nMax + 1 because looking up element n = nMax requires a size of nMax + 1 to be available
    m_R.set_size(nMax + 1);
    for(int n = 0; n <= nMax; n++) {
        m_R(n) = zeros(tMax + 1, tMax + 1, tMax + 1);
    }
}

void HermiteIntegral::set(double alpha, const Vector3 &A,
                          int t, int u, int v)
{
    m_alpha = alpha;
    m_A = A;
    setupR(t, u, v);
}

void HermiteIntegral::setupR(int tin, int uin, int vin) {
    const Vector3 &PC = m_A;
    double p = m_alpha;
    double boysArg = p * dot(m_A, m_A);
    int tuvSumMax = tin + uin + vin;
    int nMax = tin + uin + vin;
    boysFunction.set(boysArg, nMax);
    // Calculate R0_tuv
    double powResult = 1;
    m_R(0)(0,0,0) = boysFunction.result(0);
    for(int n = 1; n <= nMax; n++) {
        powResult *= -2*p;
        m_R(n)(0,0,0) = powResult * boysFunction.result(n);
    }
    // Iterate over all elements with t + u + v = tuvSum (slide 36)
    for(int tuvSum = 1; tuvSum <= tuvSumMax; tuvSum++) {
        for(int n = 0; n <= nMax - tuvSum; n++) {
            for(int t = 0; t <= tin; t++) {
                for(int u = 0; u <= uin; u++) {
                    for(int v = 0; v <= vin; v++) {
                        if(t + u + v != tuvSum || t + u + v == 0) {
                            continue;
                        }
                        int largestElement = max(t,max(u,v));
                        // Needed indices
                        int t2 = t;
                        int u2 = u;
                        int v2 = v;
                        int t3 = t;
                        int u3 = u;
                        int v3 = v;
                        // Resulting factors
                        int factor2 = t;
                        double factor3 = PC(0);
                        // Choose indices based on method (slide 34)
                        if(largestElement == t) {
                            t2 = t - 2;
                            t3 = t - 1;
                            factor2 = t - 1;
                            factor3 = PC(0);
                        } else if(largestElement == u) {
                            u2 = u - 2;
                            u3 = u - 1;
                            factor2 = u - 1;
                            factor3 = PC(1);
                        } else {
                            v2 = v - 2;
                            v3 = v - 1;
                            factor2 = v - 1;
                            factor3 = PC(2);
                        }
                        double currentValue = 0;
                        if(t2 >= 0 && u2 >= 0 && v2 >= 0) {
                            currentValue += factor2 * m_R(n+1)(t2, u2, v2);
                        }
                        if(t3 >= 0 && u3 >= 0 && v3 >= 0) {
                            currentValue += factor3 * m_R(n+1)(t3, u3, v3);
                        }
                        m_R(n)(t, u, v) = currentValue;
                    }
                }
            }
        }
    }
}
