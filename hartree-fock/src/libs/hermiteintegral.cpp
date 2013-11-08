#include "hermiteintegral.h"

#include <math/boysfunction.h>
#include <math/boysfunctionintermediate.h>

HermiteIntegral::HermiteIntegral(int dimension) :
    m_alpha(0),
    m_A(zeros<rowvec>(3)),
    m_dimension(dimension)
{
    reset(dimension);
}

HermiteIntegral::HermiteIntegral(double alpha, const rowvec &A, int dimension, bool setupImmediately) :
    HermiteIntegral(dimension)
{
    set(alpha, A, setupImmediately);
}

void HermiteIntegral::reset(int dimension)
{
    m_dimension = dimension;
    int tMax = m_dimension;
    int nMax = m_dimension;
    m_R.reset();
    m_R.set_size(nMax + 1);
    // Initialize and allocate the cubes for R
    for(int n = 0; n < nMax + 1; n++) {
        m_R(n) = zeros(tMax + 1, tMax + 1, tMax + 1);
    }
}

void HermiteIntegral::set(double alpha, const rowvec &A, bool setupImmediately)
{
    m_alpha = alpha;
    m_A = A;
    if(setupImmediately) {
        setupR();
    }
}

void HermiteIntegral::setupR() {
//    const rowvec &PA = m_centerOfMassDiffA;
    const rowvec &PC = m_A;
    double p = m_alpha;
    double boysArg = p * dot(m_A, m_A);
//    int lMax = 2;
    int tMax = m_dimension;
    int tuvSumMax = m_dimension;
    int nMax = m_dimension;
    BoysFunctionIntermediate &boysFunctionIntermediate = BoysFunctionIntermediate::getInstance();
    BoysFunction boysFunction(boysArg, nMax + 1, &boysFunctionIntermediate);
//    vector<vector<int>> calculatedIndices;
    // Calculate R0_tuv
    for(int n = 0; n < nMax + 1; n++) {
        m_R(n)(0,0,0) = pow(-2*p, n) * boysFunction.result(n);
//        cout << m_R[n](0,0,0) << endl;
//        vector<int> indices = {n,0,0,0};
//        calculatedIndices.push_back(indices);
    }
    for(int tuvSum = 1; tuvSum < tuvSumMax + 1; tuvSum++) {
//        cout << "All elements for which t + u + v = " << tuvSum << endl;
        for(int n = 0; n < nMax - tuvSum; n++) {
            for(int t = 0; t < tMax; t++) {
                for(int u = 0; u < tMax; u++) {
                    for(int v = 0; v < tMax; v++) {
                        if(t + u + v != tuvSum || t + u + v == 0) {
                            continue;
                        }
//                        if(n == 0 && t == 0 && u == 1 && v == 0) {
//                            cout << "woot" << endl;
//                        }
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
//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "t biggest" << endl;
//                            }
                            t2 = t - 2;
                            t3 = t - 1;
                            factor2 = t - 1;
                            factor3 = PC(0);
                        } else if(largestElement == u) {
//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "u biggest" << endl;
//                            }
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
//                        vector<int> currentIndices = {n, tRec1, uRec1, vRec1};
                        if(t2 >= 0 && u2 >= 0 && v2 >= 0) {

//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "Adding factor2 " << endl;
//                            }
                            currentValue += factor2 * m_R[n+1](t2, u2, v2);
//                            vector<int> neededIndices = {n+1,tRec2,uRec2,vRec2};
//                            if(!hasIndices(neededIndices, calculatedIndices)) {
//                                cout << "Target: " <<  currentIndices[0] << " " << currentIndices[1] << " " << currentIndices[2] << " " << currentIndices[3] << endl;
//                                cout << "Needs: " <<  neededIndices[0] << " " << neededIndices[1] << " " << neededIndices[2] << " " << neededIndices[3] << endl;
//                                break;
//                            }
                        }
                        if(t3 >= 0 && u3 >= 0 && v3 >= 0) {
//                            if(n == 0 && t == 0 && u == 1 && v == 0) {
//                                cout << "Adding factor3 " << endl;
//                                cout << n + 1 << " " << t3 << " " << u3 << " " << v3 << endl;
//                                cout << factor3 << endl;
//                                cout << PC(1) << endl;
//                            }
                            currentValue += factor3 * m_R[n+1](t3, u3, v3);
//                            vector<int> neededIndices2 = {n+1,tRec3,uRec3,vRec3};
//                            if(!hasIndices(neededIndices2, calculatedIndices)) {
//                                cout << "Target: " <<  currentIndices[0] << " " << currentIndices[1] << " " << currentIndices[2] << " " << currentIndices[3] << endl;
//                                cout << "Needs: " <<  neededIndices2[0] << " " << neededIndices2[1] << " " << neededIndices2[2] << " " << neededIndices2[3] << endl;
//                                break;
//                            }
                        }

//                        if(n == 0 && t == 1 && u == 0 && v == 0) {
//                            cout << "currentValue " << currentValue << endl;
//                        }
                        m_R[n](t, u, v) = currentValue;
//                        cout << "Calculated: " <<  currentIndices[0] << " " << currentIndices[1] << " " << currentIndices[2] << " " << currentIndices[3] << endl;
//                        calculatedIndices.push_back(currentIndices);
                    }
                }
            }
        }
    }
//    cout << banana.find(elementa) << endl;

//    double F0 = boysFunction.result(0);
//    double F1 = boysFunction.result(1);
//    double F2 = boysFunction.result(2);
//    double R0000 = F0;
//    double R1000 = -2 * p * F1;
//    double R2000 = (-2 * p)*(-2*p) * F2;
//    cout << "R0000 = " << R0000 << endl;
//    cout << "R1000 = " << R1000 << endl;
//    cout << "R2000 = " << R2000 << endl;
//    cube *E = m_E;
//    cout << "Result = " << 2 * M_PI / p * R0000 * E[0](0,0,0) * E[1](0,0,0) * E[2](0,0,0) << endl;
}
