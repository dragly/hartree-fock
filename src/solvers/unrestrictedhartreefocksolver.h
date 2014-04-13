#ifndef UNRESTRICTEDHARTREEFOCKSOLVER_H
#define UNRESTRICTEDHARTREEFOCKSOLVER_H

#include "solvers/hartreefocksolver.h"

class UnrestrictedHartreeFockSolver : public HartreeFockSolver
{
public:
    explicit UnrestrictedHartreeFockSolver(ElectronSystem *system);

    virtual void setup();
    virtual void advance();
    virtual double energy();
    virtual void solve();

    const mat &coeffcientMatrixUp();
    const mat &coeffcientMatrixDown();
private:
    void resetCoefficientMatrices();

    mat m_fockMatrixUp;
    mat m_fockMatrixDown;

    mat m_coefficientMatrixUp;
    mat m_coefficientMatrixDown;

    mat m_densityMatrixUp;
    mat m_densityMatrixDown;

    vec m_fockEnergiesUp;
    vec m_fockEnergiesDown;

    double m_energyUHF;

    void setupFockMatrices();
    void setupDensityMatrices();
    void resetFockMatrices();
};

#endif // UNRESTRICTEDHARTREEFOCKSOLVER_H
