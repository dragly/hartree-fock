#ifndef UNRESTRICTEDHARTREEFOCKSOLVER_H
#define UNRESTRICTEDHARTREEFOCKSOLVER_H

#include <hartreefocksolver.h>

class UnrestrictedHartreeFockSolver : public HartreeFockSolver
{
public:
    explicit UnrestrictedHartreeFockSolver(ElectronSystem *system);

    virtual void reset();
    virtual void advance();

    virtual double energy();
    virtual void solve();
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
};

#endif // UNRESTRICTEDHARTREEFOCKSOLVER_H
