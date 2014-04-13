#ifndef RESTRICTEDHARTREEFOCKSOLVER_H
#define RESTRICTEDHARTREEFOCKSOLVER_H

#include "solvers/hartreefocksolver.h"
class ElectronSystem;

class RestrictedHartreeFockSolver : public HartreeFockSolver
{
public:
    RestrictedHartreeFockSolver(ElectronSystem *electronSystem);

    virtual ~RestrictedHartreeFockSolver();

    virtual void setup();
    virtual void advance();
    virtual void solve();

    virtual double energy();
    const mat &densityMatrix() const;
    const mat& coefficientMatrix() const;
private:
    void setupFockMatrix();
    void setupDensityMatrix();
    void resetCoefficientMatrix();
    void resetFockMatrix();
    double coupledMatrixTilde(int p, int q, int r, int s);

    // Matrices
    mat m_fockMatrix;
    mat m_densityMatrix;
    mat m_coefficientMatrix;
    vec m_fockEnergies;

    // Variables
    double m_energy = 0;
};

#endif // RESTRICTEDHARTREEFOCKSOLVER_H
