#ifndef UNRESTRICTEDHARTREEFOCKSOLVER_H
#define UNRESTRICTEDHARTREEFOCKSOLVER_H

#include "solvers/hartreefocksolver.h"
#include <vector>

class UnrestrictedHartreeFockSolver : public HartreeFockSolver
{
public:
    explicit UnrestrictedHartreeFockSolver(ElectronSystem *system);

    virtual void setup();
    virtual void advance();
    virtual double energy();
    virtual void solve();

    const mat &coeffcientMatrixUp() const;
    const mat &coeffcientMatrixDown() const;

    const mat &densityMatrixUp() const;
    const mat &densityMatrixDown() const;

    void setInitialCoefficientMatrices(const mat& up, const mat &down);

    bool isDiisEnabled() const;
    void setDiisEnabled(bool isDiisEnabled);

    int diisSampleCount() const;
    void setDiisSampleCount(int diisSampleCount);

    int diisStartingIteration() const;
    void setDiisStartingIteration(int diisStartingIteration);

private:
    void resetCoefficientMatrices();
    void setupFockMatrices();
    void setupDensityMatrices();
    void resetFockMatrices();
    void randomizeCoefficientMatrices();
    void calculateEnergy();
    void performDIIS();

    mat m_fockMatrixUp;
    mat m_fockMatrixDown;

    mat m_coefficientMatrixUp;
    mat m_coefficientMatrixDown;

    mat m_initialCoefficientMatrixUp;
    mat m_initialCoefficientMatrixDown;

    mat m_densityMatrixUp;
    mat m_densityMatrixDown;

    vec m_fockEnergiesUp;
    vec m_fockEnergiesDown;

    double m_energyUHF;
    bool m_initialCoefficientMatricesSetManually;

    int m_diisSampleCount;
    int m_diisStartingIteration;

    std::vector<mat> m_errorsU, m_fockMatricesU;
    std::vector<mat> m_errorsD, m_fockMatricesD;

    bool m_isDiisEnabled;
};

#endif // UNRESTRICTEDHARTREEFOCKSOLVER_H
