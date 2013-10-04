#include <../tests/helium-tests.cpp>

#include <iostream>
#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <src/hartreesolver.h>
#include <src/hartreefocksolver.h>
#include <src/electronsystems/helium/heliumhartree.h>
#include <src/electronsystems/hydrogen/hydrogenmolecule.h>
#include <src/basisfunctions/gaussiantypeorbital.h>
#include <armadillo>
#include <iostream>

#include <fstream>

using namespace std;
using namespace std;
using namespace arma;

SUITE(Complete) {
    TEST(HydrogenOverlapIntegral)  {
        HydrogenMolecule basisFunction;
        double value = basisFunction.overlapIntegral(6,2);
        cout << "value: " << value << endl;
    }

    TEST(HydrogenElectronInteractionIntegral)  {
        HydrogenMolecule basisFunction;
        double value = basisFunction.electronInteractionIntegral(1,5,7,2);
        cout << "value: " << value << endl;
    }

    TEST(HydrogenNuclearAttractionIntegral)  {
        HydrogenMolecule basisFunction;
        double value = basisFunction.nuclearAttractionIntegral(3,7);
        cout << "value: " << value << endl;
    }

    TEST(HydrogenAdvanceMany) {
        HydrogenMolecule basisFunction(1.0);
        HartreeSolver solver(&basisFunction);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        cout << "Energy: " << solver.energy() << endl;
        cout << "Energy with repulsion: " << solver.energy() << endl;
    }

    TEST(HydrogenAdvanceManyHF) {
        cout << "Hydrogen with HF:" << endl;
        HydrogenMolecule basisFunction(1.0);
        HartreeFockSolver solver(&basisFunction);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        cout << "Energy: " << solver.energy() << endl;
        cout << "Energy with repulsion: " << solver.energy() << endl;
    }

    TEST(MultiHydrogenAdvanceManyHF) {
#include <src/electronsystems/hydrogen/multihydrogen.h>
        cout << "MultiHydrogen with HF:" << endl;
        mat nucleiPositions {0,0,0,
                            1,0,0,
                            0,1,0};
        nucleiPositions.reshape(3,3);
        nucleiPositions = nucleiPositions.t();
        MultiHydrogen basisFunction(nucleiPositions);
        HartreeFockSolver solver(&basisFunction);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        cout << "Energy: " << solver.energy() << endl;
        cout << "Energy with repulsion: " << solver.energy() << endl;
    }

    TEST(MultiHydrogenAdvanceManyHFPlot) {
#include <src/electronsystems/hydrogen/multihydrogen.h>
        cout << "MultiHydrogen with HF:" << endl;

        rowvec xVec = linspace(-4, 4, 207);
        rowvec yVec = linspace(-4, 4, 207);

        ofstream energyFile;
        ofstream xFile;
        ofstream yFile;
        energyFile.open("energymap.dat");
        xFile.open("energymapx.dat");
        yFile.open("energymapy.dat");

        for(uint i = 0; i < xVec.n_rows; i++) {
            double x = xVec(i);
            cout << "x: " << x << endl;
            for(uint j = 0; j < yVec.n_rows; j++) {
                double y = yVec(j);
                double R = 1.4;
                double Rx = R / 2;
                double Ry = sqrt(0.75)*R / 2;
                mat nucleiPositions {-Rx, -Ry, 0.0,
                                      Rx, -Ry, 0.0,
                                      x,   y,   0.0};
                nucleiPositions.reshape(3,3);
                nucleiPositions = nucleiPositions.t();

                bool skipCalculation = false;
                for(uint k = 0; k < nucleiPositions.n_rows; k++) {
                    for(uint l = k + 1; l < nucleiPositions.n_rows; l++) {
                        rowvec rDiff = nucleiPositions.row(k) - nucleiPositions.row(l);
                        if(dot(rDiff, rDiff) < 1e-4) {
                            skipCalculation = true;
                        }
                    }
                }
                if(skipCalculation) {
                    energyFile << 0 << " ";
                    xFile << x << " ";
                    yFile << y << " ";
                    continue;
                }

                MultiHydrogen basisFunction(nucleiPositions);
                HartreeFockSolver solver(&basisFunction);
                double previousValue = 0;
                bool tresholdReached = false;
                for(int i = 0; i < 30 && !tresholdReached; i++) {
                    solver.advance();
                    if(i > 0) {
                        double energyDiff = fabs(solver.energy() - previousValue);
                        if(energyDiff < 1e-3) {
                            tresholdReached = true;
                        }
                    }
                    previousValue = solver.energy();
                }
                energyFile << solver.energy() << " ";
                xFile << x << " ";
                yFile << y << " ";
            }
            energyFile << endl;
            xFile << endl;
            yFile << endl;
        }
        energyFile.close();
        xFile.close();
        yFile.close();
    }
}

SUITE(Old) {
    TEST(HydrogenPlot) {
        rowvec distances = linspace(1.0, 4.0, 100);
        std::ofstream energyFile;
        energyFile.open("energies.dat");
        for(double distance : distances) {
            cout << distance << endl;
            HydrogenMolecule basisFunction(distance);
            HartreeSolver solver(&basisFunction);
            solver.reset();
            for(int i = 0; i < 100; i++) {
                solver.advance();
            }
            double energy = solver.energy() + basisFunction.additionalEnergyTerms();
            energyFile << distance << " " << energy << "\n";
        }
        energyFile.close();
    }
}

SUITE(Development) {
    TEST(GaussianTypeOrbital) {
        GaussianTypeOrbitalIntegrator integrator;
        rowvec posA = {1,2,3};
        integrator.setCorePositionA(posA);
        integrator.setOrbitalExponentA(0.1);
        integrator.setOrbitalExponentB(0.2);
        for(const urowvec& combination : integrator.combinationsA()) {
            cout << combination << endl;
        }
        integrator.reset();
    }
}

int main()
{
    int result = 0;
#ifdef RUN_DEVELOPMENT_TESTS
    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);
    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Development", UnitTest::True(), 0);
#else
    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);
    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Complete", UnitTest::True(), 0);
#endif
    return result;
}

