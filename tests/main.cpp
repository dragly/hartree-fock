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
#include <src/basisfunctions/gaussiantypeorbitalintegrator.h>
#include <src/math/boysfunction.h>
#include <armadillo>
#include <iostream>
#include <src/math/boysfunctionintermediate.h>

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

class MyException : public exception {
public:
    MyException(string errorString) : exception() {
        m_errorString = errorString;
    }

    virtual const char* what() const throw() {
        return m_errorString.c_str();
    }

private:
    string m_errorString;
};

SUITE(Development) {
    TEST(GaussianTypeOrbitalIntegration) {
        GaussianTypeOrbitalIntegrator integrator;
        rowvec posA = {1.2,2.3,3.4};
        integrator.setCorePositionA(posA);
        rowvec posB = {-1.3,1.4,-2.4};
        integrator.setCorePositionB(posB);
        integrator.setExponentA(0.2);
        integrator.setExponentB(0.3);
        integrator.setMaxAngularMomentumA(4);
        integrator.setMaxAngularMomentumB(4);
//        for(const urowvec& combination : integrator.combinationsA()) {
//            cout << combination << endl;
//        }
        integrator.reset();
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,0), 0.119172363580852, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,1), 0.276479883507577, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,2), 0.760605693318432, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,0), 0.0429020508891068, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,1), 0.0995327580627279, 0.00001);

        CHECK_CLOSE(-0.0967870268058250, integrator.kineticIntegral(0,0,0,0,0,0), 0.00001);
        CHECK_CLOSE(-0.158190730147696, integrator.kineticIntegral(0,0,0,0,0,1), 0.00001);
        CHECK_CLOSE(-0.0245468374367114, integrator.kineticIntegral(0,0,0,0,1,0), 0.00001);
        CHECK_CLOSE(-0.0330608009181156, integrator.kineticIntegral(0,0,0,0,1,1), 0.00001);
        CHECK_CLOSE(-0.0681856595464203, integrator.kineticIntegral(0,0,0,1,0,0), 0.00001);
        CHECK_CLOSE(-0.0918355581058773, integrator.kineticIntegral(0,0,0,1,0,1), 0.00001);
        CHECK_CLOSE(-0.0142503452233260, integrator.kineticIntegral(0,0,0,1,1,0), 0.00001);
        CHECK_CLOSE(-0.00917293898306074, integrator.kineticIntegral(0,0,0,1,1,1), 0.00001);
        CHECK_CLOSE(0.251402082662018, integrator.kineticIntegral(0,0,1,0,0,1), 0.00001);
        CHECK_CLOSE(0.0495912013771734, integrator.kineticIntegral(0,0,1,0,1,0), 0.00001);
        CHECK_CLOSE(0.0176714824377215, integrator.kineticIntegral(0,0,1,0,1,1), 0.00001);
        CHECK_CLOSE(0.137753337158815, integrator.kineticIntegral(0,0,1,1,0,0), 0.00001);
        CHECK_CLOSE(0.0490874512158928, integrator.kineticIntegral(0,0,1,1,0,1), 0.00001);
        CHECK_CLOSE(0.0137594084745900, integrator.kineticIntegral(0,0,1,1,1,0), 0.00001);
        CHECK_CLOSE(-0.0551617848828804, integrator.kineticIntegral(0,0,1,1,1,1), 0.00001);
        CHECK_CLOSE(-0.0604904731258242, integrator.kineticIntegral(0,1,0,0,1,0), 0.00001);
        CHECK_CLOSE(-0.0868821710550227, integrator.kineticIntegral(0,1,0,0,1,1), 0.00001);
        CHECK_CLOSE(0.0213755178349884, integrator.kineticIntegral(0,1,0,1,0,0), 0.00001);
        CHECK_CLOSE(0.0137594084745909, integrator.kineticIntegral(0,1,0,1,0,1), 0.00001);
        CHECK_CLOSE(-0.0374492116616455, integrator.kineticIntegral(0,1,0,1,1,0), 0.00001);
        CHECK_CLOSE(-0.0334264444581330, integrator.kineticIntegral(0,1,0,1,1,1), 0.00001);
        CHECK_CLOSE(0.0788748150526478, integrator.kineticIntegral(0,1,1,0,1,1), 0.00001);
        CHECK_CLOSE(-0.0206391127118871, integrator.kineticIntegral(0,1,1,1,0,0), 0.00001);
        CHECK_CLOSE(0.0827426773243232, integrator.kineticIntegral(0,1,1,1,0,1), 0.00001);

        cout << "Tall: " << integrator.kineticIntegral(2,0,0,0,0,0) << endl;
        cout << "Tall: " << integrator.kineticIntegral(0,0,0,1,0,0) << endl;
        cout << "Tall: " << integrator.kineticIntegral(0,1,0,0,0,0) << endl;
        cout << "Tall: " << integrator.kineticIntegral(0,0,0,0,1,0) << endl;
        cout << "Tall: " << integrator.kineticIntegral(0,0,1,0,0,0) << endl;
        cout << "Tall: " << integrator.kineticIntegral(0,0,0,0,0,1) << endl;
    }

    TEST(BoysTest) {
        int maxLevel = 21;
        int testLevel = 10;
        vec x = linspace(0.001, 35, 1000);
        BoysFunctionIntermediate intermediate(maxLevel, 1000, 0, 30);
        for(int n = testLevel; n < maxLevel + 1; n++) {
            cout << "Building level " << n << endl;
            stringstream fileName;
            fileName << "results" << n << ".txt";
            ofstream results(fileName.str());
            for(uint i = 0; i < x.n_elem; i++) {
                BoysFunction boys(x(i), n, &intermediate);
                results << x(i) << " " << boys.result(testLevel) << endl;
            }
            results.close();
        }
//        BoysFunction boys(0,0);
//        cout << "Result = " << boys.calculateAsymptopticForm(12, 11) << endl;
    }

    TEST(BoysIntermediateTest) {
        BoysFunctionIntermediate boysIntermediate(20, 1000);
        boysIntermediate.updateResults();
        cout << "Boys, boys, boys = " << boysIntermediate.result(6, 10) << endl;
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

