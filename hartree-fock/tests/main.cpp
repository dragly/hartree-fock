//#include <../tests/helium-tests.cpp>

#include <iostream>
#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <hartreesolver.h>
#include <hartreefocksolver.h>
#include <electronsystems/helium/heliumhartree.h>
#include <electronsystems/hydrogen/hydrogenmolecule.h>
#include <electronsystems/hydrogen/multihydrogen.h>
#include <basisfunctions/gaussiantypeorbital.h>
#include <basisfunctions/gaussiantypeorbitalintegrator.h>
#include <math/boysfunction.h>
#include <basisfunctions/gaussiantypeorbitalintegrator.h>
#include <math/gaussiantypeoverlapintegral.h>
#include <math/gaussiantypekineticintegral.h>
#include <math/gaussiantypecoloumbattractionintegral.h>
#include <math/gaussiantypeelectroninteractionintegral.h>
#include <armadillo>
#include <iostream>
#include <math/boysfunctionintermediate.h>

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

SUITE(Development2) {
    TEST(BoysPlotTest) {
        int maxLevel = 21;
        int testLevel = 0;
        vec x = linspace(0.0, 60, 1000);
        BoysFunctionIntermediate intermediate(maxLevel, 1000, 0, 50, 5e3);
        for(int n = testLevel; n < maxLevel + 1; n++) {
            stringstream fileName;
            fileName << "results" << n << ".txt";
            ofstream results(fileName.str());
            for(uint i = 0; i < x.n_elem; i++) {
                BoysFunction boys(x(i), n, &intermediate);
                results << x(i) << " " << boys.result(testLevel) << endl;
            }
            results.close();
        }
    }

    TEST(BoysSingleResult) {
        BoysFunctionIntermediate intermediate(5, 1000, 0, 30, 5e3);
        BoysFunction boys(3.13137, 5, &intermediate);
        cout << boys.result(5) << endl;
    }

    TEST(BoysIntermediateTest) {
        BoysFunctionIntermediate boysIntermediate(20, 1000);
        boysIntermediate.updateResults();
        cout << "Boys, boys, boys = " << boysIntermediate.result(6, 10) << endl;
    }
}

SUITE(Development) {
    TEST(GaussianTypeOverlapIntegralTest) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};
        GaussianTypeOverlapIntegral integrator(posA, posB, 0.2, 0.3, 3);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,0), 0.119172363580852, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,1), 0.276479883507577, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,2), 0.760605693318432, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,0), 0.0429020508891068, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,1), 0.0995327580627279, 0.00001);

    }


    TEST(GaussianTypeKineticIntegralTest) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};
        GaussianTypeKineticIntegral integrator(posA, posB, 0.2, 0.3, 3);

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
    }

    TEST(GaussianTypeColoumbAttractionIntegralTest) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};
        rowvec posC = {2.3, 0.9, 3.2};
        double a = 0.2;
        double b = 0.3;
        GaussianTypeColoumbAttractionIntegral integrator(posA, posB, posC, a, b, 3);
        CHECK_CLOSE(2.788948987251e-02, integrator.coloumbAttractionIntegral(0,0,0,0,0,0), 1e-4);
        CHECK_CLOSE(6.971203468743e-02, integrator.coloumbAttractionIntegral(0,0,0,0,0,1), 1e-4);
        CHECK_CLOSE(2.024071525839e-01, integrator.coloumbAttractionIntegral(0,0,0,0,0,2), 1e-4);
        CHECK_CLOSE(8.727033700014e-03, integrator.coloumbAttractionIntegral(0,0,0,0,1,0), 1e-4);
        CHECK_CLOSE(2.134361291529e-02, integrator.coloumbAttractionIntegral(0,0,0,0,1,1), 1e-4);
        CHECK_CLOSE(2.921666495443e-02, integrator.coloumbAttractionIntegral(0,0,0,0,2,0), 1e-4);
        CHECK_CLOSE(3.185957751329e-02, integrator.coloumbAttractionIntegral(0,0,0,1,0,0), 1e-4);
        CHECK_CLOSE(8.105746642202e-02, integrator.coloumbAttractionIntegral(0,0,0,1,0,1), 1e-4);
        CHECK_CLOSE(9.596523510045e-03, integrator.coloumbAttractionIntegral(0,0,0,1,1,0), 1e-4);
        CHECK_CLOSE(6.388444040338e-02, integrator.coloumbAttractionIntegral(0,0,0,2,0,0), 1e-4);
        CHECK_CLOSE(-9.204700718547e-02, integrator.coloumbAttractionIntegral(0,0,1,0,0,0), 1e-4);
        CHECK_CLOSE(-2.019226499202e-01, integrator.coloumbAttractionIntegral(0,0,1,0,0,1), 1e-4);
        CHECK_CLOSE(-5.276399683274e-01, integrator.coloumbAttractionIntegral(0,0,1,0,0,2), 1e-4);
        CHECK_CLOSE(-2.927318225377e-02, integrator.coloumbAttractionIntegral(0,0,1,0,1,0), 1e-4);
        CHECK_CLOSE(-6.299787002318e-02, integrator.coloumbAttractionIntegral(0,0,1,0,1,1), 1e-4);
        CHECK_CLOSE(-9.718105595370e-02, integrator.coloumbAttractionIntegral(0,0,1,0,2,0), 1e-4);
        CHECK_CLOSE(-1.037280861539e-01, integrator.coloumbAttractionIntegral(0,0,1,1,0,0), 1e-4);
        CHECK_CLOSE(-2.312309453843e-01, integrator.coloumbAttractionIntegral(0,0,1,1,0,1), 1e-4);
        CHECK_CLOSE(-3.202910466576e-02, integrator.coloumbAttractionIntegral(0,0,1,1,1,0), 1e-4);
        CHECK_CLOSE(-2.073449397904e-01, integrator.coloumbAttractionIntegral(0,0,1,2,0,0), 1e-4);
        CHECK_CLOSE(3.319499900436e-01, integrator.coloumbAttractionIntegral(0,0,2,0,0,0), 1e-4);
        CHECK_CLOSE(6.435114042344e-01, integrator.coloumbAttractionIntegral(0,0,2,0,0,1), 1e-4);
        CHECK_CLOSE(1.536931448007e+00, integrator.coloumbAttractionIntegral(0,0,2,0,0,2), 1e-4);
        CHECK_CLOSE(1.067865861209e-01, integrator.coloumbAttractionIntegral(0,0,2,0,1,0), 1e-4);
        CHECK_CLOSE(2.033153544029e-01, integrator.coloumbAttractionIntegral(0,0,2,0,1,1), 1e-4);
        CHECK_CLOSE(3.524622701603e-01, integrator.coloumbAttractionIntegral(0,0,2,0,2,0), 1e-4);
        CHECK_CLOSE(3.703919381580e-01, integrator.coloumbAttractionIntegral(0,0,2,1,0,0), 1e-4);
        CHECK_CLOSE(7.292169308884e-01, integrator.coloumbAttractionIntegral(0,0,2,1,0,1), 1e-4);
        CHECK_CLOSE(1.162963233448e-01, integrator.coloumbAttractionIntegral(0,0,2,1,1,0), 1e-4);
        CHECK_CLOSE(7.390872806284e-01, integrator.coloumbAttractionIntegral(0,0,2,2,0,0), 1e-4);
        CHECK_CLOSE(-1.637350724302e-02, integrator.coloumbAttractionIntegral(0,1,0,0,0,0), 1e-4);
        CHECK_CLOSE(-4.139721853567e-02, integrator.coloumbAttractionIntegral(0,1,0,0,0,1), 1e-4);
        CHECK_CLOSE(-1.213713540367e-01, integrator.coloumbAttractionIntegral(0,1,0,0,0,2), 1e-4);
        CHECK_CLOSE(2.136233458775e-02, integrator.coloumbAttractionIntegral(0,1,0,0,1,0), 1e-4);
        CHECK_CLOSE(5.306634838230e-02, integrator.coloumbAttractionIntegral(0,1,0,0,1,1), 1e-4);
        CHECK_CLOSE(-1.697963144320e-04, integrator.coloumbAttractionIntegral(0,1,0,0,2,0), 1e-4);
        CHECK_CLOSE(-1.907709578263e-02, integrator.coloumbAttractionIntegral(0,1,0,1,0,0), 1e-4);
        CHECK_CLOSE(-4.932098923684e-02, integrator.coloumbAttractionIntegral(0,1,0,1,0,1), 1e-4);
        CHECK_CLOSE(2.414126847830e-02, integrator.coloumbAttractionIntegral(0,1,0,1,1,0), 1e-4);
        CHECK_CLOSE(-3.842342257899e-02, integrator.coloumbAttractionIntegral(0,1,0,2,0,0), 1e-4);
        CHECK_CLOSE(5.356912442721e-02, integrator.coloumbAttractionIntegral(0,1,1,0,0,0), 1e-4);
        CHECK_CLOSE(1.187325132374e-01, integrator.coloumbAttractionIntegral(0,1,1,0,0,1), 1e-4);
        CHECK_CLOSE(3.128036711345e-01, integrator.coloumbAttractionIntegral(0,1,1,0,0,2), 1e-4);
        CHECK_CLOSE(-7.083519267815e-02, integrator.coloumbAttractionIntegral(0,1,1,0,1,0), 1e-4);
        CHECK_CLOSE(-1.544897767601e-01, integrator.coloumbAttractionIntegral(0,1,1,0,1,1), 1e-4);
        CHECK_CLOSE(-3.894393296797e-04, integrator.coloumbAttractionIntegral(0,1,1,0,2,0), 1e-4);
        CHECK_CLOSE(6.132616876012e-02, integrator.coloumbAttractionIntegral(0,1,1,1,0,0), 1e-4);
        CHECK_CLOSE(1.386353605834e-01, integrator.coloumbAttractionIntegral(0,1,1,1,0,1), 1e-4);
        CHECK_CLOSE(-7.911527548440e-02, integrator.coloumbAttractionIntegral(0,1,1,1,1,0), 1e-4);
        CHECK_CLOSE(1.229240947242e-01, integrator.coloumbAttractionIntegral(0,1,1,2,0,0), 1e-4);
        CHECK_CLOSE(3.609849112824e-02, integrator.coloumbAttractionIntegral(0,2,0,0,0,0), 1e-4);
        CHECK_CLOSE(9.032384498820e-02, integrator.coloumbAttractionIntegral(0,2,0,0,0,1), 1e-4);
        CHECK_CLOSE(2.625292648498e-01, integrator.coloumbAttractionIntegral(0,2,0,0,0,2), 1e-4);
        CHECK_CLOSE(-1.939589748931e-02, integrator.coloumbAttractionIntegral(0,2,0,0,1,0), 1e-4);
        CHECK_CLOSE(-4.913397190183e-02, integrator.coloumbAttractionIntegral(0,2,0,0,1,1), 1e-4);
        CHECK_CLOSE(6.878262296370e-02, integrator.coloumbAttractionIntegral(0,2,0,0,2,0), 1e-4);
        CHECK_CLOSE(4.131065513841e-02, integrator.coloumbAttractionIntegral(0,2,0,1,0,0), 1e-4);
        CHECK_CLOSE(1.052929737663e-01, integrator.coloumbAttractionIntegral(0,2,0,1,0,1), 1e-4);
        CHECK_CLOSE(-2.267402937768e-02, integrator.coloumbAttractionIntegral(0,2,0,1,1,0), 1e-4);
        CHECK_CLOSE(8.289710831960e-02, integrator.coloumbAttractionIntegral(0,2,0,2,0,0), 1e-4);
        CHECK_CLOSE(-3.786414780578e-02, integrator.coloumbAttractionIntegral(1,0,0,0,0,0), 1e-4);
        CHECK_CLOSE(-9.322262110550e-02, integrator.coloumbAttractionIntegral(1,0,0,0,0,1), 1e-4);
        CHECK_CLOSE(-2.671155215998e-01, integrator.coloumbAttractionIntegral(1,0,0,0,0,2), 1e-4);
        CHECK_CLOSE(-1.222106053447e-02, integrator.coloumbAttractionIntegral(1,0,0,0,1,0), 1e-4);
        CHECK_CLOSE(-2.972830178046e-02, integrator.coloumbAttractionIntegral(1,0,0,0,1,1), 1e-4);
        CHECK_CLOSE(-4.026352276293e-02, integrator.coloumbAttractionIntegral(1,0,0,0,2,0), 1e-4);
        CHECK_CLOSE(-1.576450345257e-02, integrator.coloumbAttractionIntegral(1,0,0,1,0,0), 1e-4);
        CHECK_CLOSE(-3.945885129414e-02, integrator.coloumbAttractionIntegral(1,0,0,1,0,1), 1e-4);
        CHECK_CLOSE(-4.918734877201e-03, integrator.coloumbAttractionIntegral(1,0,0,1,1,0), 1e-4);
        CHECK_CLOSE(-2.459437143524e-02, integrator.coloumbAttractionIntegral(1,0,0,2,0,0), 1e-4);
        CHECK_CLOSE(1.263894353489e-01, integrator.coloumbAttractionIntegral(1,0,1,0,0,0), 1e-4);
        CHECK_CLOSE(2.735756798558e-01, integrator.coloumbAttractionIntegral(1,0,1,0,0,1), 1e-4);
        CHECK_CLOSE(7.071773603054e-01, integrator.coloumbAttractionIntegral(1,0,1,0,0,2), 1e-4);
        CHECK_CLOSE(4.115384967585e-02, integrator.coloumbAttractionIntegral(1,0,1,0,1,0), 1e-4);
        CHECK_CLOSE(8.802219023191e-02, integrator.coloumbAttractionIntegral(1,0,1,0,1,1), 1e-4);
        CHECK_CLOSE(1.350111738398e-01, integrator.coloumbAttractionIntegral(1,0,1,0,2,0), 1e-4);
        CHECK_CLOSE(5.197526800952e-02, integrator.coloumbAttractionIntegral(1,0,1,1,0,0), 1e-4);
        CHECK_CLOSE(1.145639876363e-01, integrator.coloumbAttractionIntegral(1,0,1,1,0,1), 1e-4);
        CHECK_CLOSE(1.638641395940e-02, integrator.coloumbAttractionIntegral(1,0,1,1,1,0), 1e-4);
        CHECK_CLOSE(8.278875254192e-02, integrator.coloumbAttractionIntegral(1,0,1,2,0,0), 1e-4);
        CHECK_CLOSE(2.185667243163e-02, integrator.coloumbAttractionIntegral(1,1,0,0,0,0), 1e-4);
        CHECK_CLOSE(5.417205698627e-02, integrator.coloumbAttractionIntegral(1,1,0,0,0,1), 1e-4);
        CHECK_CLOSE(1.560020091608e-01, integrator.coloumbAttractionIntegral(1,1,0,0,0,2), 1e-4);
        CHECK_CLOSE(-2.926456829930e-02, integrator.coloumbAttractionIntegral(1,1,0,0,1,0), 1e-4);
        CHECK_CLOSE(-7.176178735649e-02, integrator.coloumbAttractionIntegral(1,1,0,0,1,1), 1e-4);
        CHECK_CLOSE(-5.223967979758e-04, integrator.coloumbAttractionIntegral(1,1,0,0,2,0), 1e-4);
        CHECK_CLOSE(9.269318129877e-03, integrator.coloumbAttractionIntegral(1,1,0,1,0,0), 1e-4);
        CHECK_CLOSE(2.337071697343e-02, integrator.coloumbAttractionIntegral(1,1,0,1,0,1), 1e-4);
        CHECK_CLOSE(-1.203714316117e-02, integrator.coloumbAttractionIntegral(1,1,0,1,1,0), 1e-4);
        CHECK_CLOSE(1.401501778682e-02, integrator.coloumbAttractionIntegral(1,1,0,2,0,0), 1e-4);
        CHECK_CLOSE(7.889586550718e-02, integrator.coloumbAttractionIntegral(2,0,0,0,0,0), 1e-4);
        CHECK_CLOSE(1.935977010010e-01, integrator.coloumbAttractionIntegral(2,0,0,0,0,1), 1e-4);
        CHECK_CLOSE(5.534914541236e-01, integrator.coloumbAttractionIntegral(2,0,0,0,0,2), 1e-4);
        CHECK_CLOSE(2.563391673303e-02, integrator.coloumbAttractionIntegral(2,0,0,0,1,0), 1e-4);
        CHECK_CLOSE(6.217850538435e-02, integrator.coloumbAttractionIntegral(2,0,0,0,1,1), 1e-4);
        CHECK_CLOSE(8.419480232293e-02, integrator.coloumbAttractionIntegral(2,0,0,0,2,0), 1e-4);
        CHECK_CLOSE(1.481688684288e-02, integrator.coloumbAttractionIntegral(2,0,0,1,0,0), 1e-4);
        CHECK_CLOSE(3.878852644576e-02, integrator.coloumbAttractionIntegral(2,0,0,1,0,1), 1e-4);
        CHECK_CLOSE(4.176920693786e-03, integrator.coloumbAttractionIntegral(2,0,0,1,1,0), 1e-4);
        CHECK_CLOSE(6.422210627967e-02, integrator.coloumbAttractionIntegral(2,0,0,2,0,0), 1e-4);
    }

    TEST(GaussianTypeElectronInteractionIntegralTest1) {
        cout << "Number one!" << endl;
        rowvec posA = {-0.5, 0, 0};
        rowvec posB = {-0.5, 0, 0};
        rowvec posC = {-0.5, 0, 0};
        rowvec posD = {-0.5, 0, 0};
        double a = 13.0077;
        double b = 13.0077;
        double c = 13.0077;
        double d = 13.0077;
        GaussianTypeElectronInteractionIntegral integrator(posA, posB, posC, posD, a, b, c, d, 3);
        cout << integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0) << endl;
    }

    TEST(GaussianTypeElectronInteractionIntegralTest2) {
        cout << "Number two!" << endl;
        rowvec posA = {0.5, 0, 0};
        rowvec posB = {-0.5, 0, 0};
        rowvec posC = {-0.5, 0, 0};
        rowvec posD = {0.5, 0, 0};
        double a = 13.0077;
        double b = 0.121949;
        double c = 0.444529;
        double d = 13.0077;
        GaussianTypeElectronInteractionIntegral integrator(posA, posB, posC, posD, a, b, c, d, 3);
        cout << integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0) << endl;
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

