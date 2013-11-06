#include <unittest++/UnitTest++.h>

#include <hartreesolver.h>
#include <hartreefocksolver.h>

#include <electronsystems/helium/heliumhartree.h>
#include <electronsystems/hydrogen/hydrogenmolecule.h>
#include <electronsystems/hydrogen/multihydrogen.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

SUITE(Hydrogen) {
    TEST(HydrogenOverlapIntegral)  {
        HydrogenMolecule basisFunction;
        double value = basisFunction.overlapIntegral(6,2);
        // Recursion test
        CHECK_CLOSE(4.296692816768918, value, 1e-9);
    }

    TEST(HydrogenElectronInteractionIntegral)  {
        HydrogenMolecule basisFunction;
        double value = basisFunction.electronInteractionIntegral(1,5,7,2);
        // Recursion test
        CHECK_CLOSE(0.9066377680574234, value, 1e-9);
    }

    TEST(HydrogenNuclearAttractionIntegral)  {
        HydrogenMolecule basisFunction;
        double value = basisFunction.nuclearAttractionIntegral(3,7);
        // Recursion test
        CHECK_CLOSE(-43.96134564895876, value, 1e-9);
    }

    TEST(HydrogenAdvanceMany) {
        HydrogenMolecule basisFunction(1.0);
        HartreeSolver solver(&basisFunction);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        // Recursion test
        CHECK_CLOSE(-2.0785476087914549481, solver.energy(), 1e-9);
    }

    TEST(HydrogenAdvanceManyHF) {
        HydrogenMolecule basisFunction(1.0);
        HartreeFockSolver solver(&basisFunction);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        // Recursion test
        CHECK_CLOSE(-1.0785476087914718235, solver.energy(), 1e-9);
    }

    TEST(MultiHydrogenAdvanceManyHF) {
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
        // Recursion test
        CHECK_CLOSE(-1.1186501652973053211, solver.energy(), 1e-9);
    }
}
