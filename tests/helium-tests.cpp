#include <iostream>
#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <src/hartreesolver.h>
#include <src/hartreefocksolver.h>
#include <src/electronsystems/helium/heliumhartree.h>
#include <src/electronsystems/hydrogen/hydrogenmolecule.h>
#include <src/electronsystems/hydrogen/multihydrogen.h>

#include <fstream>

using namespace std;

SUITE(Helium) {

    TEST(MatrixElement) {
        Helium basisFunction;
        int p = 0;
        int q = 1;
        int r = 2;
        int s = 3;

        double value = basisFunction.electronInteractionIntegral(p,r,q,s);

        CHECK_CLOSE(0.075820496873881579, value, 1e-9);
    }

    TEST(KineticIntegral) {
        Helium basisFunction;
        int p = 0;
        int q = 3;

        double value = basisFunction.kineticIntegral(p,q);

        CHECK_CLOSE(0.0204654960291562, value, 1e-9);
    }

    TEST(NuclearAttractionIntegral) {
        Helium basisFunction;
        int p = 1;
        int q = 2;

        double value = basisFunction.nuclearAttractionIntegral(p,q);

        CHECK_CLOSE(-1.78867607774792, value, 1e-9);
    }

    TEST(OverlapIntegral) {
        Helium basisFunction;
        int p = 3;
        int q = 2;

        double value = basisFunction.overlapIntegral(p, q);

        CHECK_CLOSE(0.01891203832678, value, 1e-9);
    }

    TEST(HeliumAdvanceMany) {
        Helium basisFunction;
        HartreeSolver solver(&basisFunction);
        solver.reset();
        for(int i = 0; i < 500; i++) {
            solver.advance();
        }
        double value = solver.energy();
        CHECK_CLOSE(-2.855160382370257377, value, 1e-7);
    }

}
