#include <unittest++/UnitTest++.h>

#include <electronsystems/gaussian/gaussianoxygen321g.h>
#include <hartreefocksolver.h>

SUITE(GaussianOxygen) {
    TEST(GaussianOxygen) {
        GaussianOxygen431G system;
        cout << system.overlapIntegral(0,0) << endl;
        cout << system.uncoupledIntegral(0,0) << endl;
        cout << system.coupledIntegral(0,0,0,0) << endl;
    }
    TEST(GaussianOxygenHF) {
        GaussianOxygen431G system;
        cout << "Setting up solver" << endl;
        HartreeFockSolver solver(&system);
        cout << "Starting loop" << endl;
        for(int i = 0; i < 100; i++) {
//            cout << i << " " << solver.energy() << endl;
            solver.advance();
        }
        cout << "Result energy: " << solver.energy() << endl;
    }
}
