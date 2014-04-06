#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <solvers/unrestrictedhartreefocksolver.h>

#include <iostream>

#include <unittest++/UnitTest++.h>

using std::cout;
using std::endl;

SUITE(Development) {
    TEST(HydrogenMolecule2) {
        for(int i = 0; i < 1; i++) {
            vector<GaussianCore> cores;
            cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_4-31G.tm"));
            cores.push_back(GaussianCore({100.0,0,0}, "atom_1_basis_4-31G.tm"));
            GaussianSystem system;
            for(const GaussianCore &core : cores) {
                system.addCore(core);
            }
            UnrestrictedHartreeFockSolver solver(&system);
            solver.setConvergenceTreshold(1e-9);
            solver.setNIterationsMax(1e3);
            solver.setDensityMixFactor(0.5);
            solver.solve();
            cout << "H2 energy at high distance: " << solver.energy() << endl;
            cout << "Iterations: " << solver.iterationsUsed() << endl;
        }
    }
}
