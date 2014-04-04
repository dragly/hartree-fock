#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <solvers/unrestrictedhartreefocksolver.h>

#include <iostream>

#include <unittest++/UnitTest++.h>

using std::cout;
using std::endl;

SUITE(Unrestricted) {
    TEST(HydrogenMolecule) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_3-21G.tm"));
        cores.push_back(GaussianCore({1.4,0.0,0}, "atom_1_basis_3-21G.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(-1.122933363617109, solver.energy(), 1e-6);
    }
}
