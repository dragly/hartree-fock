#include <unittest++/UnitTest++.h>

#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include "solvers/restrictedhartreefocksolver.h"
#include <solvers/unrestrictedhartreefocksolver.h>
#include <iomanip>

using namespace std;

SUITE(Systems) {
    TEST(Dummy) {

    }

    TEST(Water) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_4-31G.tm"));
        cores.push_back(GaussianCore({-1.43,1.108,0}, "atom_1_basis_4-31G.tm"));
        cores.push_back(GaussianCore({1.43,1.108,0}, "atom_1_basis_4-31G.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(solver.energy(), -75.90736859918989, 1e-6);
    }
    TEST(HydrogenMolecule) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_3-21G.tm"));
        cores.push_back(GaussianCore({1.4,0.0,0}, "atom_1_basis_3-21G.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(-1.122933363617109, solver.energy(), 1e-6);
    }
    TEST(HydrogenMolecule631Gdsds) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_6-31Gdsds.tm"));
        cores.push_back(GaussianCore({1.4,0.0,0}, "atom_1_basis_6-31Gdsds.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(-1.13128434930047, solver.energy(), 1e-6);
    }
    TEST(Neon321) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_10_basis_3-21G.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(-127.8038245281864, solver.energy(), 1e-5);
    }
    TEST(OxygenSix) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-311G.tm"));
        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-311G.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e4);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
    }
    TEST(OxygenSixAsterisk) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gds.tm"));
        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-31Gds.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-8);
        solver.setNIterationsMax(1e4);
        solver.setDensityMixFactor(0.5);
        solver.setDiisEnabled(false);
        solver.solve();
        CHECK_CLOSE(-149.5876095851103, solver.energy(), 1e-5);
    }
}
