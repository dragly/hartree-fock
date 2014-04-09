#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <solvers/unrestrictedhartreefocksolver.h>

#include <iostream>

#include <unittest++/UnitTest++.h>

using std::cout;
using std::endl;

SUITE(Development) {
    TEST(Variations) {
        for(int i = 0; i < 100; i++) {
            vector<GaussianCore> cores;
            cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_4-31G.tm"));
            cores.push_back(GaussianCore({0.5 + i * 0.1,0,0}, "atom_1_basis_4-31G.tm"));
            GaussianSystem system;
            for(const GaussianCore &core : cores) {
                system.addCore(core);
            }
            UnrestrictedHartreeFockSolver solver(&system);
            solver.setConvergenceTreshold(1e-9);
            solver.setNIterationsMax(1e3);
            solver.setDensityMixFactor(0.5);
            solver.solve();
            cout << solver.energy() << endl;
        }
    }
//    TEST(HydrogenMolecule) {
//            vector<GaussianCore> cores;
////            cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_STO-3G.tm"));
////            cores.push_back(GaussianCore({1.4,0,0}, "atom_1_basis_STO-3G.tm"));
//            cores.push_back(GaussianCore({-0.7,0,0}, "atom_1_basis_3-21G.tm"));
//            cores.push_back(GaussianCore({ 0.7,0,0}, "atom_1_basis_3-21G.tm"));
//            GaussianSystem system;
//            for(const GaussianCore &core : cores) {
//                system.addCore(core);
//            }
//            HartreeFockSolver solver(&system);
////            solver.setConvergenceTreshold(1e-9);
////            solver.setNIterationsMax(1e3);
////            solver.advance();
////            cout << "----" << endl;
////            solver.advance();
//            solver.setDensityMixFactor(0.0);
//            solver.solve();
//            cout << "Number of iterations: " << solver.iterationsUsed() << endl;
//            cout << "Energy: " << solver.energy() << endl;
//        }
}
