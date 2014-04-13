#include "electronsystems/gaussian/gaussiancore.h"
#include "electronsystems/gaussian/gaussiansystem.h"
#include "solvers/unrestrictedhartreefocksolver.h"
#include "basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h"
#include "basisfunctions/gaussian/integrals/gaussianoverlapintegral.h"
#include "math/boysfunction.h"

#include <iostream>

#include <unittest++/UnitTest++.h>

using std::cout;
using std::endl;

SUITE(Development) {
//    TEST(OxygenSix) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_4-31G.tm"));
//        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_4-31G.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        HartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e4);
//        solver.setDensityMixFactor(0.5);
//        solver.solve();
//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
////        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
//    TEST(Water) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_4-31G.tm"));
//        cores.push_back(GaussianCore({-1.43,1.108,0}, "atom_1_basis_4-31G.tm"));
//        cores.push_back(GaussianCore({1.43,1.108,0}, "atom_1_basis_4-31G.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        HartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e3);
//        solver.setDensityMixFactor(0.5);
//        solver.solve();
//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
//        CHECK_CLOSE(-75.90736859918989, solver.energy(), 1e-6);
//    }
    TEST(OxygenSixAsterisk) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gs.tm"));
        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-31Gs.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-8);
        solver.setNIterationsMax(1e4);
        solver.setDensityMixFactor(0.5);
        solver.solve();
        CHECK_CLOSE(-149.5876095851103, solver.energy(), 1e-5);
    }

//        TEST(CheapOxygen) {
//            vector<GaussianCore> cores;
//            cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_cheap.tm"));
//            cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_cheap.tm"));
//            GaussianSystem system;
//            for(const GaussianCore &core : cores) {
//                system.addCore(core);
//            }
//            mat C;
//            HartreeFockSolver solver(&system);
//            solver.setConvergenceTreshold(1e-12);
//            solver.setNIterationsMax(1e4);
//            solver.setDensityMixFactor(0.0);
//            solver.advance();
//    //        solver.solve();
//            cout << "Used " << solver.iterationsUsed() << endl;
//            cout << "Oxygen cheap: " << solver.energy() << endl;
//        }

    //    TEST(Variations) {
    //        for(int i = 0; i < 100; i++) {
    //            vector<GaussianCore> cores;
    //            cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_4-31G.tm"));
    //            cores.push_back(GaussianCore({0.5 + i * 0.1,0,0}, "atom_1_basis_4-31G.tm"));
    //            GaussianSystem system;
    //            for(const GaussianCore &core : cores) {
    //                system.addCore(core);
    //            }
    //            UnrestrictedHartreeFockSolver solver(&system);
    //            solver.setConvergenceTreshold(1e-9);
    //            solver.setNIterationsMax(1e3);
    //            solver.setDensityMixFactor(0.5);
    //            solver.solve();
    //            cout << solver.energy() << endl;
    //        }
    //    }
    //    TEST(Water) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_3-21G.tm"));
    //        cores.push_back(GaussianCore({-1.43,1.108,0}, "atom_1_basis_3-21G.tm"));
    //        cores.push_back(GaussianCore({1.43,1.108,0}, "atom_1_basis_3-21G.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        mat C;
    //        HartreeFockSolver solver(&system);
    //        solver.setConvergenceTreshold(1e-12);
    //        solver.setNIterationsMax(1e3);
    //        solver.setDensityMixFactor(0.0);
    //        solver.solve();
    //        CHECK_CLOSE(solver.energy(), -75.90736859918989, 1e-6);
    //    }
    //    TEST(HydrogenMolecule) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_6-31Gds.tm"));
    //        cores.push_back(GaussianCore({1.4,0.0,0}, "atom_1_basis_6-31Gds.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        HartreeFockSolver solver(&system);
    //        solver.setConvergenceTreshold(1e-12);
    //        solver.setNIterationsMax(1e3);
    //        solver.setDensityMixFactor(0.5);
    //        solver.solve();
    //        CHECK_CLOSE(-1.122933363617109, solver.energy(), 1e-6);
    //    }
    //    TEST(WaterAsterix) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gs.tm"));
    //        cores.push_back(GaussianCore({-1.43,1.108,0}, "atom_1_basis_3-21G.tm"));
    //        cores.push_back(GaussianCore({1.43,1.108,0}, "atom_1_basis_3-21G.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        mat C;
    //        HartreeFockSolver solver(&system);
    //        solver.setConvergenceTreshold(1e-12);
    //        solver.setNIterationsMax(1e3);
    //        solver.setDensityMixFactor(0.0);
    //        solver.solve();
    //        CHECK_CLOSE(-75.90736859918989, solver.energy() , 1e-6);
    //    }
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
