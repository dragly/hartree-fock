#include "electronsystems/gaussian/gaussiancore.h"
#include "electronsystems/gaussian/gaussiansystem.h"
#include "solvers/unrestrictedhartreefocksolver.h"
#include "solvers/restrictedhartreefocksolver.h"
#include "basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h"
#include "basisfunctions/gaussian/integrals/gaussianoverlapintegral.h"
#include "math/boysfunction.h"

#include <iostream>
#include <iomanip>

#include <unittest++/UnitTest++.h>

using std::cout;
using std::endl;

SUITE(Development) {
    TEST(Dummy) {
    }
    TEST(LinearWater) {
        vector<GaussianCore> cores;
        string basis = "6-311++Gdsds";
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_" + basis + ".tm"));
        cores.push_back(GaussianCore({1.8, 0, 0}, "atom_1_basis_" + basis + ".tm"));
        cores.push_back(GaussianCore({-0.44916, 1.743056, 0.0}, "atom_1_basis_" + basis + ".tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        UnrestrictedHartreeFockSolver solver(&system);
        mat coefficientsUp = randn(system.nBasisFunctions(), system.nParticlesUp());
        mat coefficientsDown = randn(system.nBasisFunctions(), system.nParticlesDown());
        solver.setInitialCoefficientMatrices(coefficientsUp, coefficientsDown);
        solver.setConvergenceTreshold(1e-9);
        solver.setNIterationsMax(1e4);
        solver.setDensityMixFactor(0.95);
        solver.solve();

        vector<GaussianCore> cores2;
        cores2.push_back(GaussianCore({0,0,0}, "atom_8_basis_" + basis + ".tm"));
        cores2.push_back(GaussianCore({8,0,0}, "atom_1_basis_" + basis + ".tm"));
        cores2.push_back(GaussianCore({-8,0,0}, "atom_1_basis_" + basis + ".tm"));
        GaussianSystem system2;
        for(const GaussianCore &core : cores2) {
            system2.addCore(core);
        }
        UnrestrictedHartreeFockSolver solver2(&system);
        solver2.setInitialCoefficientMatrices(solver.coeffcientMatrixUp(), solver.coeffcientMatrixDown());
        solver2.setConvergenceTreshold(1e-9);
        solver2.setNIterationsMax(1e4);
        solver2.setDensityMixFactor(0.95);
        solver2.solve();
        cout << std::setprecision(20);
        cout << "Energy: " << solver2.energy() << endl;
        cout << "Iterations: " << solver2.iterationsUsed() << endl;
        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
    }
//    TEST(LinearHydrogen) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_1_basis_6-31Gdsds.tm"));
//        cores.push_back(GaussianCore({12,0,0}, "atom_1_basis_6-31Gdsds.tm"));
//        cores.push_back(GaussianCore({24,0,0}, "atom_1_basis_6-31Gdsds.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        UnrestrictedHartreeFockSolver solver(&system);
//        mat coefficientsUp = randn(system.nBasisFunctions(), system.nParticlesUp());
//        mat coefficientsDown = randn(system.nBasisFunctions(), system.nParticlesDown());
//        solver.setInitialCoefficientMatrices(coefficientsUp, coefficientsDown);
//        solver.setConvergenceTreshold(1e-9);
//        solver.setNIterationsMax(1e3);
//        solver.setDensityMixFactor(0.95);
//        solver.solve();
//        cout << solver.overlapMatrix() << endl;
//        cout << std::setprecision(20);
//        cout << "Energy: " << solver.energy() << endl;
//        cout << "Iterations: " << solver.iterationsUsed() << endl;
//        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
//    TEST(OxygenSix) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_6_basis_6-311++Gdsds.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        UnrestrictedHartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e4);
//        solver.setDensityMixFactor(0.95);
//        solver.solve();
//        cout << std::setprecision(20);
//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
//        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
//    TEST(Water431G) {
//        vector<GaussianCore> cores;
////        1) O: [0, 0, 0]    H: [1.7961, 0, 0]    H: [-0.6506, 1.6742, 0]
////        2) O: [0, 0, 0]    H: [1.7970, 0, 0]    H: [-0.6509, 1.6750, 0]
//        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_4-31G.tm"));
//        cores.push_back(GaussianCore({1.7970, 0, 0}, "atom_1_basis_4-31G.tm"));
//        cores.push_back(GaussianCore({-0.448, 1.740, 0.0}, "atom_1_basis_4-31G.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        UnrestrictedHartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e4);
//        solver.setDensityMixFactor(0.95);
//        solver.solve();
//        cout << std::setprecision(20);
//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
//        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
//    TEST(WaterSixPlusPlus) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-311++Gdsds.tm"));
//        cores.push_back(GaussianCore({1.797, 0.0, 0.0}, "atom_1_basis_6-311++Gdsds.tm"));
//        cores.push_back(GaussianCore({-0.448, 1.740, 0.0}, "atom_1_basis_6-311++Gdsds.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        UnrestrictedHartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e4);
//        solver.setDensityMixFactor(0.95);
//        solver.solve();
//        cout << std::setprecision(20);

//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
//        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
//    TEST(OxygenSix) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gds.tm"));
//        cores.push_back(GaussianCore({2.281,0,0}, "atom_8_basis_6-31Gds.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        system.setNParticlesDown(7);
//        mat C;
//        UnrestrictedHartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e4);
//        solver.setDensityMixFactor(0.95);
//        solver.solve();
//        cout << std::setprecision(20);
//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
//        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
//    TEST(BoronHydrogen) {
//        vector<GaussianCore> cores;
//        cores.push_back(GaussianCore({0,0,0}, "atom_5_basis_STO-3G.tm"));
//        cores.push_back(GaussianCore({1.2,0,0}, "atom_1_basis_STO-3G.tm"));
//        GaussianSystem system;
//        for(const GaussianCore &core : cores) {
//            system.addCore(core);
//        }
//        mat C;
//        system.setNParticlesDown(4);
//        UnrestrictedHartreeFockSolver solver(&system);
//        solver.setConvergenceTreshold(1e-12);
//        solver.setNIterationsMax(1e4);
//        solver.setDensityMixFactor(0.95);
//        solver.solve();
//        cout << std::setprecision(20);
//        cout << solver.energy() << endl;
//        cout << solver.iterationsUsed() << endl;
//        //        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//    }
    //    TEST(OxygenSix) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_4-31G.tm"));
    //        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_4-31G.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        mat C;
    //        UnrestrictedHartreeFockSolver solver(&system);
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

    //    TEST(OxygenSixAsterisk) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gs.tm"));
    //        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-31Gs.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        mat C;
    //        UnrestrictedHartreeFockSolver solver(&system);
    //        solver.setConvergenceTreshold(1e-8);
    //        solver.setNIterationsMax(1e4);
    //        solver.setDensityMixFactor(0.5);
    //        solver.solve();
    //        cout << solver.iterationsUsed() << endl;
    //        CHECK_CLOSE(-149.5876095851103, solver.energy(), 1e-5);
    //    }

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

}
