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
        HartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.5);
        cout << solver.coupledMatrix() << endl;
        solver.solve();
        CHECK_CLOSE(-75.90736859918989, solver.energy(), 1e-6);
    }

    //    TEST(OxygenSixAsterix) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gs.tm"));
    //        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-31Gs.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        mat C;
    //        UnrestrictedHartreeFockSolver solver(&system);
    //        solver.setConvergenceTreshold(1e-12);
    //        solver.setNIterationsMax(1e4);
    //        solver.setDensityMixFactor(0.0);
    //        solver.solve();
    //        cout << "Used: " << solver.iterationsUsed() << endl;
    //        cout << "Energy: " << solver.energy() << endl;
    ////            CHECK_CLOSE(-149.5295405478762, solver.energy(), 1e-5);
    //        CHECK_CLOSE(-149.5876095851103, solver.energy(), 1e-5);
    ////        cout << "Oxygen 6-311G: " << solver.energy() << endl;
    //    }

    //    TEST(CheapOxygen) {
    //        vector<GaussianCore> cores;
    //        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_cheap.tm"));
    //        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_cheap.tm"));
    //        GaussianSystem system;
    //        for(const GaussianCore &core : cores) {
    //            system.addCore(core);
    //        }
    //        mat C;
    //        HartreeFockSolver solver(&system);
    //        solver.setConvergenceTreshold(1e-12);
    //        solver.setNIterationsMax(1e4);
    //        solver.setDensityMixFactor(0.0);
    //        solver.advance();
    ////        solver.solve();
    //        cout << "Used " << solver.iterationsUsed() << endl;
    //        cout << "Oxygen cheap: " << solver.energy() << endl;
    //    }

    TEST(GTOelectronElectronIntegral)
    {
        GaussianElectronInteractionIntegral integrator(2);

        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};
        rowvec posC = {2.3,0.9,3.2};
        rowvec posD = {5.0,1.9,1.2};

        GaussianPrimitiveOrbital primitiveA(1.0, 2, 2, 2, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 2, 2, 2, 0.3);
        GaussianPrimitiveOrbital primitiveC(1.0, 2, 2, 2, 0.4);
        GaussianPrimitiveOrbital primitiveD(1.0, 2, 2, 2, 0.1);

        integrator.set(posA, posB, posC, posD, primitiveA, primitiveB, primitiveC, primitiveD);

        CHECK_CLOSE(1.624848e-01, integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0), 1e-5);
        CHECK_CLOSE(0.2667434785828074, integrator.electronInteractionIntegral(0,0,0, 1,0,0, 0,0,0, 0,0,1), 1e-7);
        CHECK_CLOSE(0.2681206720738772, integrator.electronInteractionIntegral(0,0,0, 1,0,0, 0,2,0, 0,0,1), 1e-7);
        CHECK_CLOSE(0.5266872995197744, integrator.electronInteractionIntegral(1,1,0, 2,0,0, 2,0,0, 2,0,0), 1e-7);
        CHECK_CLOSE(-0.1273045183436938, integrator.electronInteractionIntegral(1,1,0, 0,2,0, 0,2,0, 2,0,0), 1e-7);
    }

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
