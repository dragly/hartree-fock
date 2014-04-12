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
        solver.solve();
        CHECK_CLOSE(solver.energy(), -75.90736859918989, 1e-6);
    }

    TEST(OxygenSixAsterix) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gs.tm"));
        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-31Gs.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e4);
        solver.setDensityMixFactor(0.0);
        solver.solve();
        cout << "Used: " << solver.iterationsUsed() << endl;
        cout << "Energy: " << solver.energy() << endl;
//            CHECK_CLOSE(-149.5295405478762, solver.energy(), 1e-5);
        CHECK_CLOSE(-149.5876095851103, solver.energy(), 1e-5);
//        cout << "Oxygen 6-311G: " << solver.energy() << endl;
    }
    //    TEST(OverlapNew) {
    //        rowvec posA = {1.2,2.3,3.4};
    //        rowvec posB = {-1.3,1.4,-2.4};
    //        GaussianOverlapIntegral integrator(posA, posB, 0.2, 0.3, 3);
    //        CHECK_CLOSE(2.979309089521e-01, integrator.overlapIntegral(2,0,0,2,0,0), 1e-7);
    //        CHECK_CLOSE(1.072551272228e-02, integrator.overlapIntegral(2,0,0,1,1,0), 1e-7);
    //        CHECK_CLOSE(6.911997087690e-02, integrator.overlapIntegral(2,0,0,1,0,1), 1e-7);
    //    }
//    TEST(KineticNew) {
//        rowvec posA = {1.2,2.3,3.4};
//        rowvec posB = {-1.3,1.4,-2.4};
//        GaussianKineticIntegral integrator(posA, posB, 0.2, 0.3, 3);
//        CHECK_CLOSE(-3.468392469657e-01, integrator.kineticIntegral(2,0,0,2,0,0), 1e-7);
//        CHECK_CLOSE(-3.562586305833e-03, integrator.kineticIntegral(2,0,0,1,1,0), 1e-7);
//        CHECK_CLOSE(-2.295888952647e-02, integrator.kineticIntegral(2,0,0,1,0,1), 1e-7);
//        CHECK_CLOSE(2.514020826620e-01, integrator.kineticIntegral(0,0,1,0,0,1), 1e-7);
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

    //    TEST(GTOelectronElectronIntegral)
    //    {
    //        GaussianElectronInteractionIntegral integrator(2);

    //        rowvec posA = {1.2,2.3,3.4};
    //        rowvec posB = {-1.3,1.4,-2.4};
    //        rowvec posC = {2.3,0.9,3.2};
    //        rowvec posD = {5.0,1.9,1.2};

    //        integrator.set(posA, posB, posC, posD, 0.2, 0.3, 0.4, 0.1);

    //        CHECK_CLOSE(1.624848e-01, integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0), 1e-5);

    ////        primitiveA.setPowers({0,0,0});
    ////        primitiveB.setPowers({1,0,0});
    ////        primitiveC.setPowers({0,0,0});
    ////        primitiveD.setPowers({0,0,1});
    //        cout << integrator.electronInteractionIntegral(0,0,0, 1,0,0, 0,0,0, 0,0,1) << endl;

    ////        primitiveA.setPowers({0,0,0});
    ////        primitiveB.setPowers({1,0,0});
    ////        primitiveC.setPowers({0,2,0});
    ////        primitiveD.setPowers({0,0,1});
    ////        integrator.setPrimitiveA(primitiveA);
    ////        integrator.setPrimitiveB(primitiveB);
    ////        integrator.setPrimitiveC(primitiveC);
    ////        integrator.setPrimitiveD(primitiveD);
    //        cout << integrator.electronInteractionIntegral(0,0,0, 1,0,0, 0,2,0, 0,0,1) << endl;
    //        cout << integrator.electronInteractionIntegral(1,1,0, 2,0,0, 2,0,0, 2,0,0) << endl;
    //        cout << integrator.electronInteractionIntegral(1,1,0, 0,2,0, 0,2,0, 2,0,0) << endl;


    ////        primitiveA.setPowers({0,0,1});
    ////        primitiveB.setPowers({1,0,0});
    ////        primitiveC.setPowers({0,2,0});
    ////        primitiveD.setPowers({0,0,1});
    ////        integrator.setPrimitiveA(primitiveA);
    ////        integrator.setPrimitiveB(primitiveB);
    ////        integrator.setPrimitiveC(primitiveC);
    ////        integrator.setPrimitiveD(primitiveD);
    ////        CHECK_CLOSE(-1.109942e+02, integrator.electronRepulsionIntegral(),
    ////    1e-5);



    //    }

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
