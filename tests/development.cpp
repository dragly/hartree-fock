#include "electronsystems/gaussian/gaussiancore.h"
#include "electronsystems/gaussian/gaussiansystem.h"
#include "solvers/unrestrictedhartreefocksolver.h"
#include "basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h"

#include <iostream>

#include <unittest++/UnitTest++.h>

using std::cout;
using std::endl;

SUITE(Development) {
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
//        CHECK_CLOSE(3.297371e+01, integrator.electronInteractionIntegral(0,0,0, 1,0,0, 0,0,0, 0,0,1), 1e-5);

////        primitiveA.setPowers({0,0,0});
////        primitiveB.setPowers({1,0,0});
////        primitiveC.setPowers({0,2,0});
////        primitiveD.setPowers({0,0,1});
////        integrator.setPrimitiveA(primitiveA);
////        integrator.setPrimitiveB(primitiveB);
////        integrator.setPrimitiveC(primitiveC);
////        integrator.setPrimitiveD(primitiveD);
////        CHECK_CLOSE(3.333488e+01, integrator.electronRepulsionIntegral(),
////    1e-5);


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

    TEST(OxygenSixAsterix) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_8_basis_6-31Gs.tm"));
        cores.push_back(GaussianCore({2.282,0,0}, "atom_8_basis_6-31Gs.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setNIterationsMax(1e4);
        solver.setDensityMixFactor(0.0);
        solver.solve();
        CHECK_CLOSE(-149.5117583638509, solver.energy(), 1e-5);
//        cout << "Oxygen 6-311G: " << solver.energy() << endl;
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
