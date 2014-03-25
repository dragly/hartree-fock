#include <unittest++/UnitTest++.h>

#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <hartreefocksolver.h>

SUITE(Systems) {
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
        solver.setConvergenceTreshold(1e-6);
        solver.setNIterationsMax(1e3);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
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
        mat C;
        HartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-6);
        solver.setNIterationsMax(1e3);
        for(int i = 0; i < 1000; i++) {
            solver.advance();
        }
        CHECK_CLOSE(-1.122933363617109, solver.energy(), 1e-6);
    }
    TEST(Neon321) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "atom_10_basis_3-21G.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-6);
        solver.setNIterationsMax(1e3);
//        for(int i = 0; i < 100; i++) {
//            solver.advance();
//        }
        solver.solve();
        CHECK_CLOSE(-127.8038139236938, solver.energy(), 1e-5);
//        cout << "Neon 3-21G: " << solver.energy() << endl;
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
        HartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-6);
        solver.setNIterationsMax(1e3);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        solver.overlapMatrix().save("S.mat", raw_ascii);
//        cout << "Oxygen energy: " << solver.energy() << endl;
        CHECK_CLOSE(-149.5111286894001, solver.energy(), 1e-5);
//        cout << "Oxygen 6-311G: " << solver.energy() << endl;
    }
}
