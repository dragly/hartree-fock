#include <unittest++/UnitTest++.h>

#include <electronsystems/gaussian/gaussiancore.h>
#include <hartreefocksolver.h>

SUITE(Systems) {
    TEST(Water) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "oxygen431g.tm"));
        cores.push_back(GaussianCore({-1.43,1.108,0}, "hydrogen431g.tm"));
        cores.push_back(GaussianCore({1.43,1.108,0}, "hydrogen431g.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        CHECK_CLOSE(solver.energy(), -75.90736859918989, 1e-6);
    }
    TEST(NeonTest) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "neontest321g.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
//        for(int i = 0; i < 100; i++) {
//            solver.advance();
//        }
        solver.solve();
        cout << "Neon 3-21G s-limited: " << solver.energy() << endl;
    }
    TEST(Neon321) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "neon321g.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
//        for(int i = 0; i < 100; i++) {
//            solver.advance();
//        }
        solver.solve();
        cout << "Neon 3-21G: " << solver.energy() << endl;
    }
    TEST(OxygenSix) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "oxygen6311g.tm"));
        cores.push_back(GaussianCore({2.282,0,0}, "oxygen6311g.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        solver.overlapMatrix().save("S.mat", raw_ascii);
        cout << "Oxygen 6-311G: " << solver.energy() << endl;
    }
}
