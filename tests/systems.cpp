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
}
