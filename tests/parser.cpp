#include <unittest++/UnitTest++.h>

#include <parsers/turbomoleparser.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include "solvers/restrictedhartreefocksolver.h"

SUITE(Parser) {
    TEST(Parser) {
        TurboMoleParser parser;
        parser.load("oxygen431g.tm");
//        cout << " --------------- " << endl;
//        GaussianNitrogen431G ox;
    }
}
