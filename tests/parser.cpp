#include <unittest++/UnitTest++.h>

#include <parsers/turbomoleparser.h>
#include <electronsystems/gaussian/gaussianoxygen431g.h>
#include <electronsystems/gaussian/gaussiannitrogen431g.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <hartreefocksolver.h>

SUITE(Parser) {
    TEST(Parser) {
        TurboMoleParser parser;
        parser.load("oxygen431g.tm");
//        cout << " --------------- " << endl;
//        GaussianNitrogen431G ox;
    }
}
