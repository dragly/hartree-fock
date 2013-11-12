#include <unittest++/UnitTest++.h>

#include <parsers/turbomoleparser.h>

SUITE(Parser) {
    TEST(Parser) {
        TurboMoleParser parser;
        parser.read("oxygen431g.tm");
    }
}
