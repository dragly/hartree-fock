#include <unittest++/UnitTest++.h>
#include <math/boysfunction.h>
#include <math/boysfunctionintermediate.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

SUITE(BoysFunction) {
    TEST(BoysSingleResult) {
        BoysFunctionIntermediate intermediate(5);
        BoysFunction boys(3.13137, 5, &intermediate);
        // Regression test
        CHECK_CLOSE(0.007109579000859712, boys.result(5), 1e-9);
    }
    TEST(BoysIntermediateTest) {
        BoysFunctionIntermediate boysIntermediate(20);
        boysIntermediate.updateResults();
        // Regression test
        CHECK_CLOSE(0.0002310081138823203, boysIntermediate.result(6, 10), 1e-9);
    }
}
