#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include "solvers/restrictedhartreefocksolver.h"

#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>

#include <memory>

int main(int argc, char* argv[])
{
    int result = 0;
    if(argc > 1 && !strcmp(argv[1], "dev")) {
        UnitTest::TestReporterStdout reporter;
        UnitTest::TestRunner runner(reporter);

        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Development", UnitTest::True(), 0);
    } else {
        UnitTest::TestReporterStdout reporter;
        UnitTest::TestRunner runner(reporter);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "CoulombIntegrals", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "BoysFunction", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "GaussianIntegral", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Hydrogen", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "GaussianOxygen", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Parser", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Systems", UnitTest::True(), 0);
        result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Unrestricted", UnitTest::True(), 0);
    }
    if(result != 0) {
        std::cerr << "FAILURE: Some tests failed! See above log for details..." << std::endl;
    }
    return result;
}

