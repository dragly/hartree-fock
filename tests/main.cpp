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
    std::vector<std::string> tests;
    if(argc > 1 && !strcmp(argv[1], "dev")) {
        tests.push_back("Development");
    } else {
        tests.push_back("CoulombIntegrals");
        tests.push_back("BoysFunction");
        tests.push_back("GaussianIntegral");
        tests.push_back("Hydrogen");
        tests.push_back("GaussianOxygen");
        tests.push_back("Parser");
        tests.push_back("Systems");
        tests.push_back("Unrestricted");
    }
    int result = 0;
    for(const std::string& testName : tests) {
        std::cout << "Running " << testName << std::endl;
        UnitTest::TestReporterStdout reporter;
        UnitTest::TestRunner runner(reporter);
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), testName.c_str(), UnitTest::True(), 0);
    }
    if(result != 0) {
        std::cerr << "FAILURE: Some tests failed! See above log for details..." << std::endl;
    }
    return result;
}

