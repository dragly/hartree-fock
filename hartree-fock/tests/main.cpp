#include <boysfunction.cpp>
#include <gaussiantypeintegrals.cpp>
#include <hydrogen.cpp>
#include <gaussianoxygen.cpp>

#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>

#include <memory>

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

class TestMe {
public:
    void banana() {
        std::cout << "Banana!" << std::endl;
    }
};

int main()
{
//    mat test;
//    test.load("boys_function_data_nx1000_limmin_0_limmax_50_nt1000000.arma");
//    test.save("boys_function_data_nx1000_limmin_0_limmax_50_nt1000000.dat", raw_ascii);

    int result = 0;
#ifdef RUN_DEVELOPMENT_TESTS
    ::testing::InitGoogleMock(&argc, argv);
    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);
    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Development", UnitTest::True(), 0);
#else
    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);
    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "BoysFunction", UnitTest::True(), 0);
//    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "GaussianIntegral", UnitTest::True(), 0);
//    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "Hydrogen", UnitTest::True(), 0);
//    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "GaussianOxygen", UnitTest::True(), 0);
#endif

    return result;
}

