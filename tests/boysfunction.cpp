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
    TEST(BoysfactorialFunctions)
    {
//        Boys boys(0);
        CHECK_EQUAL(BoysFunction::factorial(0),1);
        CHECK_EQUAL(BoysFunction::factorial(1),1);
        CHECK_EQUAL(BoysFunction::factorial(5),120);
        CHECK_EQUAL(BoysFunction::factorial(20),2432902008176640000);
        CHECK_EQUAL(BoysFunction::factorial(-1),1);

        CHECK_EQUAL(BoysFunction::doubleFactorial(0),1);
        CHECK_EQUAL(BoysFunction::doubleFactorial(1),1);
        CHECK_EQUAL(BoysFunction::doubleFactorial(5),15);
        CHECK_EQUAL(BoysFunction::doubleFactorial(-1),1);
    }
    TEST(boysFunction)
    {
        BoysFunctionIntermediate boysIntermediate;
        rowvec xvec = linspace<rowvec>(0.01,50,20);
        mat F0 = zeros(20,1);
        for(uint i = 0; i < 20; i++){
            BoysFunction boys(xvec[i], 0, &boysIntermediate);
            F0(i,0) = boys.result(0);
        }

        CHECK_CLOSE(F0(0,0), 9.9667664290336333e-01, 1e-10);
        CHECK_CLOSE(F0(1,0), 5.3357683580906246e-01, 1e-10);
        CHECK_CLOSE(F0(2,0), 3.8551956882369670e-01, 1e-10);
        CHECK_CLOSE(F0(3,0), 3.1522027004511449e-01, 1e-10);
        CHECK_CLOSE(F0(4,0), 2.7304989839411259e-01, 1e-10);
        CHECK_CLOSE(F0(5,0), 2.4424745291117503e-01, 1e-10);
        CHECK_CLOSE(F0(6,0), 2.2298057384240719e-01, 1e-10);
        CHECK_CLOSE(F0(7,0), 2.0644923632081960e-01, 1e-10);
        CHECK_CLOSE(F0(8,0), 1.9312212797082015e-01, 1e-10);
        CHECK_CLOSE(F0(9,0), 1.8208209208151274e-01, 1e-10);
        CHECK_CLOSE(F0(10,0), 1.7274188563417292e-01, 1e-10);
        CHECK_CLOSE(F0(11,0), 1.6470576997839928e-01, 1e-10);
        CHECK_CLOSE(F0(12,0), 1.5769603853464209e-01, 1e-10);
        CHECK_CLOSE(F0(13,0), 1.5151129820345471e-01, 1e-10);
        CHECK_CLOSE(F0(14,0), 1.4600146419605400e-01, 1e-10);
        CHECK_CLOSE(F0(15,0), 1.4105209097756166e-01, 1e-10);
        CHECK_CLOSE(F0(16,0), 1.3657418098484136e-01, 1e-10);
        CHECK_CLOSE(F0(17,0), 1.3249734289737117e-01, 1e-10);
        CHECK_CLOSE(F0(18,0), 1.2876507161025572e-01, 1e-10);
        CHECK_CLOSE(F0(19,0), 1.2533141373155002e-01, 1e-10);




        //test for F0(x) for large x > 50:
        //Recursion relation is not used!
//        Boys boysF0_large(0);
        xvec = linspace<rowvec>(51,100,20);
        F0 = zeros(20,1);
        for(uint i = 0; i < 20; i++){
            BoysFunction boys(xvec[i], 0, &boysIntermediate);
//            boysF0_large.evaluateBoysFunctions(xvec[i]);
            F0(i,0) = boys.result(0);
        }

        CHECK_CLOSE(F0(0,0), 1.2409659136408727e-01, 1e-10);
        CHECK_CLOSE(F0(1,0), 1.2107315290425182e-01, 1e-10);
        CHECK_CLOSE(F0(2,0), 1.1826045114934411e-01, 1e-10);
        CHECK_CLOSE(F0(3,0), 1.1563509029711123e-01, 1e-10);
        CHECK_CLOSE(F0(4,0), 1.1317715652582763e-01, 1e-10);
        CHECK_CLOSE(F0(5,0), 1.1086957884293665e-01, 1e-10);
        CHECK_CLOSE(F0(6,0), 1.0869762771661817e-01, 1e-10);
        CHECK_CLOSE(F0(7,0), 1.0664851770975842e-01, 1e-10);
        CHECK_CLOSE(F0(8,0), 1.0471108954515591e-01, 1e-10);
        CHECK_CLOSE(F0(9,0), 1.0287555349437555e-01, 1e-10);
        CHECK_CLOSE(F0(10,0), 1.0113328058446551e-01, 1e-10);
        CHECK_CLOSE(F0(11,0), 9.9476631436519997e-02, 1e-10);
        CHECK_CLOSE(F0(12,0), 9.7898814974335960e-02, 1e-10);
        CHECK_CLOSE(F0(13,0), 9.6393771031862641e-02, 1e-10);
        CHECK_CLOSE(F0(14,0), 9.4956072224460661e-02, 1e-10);
        CHECK_CLOSE(F0(15,0), 9.3580841456184838e-02, 1e-10);
        CHECK_CLOSE(F0(16,0), 9.2263682201430275e-02, 1e-10);
        CHECK_CLOSE(F0(17,0), 9.1000619287057175e-02, 1e-10);
        CHECK_CLOSE(F0(18,0), 8.9788048355705335e-02, 1e-10);
        CHECK_CLOSE(F0(19,0), 8.8622692545275800e-02, 1e-10);
    }

    TEST(BoysDownwardrecursionFunction)
    {

        BoysFunctionIntermediate boysIntermediate;
        //test for Fn(0):
        BoysFunction boysFn_0(0, 20, &boysIntermediate);
//        boysFn_0.evaluateBoysFunctions(0);
        rowvec Fn = zeros<rowvec>(20);
        for(int i = 0; i < 20; i++) {
            Fn(i) = boysFn_0.result(i);
        }
        CHECK_CLOSE(Fn[0], 1.0000000000000000e+00, 1e-10);
        CHECK_CLOSE(Fn[1], 3.3333333333333331e-01, 1e-10);
        CHECK_CLOSE(Fn[2], 2.0000000000000001e-01, 1e-10);
        CHECK_CLOSE(Fn[3], 1.4285714285714285e-01, 1e-10);
        CHECK_CLOSE(Fn[4], 1.1111111111111110e-01, 1e-10);
        CHECK_CLOSE(Fn[5], 9.0909090909090912e-02, 1e-10);
        CHECK_CLOSE(Fn[6], 7.6923076923076927e-02, 1e-10);
        CHECK_CLOSE(Fn[7], 6.6666666666666666e-02, 1e-10);
        CHECK_CLOSE(Fn[8], 5.8823529411764705e-02, 1e-10);
        CHECK_CLOSE(Fn[9], 5.2631578947368418e-02, 1e-10);
        CHECK_CLOSE(Fn[10], 4.7619047619047616e-02, 1e-10);
        CHECK_CLOSE(Fn[11], 4.3478260869565216e-02, 1e-10);
        CHECK_CLOSE(Fn[12], 4.0000000000000001e-02, 1e-10);
        CHECK_CLOSE(Fn[13], 3.7037037037037035e-02, 1e-10);
        CHECK_CLOSE(Fn[14], 3.4482758620689655e-02, 1e-10);
        CHECK_CLOSE(Fn[15], 3.2258064516129031e-02, 1e-10);
        CHECK_CLOSE(Fn[16], 3.0303030303030304e-02, 1e-10);
        CHECK_CLOSE(Fn[17], 2.8571428571428571e-02, 1e-10);
        CHECK_CLOSE(Fn[18], 2.7027027027027029e-02, 1e-10);
        CHECK_CLOSE(Fn[19], 2.5641025641025640e-02, 1e-10);
//        CHECK_CLOSE(Fn[20], 2.4390243902439025e-02, 1e-9);




        //test for F0(x) for small x <= 50:
        // n goes from nMax to 0
        mat F0; rowvec xvec;
        xvec = linspace<rowvec>(0.01,50,20);
        F0 = zeros(20,16);
        for(uint i = 0; i < 20; i++){
            BoysFunction boysF0_small(xvec[i], 15, &boysIntermediate);
//            boysF0_small.evaluateBoysFunctions(xvec[i]);
            F0(i,0) = boysF0_small.result(0);
        }

        CHECK_CLOSE(F0(0,0), 9.9667664290336333e-01, 1e-10);
        CHECK_CLOSE(F0(1,0), 5.3357683580906246e-01, 1e-10);
        CHECK_CLOSE(F0(2,0), 3.8551956882369670e-01, 1e-10);
        CHECK_CLOSE(F0(3,0), 3.1522027004511449e-01, 1e-10);
        CHECK_CLOSE(F0(4,0), 2.7304989839411259e-01, 1e-10);
        CHECK_CLOSE(F0(5,0), 2.4424745291117503e-01, 1e-10);
        CHECK_CLOSE(F0(6,0), 2.2298057384240719e-01, 1e-10);
        CHECK_CLOSE(F0(7,0), 2.0644923632081960e-01, 1e-10);
        CHECK_CLOSE(F0(8,0), 1.9312212797082015e-01, 1e-10);
        CHECK_CLOSE(F0(9,0), 1.8208209208151274e-01, 1e-10);
        CHECK_CLOSE(F0(10,0), 1.7274188563417292e-01, 1e-10);
        CHECK_CLOSE(F0(11,0), 1.6470576997839928e-01, 1e-10);
        CHECK_CLOSE(F0(12,0), 1.5769603853464209e-01, 1e-10);
        CHECK_CLOSE(F0(13,0), 1.5151129820345471e-01, 1e-10);
        CHECK_CLOSE(F0(14,0), 1.4600146419605400e-01, 1e-10);
        CHECK_CLOSE(F0(15,0), 1.4105209097756166e-01, 1e-10);
        CHECK_CLOSE(F0(16,0), 1.3657418098484136e-01, 1e-10);
        CHECK_CLOSE(F0(17,0), 1.3249734289737117e-01, 1e-10);
        CHECK_CLOSE(F0(18,0), 1.2876507161025572e-01, 1e-10);
        CHECK_CLOSE(F0(19,0), 1.2533141373155002e-01, 1e-10);


        //test for F0(x) for large x > 50:
        // n goes from nMax to 0
        xvec = linspace<rowvec>(51,100,20);
        F0 = zeros(20,16);
        for(uint i = 0; i < 20; i++){
            BoysFunction boysF0_large(xvec[i], 15, &boysIntermediate);
//             boysF0_large.evaluateBoysFunctions(xvec[i]);
            F0(i,0) = boysF0_large.result(0);
        }

        CHECK_CLOSE(F0(0,0), 1.2409659136408727e-01, 1e-9); //OBS!!
        CHECK_CLOSE(F0(1,0), 1.2107315290425182e-01, 1e-10);
        CHECK_CLOSE(F0(2,0), 1.1826045114934411e-01, 1e-10);
        CHECK_CLOSE(F0(3,0), 1.1563509029711123e-01, 1e-10);
        CHECK_CLOSE(F0(4,0), 1.1317715652582763e-01, 1e-10);
        CHECK_CLOSE(F0(5,0), 1.1086957884293665e-01, 1e-10);
        CHECK_CLOSE(F0(6,0), 1.0869762771661817e-01, 1e-10);
        CHECK_CLOSE(F0(7,0), 1.0664851770975842e-01, 1e-10);
        CHECK_CLOSE(F0(8,0), 1.0471108954515591e-01, 1e-10);
        CHECK_CLOSE(F0(9,0), 1.0287555349437555e-01, 1e-10);
        CHECK_CLOSE(F0(10,0), 1.0113328058446551e-01, 1e-10);
        CHECK_CLOSE(F0(11,0), 9.9476631436519997e-02, 1e-10);
        CHECK_CLOSE(F0(12,0), 9.7898814974335960e-02, 1e-10);
        CHECK_CLOSE(F0(13,0), 9.6393771031862641e-02, 1e-10);
        CHECK_CLOSE(F0(14,0), 9.4956072224460661e-02, 1e-10);
        CHECK_CLOSE(F0(15,0), 9.3580841456184838e-02, 1e-10);
        CHECK_CLOSE(F0(16,0), 9.2263682201430275e-02, 1e-10);
        CHECK_CLOSE(F0(17,0), 9.1000619287057175e-02, 1e-10);
        CHECK_CLOSE(F0(18,0), 8.9788048355705335e-02, 1e-10);
        CHECK_CLOSE(F0(19,0), 8.8622692545275800e-02, 1e-10);
    }
}
