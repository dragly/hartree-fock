#include <unittest++/UnitTest++.h>

#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

SUITE(CoulombIntegrals) {
//    TEST(GTOelectronElectronIntegral)
//    {
//        GaussianElectronInteractionIntegral integrator(2);

//        rowvec posA = {1.2,2.3,3.4};
//        rowvec posB = {-1.3,1.4,-2.4};
//        rowvec posC = {2.3,0.9,3.2};
//        rowvec posD = {5.0,1.9,1.2};

//        integrator.setAB(posA, posB, 0.2, 0.3);
//        integrator.setCD(posC, posD, 0.4, 0.1);

//        double value = integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0);
//        CHECK_CLOSE(1.624848e-01, value, 1e-5);

//    //    primitiveA.setPowers({0,0,0});
//    //    primitiveB.setPowers({1,0,0});
//    //    primitiveC.setPowers({0,0,0});
//    //    primitiveD.setPowers({0,0,1});
//    //    integrator.setPrimitiveA(primitiveA);
//    //    integrator.setPrimitiveB(primitiveB);
//    //    integrator.setPrimitiveC(primitiveC);
//    //    integrator.setPrimitiveD(primitiveD);
//    //    CHECK_CLOSE(3.297371e+01, integrator.electronRepulsionIntegral(),
//    //1e-5);

//        value = integrator.electronInteractionIntegral(0,0,0,1,0,0,0,0,0,0,0,1);
//        CHECK_CLOSE(3.297371e+01, value, 1e-5);

//    //    primitiveA.setPowers({0,0,0});
//    //    primitiveB.setPowers({1,0,0});
//    //    primitiveC.setPowers({0,2,0});
//    //    primitiveD.setPowers({0,0,1});
//    //    integrator.setPrimitiveA(primitiveA);
//    //    integrator.setPrimitiveB(primitiveB);
//    //    integrator.setPrimitiveC(primitiveC);
//    //    integrator.setPrimitiveD(primitiveD);
//    //    CHECK_CLOSE(3.333488e+01, integrator.electronRepulsionIntegral(),
//        //1e-5);

//        value = integrator.electronInteractionIntegral(0,0,0,1,0,0,0,2,0,0,0,1);
//        CHECK_CLOSE(3.333488e+01, value, 1e-5);


//    //    primitiveA.setPowers({0,0,1});
//    //    primitiveB.setPowers({1,0,0});
//    //    primitiveC.setPowers({0,2,0});
//    //    primitiveD.setPowers({0,0,1});
//    //    integrator.setPrimitiveA(primitiveA);
//    //    integrator.setPrimitiveB(primitiveB);
//    //    integrator.setPrimitiveC(primitiveC);
//    //    integrator.setPrimitiveD(primitiveD);
//    //    CHECK_CLOSE(-1.109942e+02, integrator.electronRepulsionIntegral(),
//    //1e-5);

//        value = integrator.electronInteractionIntegral(0,0,1,1,0,0,0,2,0,0,0,1);
//        CHECK_CLOSE(3.333488e+01, value, 1e-5);



//    }
}
