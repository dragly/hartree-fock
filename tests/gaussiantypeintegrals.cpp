#include <unittest++/UnitTest++.h>

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

SUITE(GaussianIntegral) {
    TEST(GaussianOverlapIntegralTest) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, 0.3);

        GaussianOverlapIntegral integrator(posA, posB, primitiveA, primitiveB);

        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,0), 0.119172363580852, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,1), 0.276479883507577, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,2), 0.760605693318432, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,0), 0.0429020508891068, 0.00001);
        CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,1), 0.0995327580627279, 0.00001);

    }

    TEST(OverlapNew) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, 0.3);

        GaussianOverlapIntegral integrator(posA, posB, primitiveA, primitiveB);

        CHECK_CLOSE(2.979309089521e-01, integrator.overlapIntegral(2,0,0,2,0,0), 1e-7);
        CHECK_CLOSE(1.072551272228e-02, integrator.overlapIntegral(2,0,0,1,1,0), 1e-7);
        CHECK_CLOSE(6.911997087690e-02, integrator.overlapIntegral(2,0,0,1,0,1), 1e-7);
    }

    TEST(KineticNew) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, 0.3);

        GaussianKineticIntegral integrator(3);

        integrator.set(posA, posB, primitiveA, primitiveB);

        CHECK_CLOSE(-3.468392469657e-01, integrator.kineticIntegral(2,0,0,2,0,0), 1e-7);
        CHECK_CLOSE(-3.562586305833e-03, integrator.kineticIntegral(2,0,0,1,1,0), 1e-7);
        CHECK_CLOSE(-2.295888952647e-02, integrator.kineticIntegral(2,0,0,1,0,1), 1e-7);
        CHECK_CLOSE(2.514020826620e-01, integrator.kineticIntegral(0,0,1,0,0,1), 1e-7);
    }

    TEST(KineticNew2) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};

        GaussianPrimitiveOrbital primitiveA(1.0, 2, 0, 0, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 2, 0, 0, 0.3);

        GaussianKineticIntegral integrator(2);

        integrator.set(posA, posB, primitiveA, primitiveB);

        CHECK_CLOSE(-3.468392469657e-01, integrator.kineticIntegral(2,0,0,2,0,0), 1e-7);
    }

    TEST(GaussianKineticIntegralTest) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, 0.3);

        GaussianKineticIntegral integrator(3);

        integrator.set(posA, posB, primitiveA, primitiveB);

        CHECK_CLOSE(-0.0967870268058250, integrator.kineticIntegral(0,0,0,0,0,0), 0.00001);
        CHECK_CLOSE(-0.158190730147696, integrator.kineticIntegral(0,0,0,0,0,1), 0.00001);
        CHECK_CLOSE(-0.0245468374367114, integrator.kineticIntegral(0,0,0,0,1,0), 0.00001);
        CHECK_CLOSE(-0.0330608009181156, integrator.kineticIntegral(0,0,0,0,1,1), 0.00001);
        CHECK_CLOSE(-0.0681856595464203, integrator.kineticIntegral(0,0,0,1,0,0), 0.00001);
        CHECK_CLOSE(-0.0918355581058773, integrator.kineticIntegral(0,0,0,1,0,1), 0.00001);
        CHECK_CLOSE(-0.0142503452233260, integrator.kineticIntegral(0,0,0,1,1,0), 0.00001);
        CHECK_CLOSE(-0.00917293898306074, integrator.kineticIntegral(0,0,0,1,1,1), 0.00001);
        CHECK_CLOSE(0.251402082662018, integrator.kineticIntegral(0,0,1,0,0,1), 0.00001);
        CHECK_CLOSE(0.0495912013771734, integrator.kineticIntegral(0,0,1,0,1,0), 0.00001);
        CHECK_CLOSE(0.0176714824377215, integrator.kineticIntegral(0,0,1,0,1,1), 0.00001);
        CHECK_CLOSE(0.137753337158815, integrator.kineticIntegral(0,0,1,1,0,0), 0.00001);
        CHECK_CLOSE(0.0490874512158928, integrator.kineticIntegral(0,0,1,1,0,1), 0.00001);
        CHECK_CLOSE(0.0137594084745900, integrator.kineticIntegral(0,0,1,1,1,0), 0.00001);
        CHECK_CLOSE(-0.0551617848828804, integrator.kineticIntegral(0,0,1,1,1,1), 0.00001);
        CHECK_CLOSE(-0.0604904731258242, integrator.kineticIntegral(0,1,0,0,1,0), 0.00001);
        CHECK_CLOSE(-0.0868821710550227, integrator.kineticIntegral(0,1,0,0,1,1), 0.00001);
        CHECK_CLOSE(0.0213755178349884, integrator.kineticIntegral(0,1,0,1,0,0), 0.00001);
        CHECK_CLOSE(0.0137594084745909, integrator.kineticIntegral(0,1,0,1,0,1), 0.00001);
        CHECK_CLOSE(-0.0374492116616455, integrator.kineticIntegral(0,1,0,1,1,0), 0.00001);
        CHECK_CLOSE(-0.0334264444581330, integrator.kineticIntegral(0,1,0,1,1,1), 0.00001);
        CHECK_CLOSE(0.0788748150526478, integrator.kineticIntegral(0,1,1,0,1,1), 0.00001);
        CHECK_CLOSE(-0.0206391127118871, integrator.kineticIntegral(0,1,1,1,0,0), 0.00001);
        CHECK_CLOSE(0.0827426773243232, integrator.kineticIntegral(0,1,1,1,0,1), 0.00001);
    }

    TEST(GaussianColoumbAttractionIntegralTest) {
        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};
        rowvec posC = {2.3, 0.9, 3.2};

        GaussianColoumbAttractionIntegral integrator(3);

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, 0.2);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, 0.3);

        integrator.set(posA, posB, posC, primitiveA, primitiveB);
        CHECK_CLOSE(2.788948987251e-02, integrator.coloumbAttractionIntegral(0,0,0,0,0,0), 1e-4);
        CHECK_CLOSE(6.971203468743e-02, integrator.coloumbAttractionIntegral(0,0,0,0,0,1), 1e-4);
        CHECK_CLOSE(2.024071525839e-01, integrator.coloumbAttractionIntegral(0,0,0,0,0,2), 1e-4);
        CHECK_CLOSE(8.727033700014e-03, integrator.coloumbAttractionIntegral(0,0,0,0,1,0), 1e-4);
        CHECK_CLOSE(2.134361291529e-02, integrator.coloumbAttractionIntegral(0,0,0,0,1,1), 1e-4);
        CHECK_CLOSE(2.921666495443e-02, integrator.coloumbAttractionIntegral(0,0,0,0,2,0), 1e-4);
        CHECK_CLOSE(3.185957751329e-02, integrator.coloumbAttractionIntegral(0,0,0,1,0,0), 1e-4);
        CHECK_CLOSE(8.105746642202e-02, integrator.coloumbAttractionIntegral(0,0,0,1,0,1), 1e-4);
        CHECK_CLOSE(9.596523510045e-03, integrator.coloumbAttractionIntegral(0,0,0,1,1,0), 1e-4);
        CHECK_CLOSE(6.388444040338e-02, integrator.coloumbAttractionIntegral(0,0,0,2,0,0), 1e-4);
        CHECK_CLOSE(-9.204700718547e-02, integrator.coloumbAttractionIntegral(0,0,1,0,0,0), 1e-4);
        CHECK_CLOSE(-2.019226499202e-01, integrator.coloumbAttractionIntegral(0,0,1,0,0,1), 1e-4);
        CHECK_CLOSE(-5.276399683274e-01, integrator.coloumbAttractionIntegral(0,0,1,0,0,2), 1e-4);
        CHECK_CLOSE(-2.927318225377e-02, integrator.coloumbAttractionIntegral(0,0,1,0,1,0), 1e-4);
        CHECK_CLOSE(-6.299787002318e-02, integrator.coloumbAttractionIntegral(0,0,1,0,1,1), 1e-4);
        CHECK_CLOSE(-9.718105595370e-02, integrator.coloumbAttractionIntegral(0,0,1,0,2,0), 1e-4);
        CHECK_CLOSE(-1.037280861539e-01, integrator.coloumbAttractionIntegral(0,0,1,1,0,0), 1e-4);
        CHECK_CLOSE(-2.312309453843e-01, integrator.coloumbAttractionIntegral(0,0,1,1,0,1), 1e-4);
        CHECK_CLOSE(-3.202910466576e-02, integrator.coloumbAttractionIntegral(0,0,1,1,1,0), 1e-4);
        CHECK_CLOSE(-2.073449397904e-01, integrator.coloumbAttractionIntegral(0,0,1,2,0,0), 1e-4);
        CHECK_CLOSE(3.319499900436e-01, integrator.coloumbAttractionIntegral(0,0,2,0,0,0), 1e-4);
        CHECK_CLOSE(6.435114042344e-01, integrator.coloumbAttractionIntegral(0,0,2,0,0,1), 1e-4);
        CHECK_CLOSE(1.536931448007e+00, integrator.coloumbAttractionIntegral(0,0,2,0,0,2), 1e-4);
        CHECK_CLOSE(1.067865861209e-01, integrator.coloumbAttractionIntegral(0,0,2,0,1,0), 1e-4);
        CHECK_CLOSE(2.033153544029e-01, integrator.coloumbAttractionIntegral(0,0,2,0,1,1), 1e-4);
        CHECK_CLOSE(3.524622701603e-01, integrator.coloumbAttractionIntegral(0,0,2,0,2,0), 1e-4);
        CHECK_CLOSE(3.703919381580e-01, integrator.coloumbAttractionIntegral(0,0,2,1,0,0), 1e-4);
        CHECK_CLOSE(7.292169308884e-01, integrator.coloumbAttractionIntegral(0,0,2,1,0,1), 1e-4);
        CHECK_CLOSE(1.162963233448e-01, integrator.coloumbAttractionIntegral(0,0,2,1,1,0), 1e-4);
        CHECK_CLOSE(7.390872806284e-01, integrator.coloumbAttractionIntegral(0,0,2,2,0,0), 1e-4);
        CHECK_CLOSE(-1.637350724302e-02, integrator.coloumbAttractionIntegral(0,1,0,0,0,0), 1e-4);
        CHECK_CLOSE(-4.139721853567e-02, integrator.coloumbAttractionIntegral(0,1,0,0,0,1), 1e-4);
        CHECK_CLOSE(-1.213713540367e-01, integrator.coloumbAttractionIntegral(0,1,0,0,0,2), 1e-4);
        CHECK_CLOSE(2.136233458775e-02, integrator.coloumbAttractionIntegral(0,1,0,0,1,0), 1e-4);
        CHECK_CLOSE(5.306634838230e-02, integrator.coloumbAttractionIntegral(0,1,0,0,1,1), 1e-4);
        CHECK_CLOSE(-1.697963144320e-04, integrator.coloumbAttractionIntegral(0,1,0,0,2,0), 1e-4);
        CHECK_CLOSE(-1.907709578263e-02, integrator.coloumbAttractionIntegral(0,1,0,1,0,0), 1e-4);
        CHECK_CLOSE(-4.932098923684e-02, integrator.coloumbAttractionIntegral(0,1,0,1,0,1), 1e-4);
        CHECK_CLOSE(2.414126847830e-02, integrator.coloumbAttractionIntegral(0,1,0,1,1,0), 1e-4);
        CHECK_CLOSE(-3.842342257899e-02, integrator.coloumbAttractionIntegral(0,1,0,2,0,0), 1e-4);
        CHECK_CLOSE(5.356912442721e-02, integrator.coloumbAttractionIntegral(0,1,1,0,0,0), 1e-4);
        CHECK_CLOSE(1.187325132374e-01, integrator.coloumbAttractionIntegral(0,1,1,0,0,1), 1e-4);
        CHECK_CLOSE(3.128036711345e-01, integrator.coloumbAttractionIntegral(0,1,1,0,0,2), 1e-4);
        CHECK_CLOSE(-7.083519267815e-02, integrator.coloumbAttractionIntegral(0,1,1,0,1,0), 1e-4);
        CHECK_CLOSE(-1.544897767601e-01, integrator.coloumbAttractionIntegral(0,1,1,0,1,1), 1e-4);
        CHECK_CLOSE(-3.894393296797e-04, integrator.coloumbAttractionIntegral(0,1,1,0,2,0), 1e-4);
        CHECK_CLOSE(6.132616876012e-02, integrator.coloumbAttractionIntegral(0,1,1,1,0,0), 1e-4);
        CHECK_CLOSE(1.386353605834e-01, integrator.coloumbAttractionIntegral(0,1,1,1,0,1), 1e-4);
        CHECK_CLOSE(-7.911527548440e-02, integrator.coloumbAttractionIntegral(0,1,1,1,1,0), 1e-4);
        CHECK_CLOSE(1.229240947242e-01, integrator.coloumbAttractionIntegral(0,1,1,2,0,0), 1e-4);
        CHECK_CLOSE(3.609849112824e-02, integrator.coloumbAttractionIntegral(0,2,0,0,0,0), 1e-4);
        CHECK_CLOSE(9.032384498820e-02, integrator.coloumbAttractionIntegral(0,2,0,0,0,1), 1e-4);
        CHECK_CLOSE(2.625292648498e-01, integrator.coloumbAttractionIntegral(0,2,0,0,0,2), 1e-4);
        CHECK_CLOSE(-1.939589748931e-02, integrator.coloumbAttractionIntegral(0,2,0,0,1,0), 1e-4);
        CHECK_CLOSE(-4.913397190183e-02, integrator.coloumbAttractionIntegral(0,2,0,0,1,1), 1e-4);
        CHECK_CLOSE(6.878262296370e-02, integrator.coloumbAttractionIntegral(0,2,0,0,2,0), 1e-4);
        CHECK_CLOSE(4.131065513841e-02, integrator.coloumbAttractionIntegral(0,2,0,1,0,0), 1e-4);
        CHECK_CLOSE(1.052929737663e-01, integrator.coloumbAttractionIntegral(0,2,0,1,0,1), 1e-4);
        CHECK_CLOSE(-2.267402937768e-02, integrator.coloumbAttractionIntegral(0,2,0,1,1,0), 1e-4);
        CHECK_CLOSE(8.289710831960e-02, integrator.coloumbAttractionIntegral(0,2,0,2,0,0), 1e-4);
        CHECK_CLOSE(-3.786414780578e-02, integrator.coloumbAttractionIntegral(1,0,0,0,0,0), 1e-4);
        CHECK_CLOSE(-9.322262110550e-02, integrator.coloumbAttractionIntegral(1,0,0,0,0,1), 1e-4);
        CHECK_CLOSE(-2.671155215998e-01, integrator.coloumbAttractionIntegral(1,0,0,0,0,2), 1e-4);
        CHECK_CLOSE(-1.222106053447e-02, integrator.coloumbAttractionIntegral(1,0,0,0,1,0), 1e-4);
        CHECK_CLOSE(-2.972830178046e-02, integrator.coloumbAttractionIntegral(1,0,0,0,1,1), 1e-4);
        CHECK_CLOSE(-4.026352276293e-02, integrator.coloumbAttractionIntegral(1,0,0,0,2,0), 1e-4);
        CHECK_CLOSE(-1.576450345257e-02, integrator.coloumbAttractionIntegral(1,0,0,1,0,0), 1e-4);
        CHECK_CLOSE(-3.945885129414e-02, integrator.coloumbAttractionIntegral(1,0,0,1,0,1), 1e-4);
        CHECK_CLOSE(-4.918734877201e-03, integrator.coloumbAttractionIntegral(1,0,0,1,1,0), 1e-4);
        CHECK_CLOSE(-2.459437143524e-02, integrator.coloumbAttractionIntegral(1,0,0,2,0,0), 1e-4);
        CHECK_CLOSE(1.263894353489e-01, integrator.coloumbAttractionIntegral(1,0,1,0,0,0), 1e-4);
        CHECK_CLOSE(2.735756798558e-01, integrator.coloumbAttractionIntegral(1,0,1,0,0,1), 1e-4);
        CHECK_CLOSE(7.071773603054e-01, integrator.coloumbAttractionIntegral(1,0,1,0,0,2), 1e-4);
        CHECK_CLOSE(4.115384967585e-02, integrator.coloumbAttractionIntegral(1,0,1,0,1,0), 1e-4);
        CHECK_CLOSE(8.802219023191e-02, integrator.coloumbAttractionIntegral(1,0,1,0,1,1), 1e-4);
        CHECK_CLOSE(1.350111738398e-01, integrator.coloumbAttractionIntegral(1,0,1,0,2,0), 1e-4);
        CHECK_CLOSE(5.197526800952e-02, integrator.coloumbAttractionIntegral(1,0,1,1,0,0), 1e-4);
        CHECK_CLOSE(1.145639876363e-01, integrator.coloumbAttractionIntegral(1,0,1,1,0,1), 1e-4);
        CHECK_CLOSE(1.638641395940e-02, integrator.coloumbAttractionIntegral(1,0,1,1,1,0), 1e-4);
        CHECK_CLOSE(8.278875254192e-02, integrator.coloumbAttractionIntegral(1,0,1,2,0,0), 1e-4);
        CHECK_CLOSE(2.185667243163e-02, integrator.coloumbAttractionIntegral(1,1,0,0,0,0), 1e-4);
        CHECK_CLOSE(5.417205698627e-02, integrator.coloumbAttractionIntegral(1,1,0,0,0,1), 1e-4);
        CHECK_CLOSE(1.560020091608e-01, integrator.coloumbAttractionIntegral(1,1,0,0,0,2), 1e-4);
        CHECK_CLOSE(-2.926456829930e-02, integrator.coloumbAttractionIntegral(1,1,0,0,1,0), 1e-4);
        CHECK_CLOSE(-7.176178735649e-02, integrator.coloumbAttractionIntegral(1,1,0,0,1,1), 1e-4);
        CHECK_CLOSE(-5.223967979758e-04, integrator.coloumbAttractionIntegral(1,1,0,0,2,0), 1e-4);
        CHECK_CLOSE(9.269318129877e-03, integrator.coloumbAttractionIntegral(1,1,0,1,0,0), 1e-4);
        CHECK_CLOSE(2.337071697343e-02, integrator.coloumbAttractionIntegral(1,1,0,1,0,1), 1e-4);
        CHECK_CLOSE(-1.203714316117e-02, integrator.coloumbAttractionIntegral(1,1,0,1,1,0), 1e-4);
        CHECK_CLOSE(1.401501778682e-02, integrator.coloumbAttractionIntegral(1,1,0,2,0,0), 1e-4);
        CHECK_CLOSE(7.889586550718e-02, integrator.coloumbAttractionIntegral(2,0,0,0,0,0), 1e-4);
        CHECK_CLOSE(1.935977010010e-01, integrator.coloumbAttractionIntegral(2,0,0,0,0,1), 1e-4);
        CHECK_CLOSE(5.534914541236e-01, integrator.coloumbAttractionIntegral(2,0,0,0,0,2), 1e-4);
        CHECK_CLOSE(2.563391673303e-02, integrator.coloumbAttractionIntegral(2,0,0,0,1,0), 1e-4);
        CHECK_CLOSE(6.217850538435e-02, integrator.coloumbAttractionIntegral(2,0,0,0,1,1), 1e-4);
        CHECK_CLOSE(8.419480232293e-02, integrator.coloumbAttractionIntegral(2,0,0,0,2,0), 1e-4);
        CHECK_CLOSE(1.481688684288e-02, integrator.coloumbAttractionIntegral(2,0,0,1,0,0), 1e-4);
        CHECK_CLOSE(3.878852644576e-02, integrator.coloumbAttractionIntegral(2,0,0,1,0,1), 1e-4);
        CHECK_CLOSE(4.176920693786e-03, integrator.coloumbAttractionIntegral(2,0,0,1,1,0), 1e-4);
        CHECK_CLOSE(6.422210627967e-02, integrator.coloumbAttractionIntegral(2,0,0,2,0,0), 1e-4);
    }

    TEST(GaussianElectronInteractionIntegralTest1) {
        rowvec posA = {-0.5, 0, 0};
        rowvec posB = {-0.5, 0, 0};
        rowvec posC = {-0.5, 0, 0};
        rowvec posD = {-0.5, 0, 0};
        double a = 13.0077;
        double b = 13.0077;
        double c = 13.0077;
        double d = 13.0077;

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, a);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, b);
        GaussianPrimitiveOrbital primitiveC(1.0, 3, 3, 3, c);
        GaussianPrimitiveOrbital primitiveD(1.0, 3, 3, 3, d);

        GaussianElectronInteractionIntegral integrator(3);
        integrator.set(posA, posB, posC, posD, primitiveA, primitiveB, primitiveC, primitiveD);
        // regression test
        CHECK_CLOSE(0.0071666040410096028615, integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0), 1e-9);
    }

    TEST(GaussianElectronInteractionIntegralTest2) {
        rowvec posA = {0.5, 0, 0};
        rowvec posB = {-0.5, 0, 0};
        rowvec posC = {-0.5, 0, 0};
        rowvec posD = {0.5, 0, 0};
        double a = 13.0077;
        double b = 0.121949;
        double c = 0.444529;
        double d = 13.0077;

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, a);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, b);
        GaussianPrimitiveOrbital primitiveC(1.0, 3, 3, 3, c);
        GaussianPrimitiveOrbital primitiveD(1.0, 3, 3, 3, d);

        GaussianElectronInteractionIntegral integrator(3);
        integrator.set(posA, posB, posC, posD, primitiveA, primitiveB, primitiveC, primitiveD);
        // regression test
        CHECK_CLOSE(0.022124581472837051566, integrator.electronInteractionIntegral(0,0,0,0,0,0,0,0,0,0,0,0), 1e-9);
    }

    TEST(GaussianElectronInteractionIntegralTest3) {
        rowvec posA = {0.5, 0, 0};
        rowvec posB = {-0.5, 0, 0};
        rowvec posC = {-0.5, 0, 0};
        rowvec posD = {0.5, 0, 0};
        double a = 13.0077;
        double b = 0.121949;
        double c = 0.444529;
        double d = 13.0077;

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, a);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, b);
        GaussianPrimitiveOrbital primitiveC(1.0, 3, 3, 3, c);
        GaussianPrimitiveOrbital primitiveD(1.0, 3, 3, 3, d);

        GaussianElectronInteractionIntegral integrator(3);
        integrator.set(posA, posB, posC, posD, primitiveA, primitiveB, primitiveC, primitiveD);
        // regression test
        CHECK_CLOSE(0.0001385810300677682, integrator.electronInteractionIntegral(0,0,0,0,1,0,0,1,0,0,0,0), 1e-9);
    }

    TEST(GaussianElectronInteractionIntegralTest4) {
        rowvec posA = {0.55, 1, 3};
        rowvec posB = {-0.52, 5, 6};
        rowvec posC = {-0.53, 1, 2};
        rowvec posD = {0.45, 2, 4};
        double a = 13.0077;
        double b = 0.121949;
        double c = 0.444529;
        double d = 10.0077;

        GaussianPrimitiveOrbital primitiveA(1.0, 3, 3, 3, a);
        GaussianPrimitiveOrbital primitiveB(1.0, 3, 3, 3, b);
        GaussianPrimitiveOrbital primitiveC(1.0, 3, 3, 3, c);
        GaussianPrimitiveOrbital primitiveD(1.0, 3, 3, 3, d);

        GaussianElectronInteractionIntegral integrator(3);
        integrator.set(posA, posB, posC, posD, primitiveA, primitiveB, primitiveC, primitiveD);
        // regression test
        CHECK_CLOSE(-6.8145328932903484228e-08, integrator.electronInteractionIntegral(1,0,0,0,1,0,0,1,0,0,1,0), 1e-9);
    }
}
