#include "gaussiannitrogen431g.h"

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

#include <iostream>

using namespace std;

GaussianNitrogen431G::GaussianNitrogen431G() :
    GaussianSystem()
{
    m_nParticles = 14;
    m_nBasisFunctions = 2 * 9;
    m_angularMomentumMax = 1;
    electronInteractionIntegral = GaussianElectronInteractionIntegral(m_angularMomentumMax);
    m_coreCharge = 7;

    rowvec corePositions = {-1.025, 0, 0,
                            1.025, 0, 0};
    m_corePositions = corePositions;
    m_corePositions.reshape(3,2);
    m_corePositions = m_corePositions.t();

    for(uint i = 0; i < m_corePositions.n_rows; i++) {
        GaussianContractedOrbital contracted1_1s(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive1_1s1(0.0175982511 * pow(2 * 671.2795000 / M_PI, 0.75), 0, 0, 0, 671.2795000);
        GaussianPrimitiveOrbital primitive1_1s2(0.1228462410 * pow(2 * 101.2017000 / M_PI, 0.75), 0, 0, 0, 101.2017000);
        GaussianPrimitiveOrbital primitive1_1s3(0.4337821410 * pow(2 * 22.6999700 / M_PI, 0.75), 0, 0, 0, 22.6999700);
        GaussianPrimitiveOrbital primitive1_1s4(0.5614182170 * pow(2 * 6.0406090 / M_PI, 0.75), 0, 0, 0, 6.0406090);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s1);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s2);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s3);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s4);
        m_basisFunctions.push_back(contracted1_1s);

        GaussianContractedOrbital contracted2_2s(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive2_2s1(-0.1174892990 * pow(2 * 12.3935997 / M_PI, 0.75), 0, 0, 0, 12.3935997);
        GaussianPrimitiveOrbital primitive2_2s2(-0.2139940160 * pow(2 * 2.9223828 / M_PI, 0.75), 0, 0, 0, 2.9223828);
        GaussianPrimitiveOrbital primitive2_2s3(1.1745021100 * pow(2 * 0.83252808 / M_PI, 0.75), 0, 0, 0, 0.83252808);
        contracted2_2s.addPrimitiveBasisFunction(primitive2_2s1);
        contracted2_2s.addPrimitiveBasisFunction(primitive2_2s2);
        contracted2_2s.addPrimitiveBasisFunction(primitive2_2s3);
        m_basisFunctions.push_back(contracted2_2s);

        for(int dim = 0; dim < 3; dim++) {
            GaussianContractedOrbital contracted3_2p(m_corePositions.row(i));
            GaussianPrimitiveOrbital primitive3_2p1(0.0640203443 * pow(2 * 12.3935997/ M_PI, 0.75) * 2 * sqrt(12.3935997), (dim == 0), (dim == 1), (dim == 2), 12.3935997);
            GaussianPrimitiveOrbital primitive3_2p2(0.3112025550 * pow(2 * 2.9223828 / M_PI, 0.75) * 2 * sqrt(2.9223828), (dim == 0), (dim == 1), (dim == 2), 2.9223828);
            GaussianPrimitiveOrbital primitive3_2p3(0.7527482390 * pow(2 * 0.83252808/ M_PI, 0.75) * 2 * sqrt(0.83252808), (dim == 0), (dim == 1), (dim == 2), 0.83252808);
            contracted3_2p.addPrimitiveBasisFunction(primitive3_2p1);
            contracted3_2p.addPrimitiveBasisFunction(primitive3_2p2);
            contracted3_2p.addPrimitiveBasisFunction(primitive3_2p3);
            m_basisFunctions.push_back(contracted3_2p);
        }

        GaussianContractedOrbital contracted4_2s(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive4_2s(1.0000000 * pow(2 * 0.2259640 / M_PI, 0.75), 0, 0, 0, 0.2259640);
        contracted4_2s.addPrimitiveBasisFunction(primitive4_2s);
        m_basisFunctions.push_back(contracted4_2s);

        for(int dim = 0; dim < 3; dim++) {
            GaussianContractedOrbital contracted5_2p(m_corePositions.row(i));
            GaussianPrimitiveOrbital primitive5_2p(1.0000000 * pow(2 * 0.2259640 / M_PI, 0.75) * 2 * sqrt(0.2259640), (dim == 0), (dim == 1), (dim == 2), 0.2259640);
            contracted5_2p.addPrimitiveBasisFunction(primitive5_2p);
            m_basisFunctions.push_back(contracted5_2p);
        }

    }
    for(GaussianContractedOrbital orbital : m_basisFunctions) {
        cout << "Contracted:" << endl;
        for(GaussianPrimitiveOrbital primitive : orbital.primitiveBasisFunctions()) {
            cout << primitive.weight() << " * " << primitive.exponent() << endl;
        }
    }
}



