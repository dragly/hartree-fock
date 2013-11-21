#include "gaussianoxygen431g.h"

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

#include <iostream>

using namespace std;

GaussianOxygen431G::GaussianOxygen431G() :
    GaussianSystem()
{
    m_nParticles = 16;
    m_nBasisFunctions = 2 * 9;
    m_angularMomentumMax = 1;
    electronInteractionIntegral = GaussianElectronInteractionIntegral(m_angularMomentumMax);
    m_coreCharge = 8;

    rowvec corePositions = {-1.14, 0, 0,
                            1.14, 0, 0};
    m_corePositions = corePositions;
    m_corePositions.reshape(3,2);
    m_corePositions = m_corePositions.t();

    for(uint i = 0; i < m_corePositions.n_rows; i++) {
        GaussianContractedOrbital contracted1_1s(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive1_1s1(0.0175506 * pow(2 * 883.2728600 / M_PI, 0.75), 0, 0, 0, 883.2728600);
        GaussianPrimitiveOrbital primitive1_1s2(0.1228292 * pow(2 * 133.1292800 / M_PI, 0.75), 0, 0, 0, 133.1292800);
        GaussianPrimitiveOrbital primitive1_1s3(0.4348836 * pow(2 * 29.9064080 / M_PI, 0.75), 0, 0, 0, 29.9064080);
        GaussianPrimitiveOrbital primitive1_1s4(0.5600108 * pow(2 * 7.9786772 / M_PI, 0.75), 0, 0, 0, 7.9786772);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s1);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s2);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s3);
        contracted1_1s.addPrimitiveBasisFunction(primitive1_1s4);
        m_basisFunctions.push_back(contracted1_1s);

        GaussianContractedOrbital contracted2_2s(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive2_2s1(-0.1134010 * pow(2 * 16.1944470 / M_PI, 0.75), 0, 0, 0, 16.1944470);
        GaussianPrimitiveOrbital primitive2_2s2(-0.1772865 * pow(2 * 3.7800860 / M_PI, 0.75), 0, 0, 0, 3.7800860);
        GaussianPrimitiveOrbital primitive2_2s3(1.1504079 * pow(2 * 1.0709836 / M_PI, 0.75), 0, 0, 0, 1.0709836);
        contracted2_2s.addPrimitiveBasisFunction(primitive2_2s1);
        contracted2_2s.addPrimitiveBasisFunction(primitive2_2s2);
        contracted2_2s.addPrimitiveBasisFunction(primitive2_2s3);
        m_basisFunctions.push_back(contracted2_2s);

        for(int dim = 0; dim < 3; dim++) {
            GaussianContractedOrbital contracted3_2p(m_corePositions.row(i));
            GaussianPrimitiveOrbital primitive3_2p1(0.0685453 * pow(2 * 16.1944470 / M_PI, 0.75) * 2 * sqrt(16.1944470), (dim == 0), (dim == 1), (dim == 2), 16.1944470);
            GaussianPrimitiveOrbital primitive3_2p2(0.3312254 * pow(2 * 3.7800860 / M_PI, 0.75) * 2 * sqrt(3.7800860), (dim == 0), (dim == 1), (dim == 2), 3.7800860);
            GaussianPrimitiveOrbital primitive3_2p3(0.7346079 * pow(2 * 1.0709836/ M_PI, 0.75) * 2 * sqrt(1.0709836), (dim == 0), (dim == 1), (dim == 2), 1.0709836);
            contracted3_2p.addPrimitiveBasisFunction(primitive3_2p1);
            contracted3_2p.addPrimitiveBasisFunction(primitive3_2p2);
            contracted3_2p.addPrimitiveBasisFunction(primitive3_2p3);
            m_basisFunctions.push_back(contracted3_2p);
        }

        GaussianContractedOrbital contracted4_2s(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive4_2s(1.0000000 * pow(2 * 0.2838798 / M_PI, 0.75), 0, 0, 0, 0.2838798);
        contracted4_2s.addPrimitiveBasisFunction(primitive4_2s);
        m_basisFunctions.push_back(contracted4_2s);

        for(int dim = 0; dim < 3; dim++) {
            GaussianContractedOrbital contracted5_2p(m_corePositions.row(i));
            GaussianPrimitiveOrbital primitive5_2p(1.0000000 * pow(2 * 0.2838798 / M_PI, 0.75) * 2 * sqrt(0.2838798), (dim == 0), (dim == 1), (dim == 2), 0.2838798);
            contracted5_2p.addPrimitiveBasisFunction(primitive5_2p);
            m_basisFunctions.push_back(contracted5_2p);
        }

    }
}


