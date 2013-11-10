#include "simplegaussian.h"

SimpleGaussian::SimpleGaussian() :
    GaussianSystem()
{
    m_nParticles = 4;
    m_nBasisFunctions = 4;
    m_angularMomentumMax = 0;
    electronInteractionIntegral = GaussianElectronInteractionIntegral(m_angularMomentumMax);
    m_coreCharge = 1;

    rowvec corePositions = {-2, 0, 0,
                            2, 0, 0,
                            0, -2, 0,
                            0, 2, 0};
    m_corePositions = corePositions;
    m_corePositions.reshape(3,corePositions.n_elem / 3);
    m_corePositions = m_corePositions.t();
    cout << m_corePositions << endl;

    for(uint i = 0; i < m_corePositions.n_rows; i++) {
        GaussianContractedOrbital contracted1(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive1(5.0, 0, 0, 0, 2.3);
        contracted1.addPrimitiveBasisFunction(primitive1);
        m_basisFunctions.push_back(contracted1);
    }
}
