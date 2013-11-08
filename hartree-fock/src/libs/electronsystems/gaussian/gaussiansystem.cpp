#include "gaussiansystem.h"

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

#include <iostream>

using namespace std;

GaussianSystem::GaussianSystem() :
    m_nParticles(0),
    m_nBasisFunctions(0),
    m_angularMomentumMax(0),
    electronInteractionIntegral(0),
    m_coreCharge(0)
{
}

double GaussianSystem::coupledIntegral(int p, int q, int r, int s)
{
    double result = 0;
    const GaussianContractedOrbital& pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    const GaussianContractedOrbital& rBF = m_basisFunctions.at(r);
    const GaussianContractedOrbital& sBF = m_basisFunctions.at(s);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
            electronInteractionIntegral.setAB(pBF.corePosition(), qBF.corePosition(), pP.exponent(), qP.exponent());
            for(const GaussianPrimitiveOrbital& rP : rBF.primitiveBasisFunctions()) {
                for(const GaussianPrimitiveOrbital& sP : sBF.primitiveBasisFunctions()) {
                    electronInteractionIntegral.setCD(rBF.corePosition(), sBF.corePosition(), rP.exponent(), sP.exponent());
                    result += pP.weight() * rP.weight() * qP.weight() * sP.weight()
                            * electronInteractionIntegral.electronInteractionIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                                      qP.xExponent(), qP.yExponent(), qP.zExponent(),
                                                                                      rP.xExponent(), rP.yExponent(), rP.zExponent(),
                                                                                      sP.xExponent(), sP.yExponent(), sP.zExponent());
                }
            }
        }
    }
    return result;
}

double GaussianSystem::uncoupledIntegral(int p, int q)
{
    double result = 0;
    const GaussianContractedOrbital& pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
            GaussianKineticIntegral kineticIntegral(pBF.corePosition(), qBF.corePosition(),
                                                    pP.exponent(), qP.exponent(),
                                                    m_angularMomentumMax);
            result += pP.weight() * qP.weight() * kineticIntegral.kineticIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                                  qP.xExponent(), qP.yExponent(), qP.zExponent());
            for(uint i = 0; i < m_corePositions.n_rows; i++) {
                const rowvec &corePositionC = m_corePositions.row(i);
                GaussianColoumbAttractionIntegral coloumbIntegral(pBF.corePosition(), qBF.corePosition(), corePositionC,
                                                                  pP.exponent(), qP.exponent(),
                                                                  m_angularMomentumMax);
                result -= m_coreCharge * pP.weight() * qP.weight() * coloumbIntegral.coloumbAttractionIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                                                               qP.xExponent(), qP.yExponent(), qP.zExponent());
            }
        }
    }
    return result;
}

double GaussianSystem::overlapIntegral(int p, int q)
{
    double result = 0;
    const GaussianContractedOrbital& pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
            GaussianOverlapIntegral integrator(pBF.corePosition(), qBF.corePosition(),pP.exponent(), qP.exponent(), m_angularMomentumMax);
            result += pP.weight() * qP.weight() * integrator.overlapIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                             qP.xExponent(), qP.yExponent(), qP.zExponent());
        }
    }
    return result;
}

uint GaussianSystem::nBasisFunctions()
{
    return m_nBasisFunctions;
}

uint GaussianSystem::nParticles()
{
    return m_nParticles;
}

double GaussianSystem::additionalEnergyTerms()
{
    double result = 0;
    for(uint i = 0; i < m_corePositions.n_rows; i++) {
        for(uint j = i + 1; j < m_corePositions.n_rows; j++) {
            result += m_coreCharge*m_coreCharge / sqrt(dot(m_corePositions.row(i) - m_corePositions.row(j), m_corePositions.row(i)- m_corePositions.row(j)));
        }
    }
    return result;
}
