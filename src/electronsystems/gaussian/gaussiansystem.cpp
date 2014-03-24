#include "gaussiansystem.h"

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

#include <iostream>
#include <cmath>

using namespace std;

GaussianSystem::GaussianSystem() :
    m_nParticles(0),
    m_nBasisFunctions(0),
    m_angularMomentumMax(0),
    electronInteractionIntegral(0),
    kineticIntegral(0),
    coloumbIntegral(0)
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
            kineticIntegral.set(pBF.corePosition(), qBF.corePosition(), pP.exponent(), qP.exponent());
            result += pP.weight() * qP.weight() * kineticIntegral.kineticIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                                  qP.xExponent(), qP.yExponent(), qP.zExponent());
            for(uint i = 0; i < m_cores.size(); i++) {
                const GaussianCore &core = m_cores.at(i);
                const rowvec &corePositionC = core.position();
                coloumbIntegral.set(pBF.corePosition(), qBF.corePosition(), corePositionC,
                                    pP.exponent(), qP.exponent());
                result -= core.charge() * pP.weight() * qP.weight() * coloumbIntegral.coloumbAttractionIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
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
    for(uint i = 0; i < m_cores.size(); i++) {
        const GaussianCore &core1 = m_cores.at(i);
        for(uint j = i + 1; j < m_cores.size(); j++) {
            const GaussianCore &core2 = m_cores.at(j);
            result += core1.charge()*core2.charge() / sqrt(dot(core1.position() - core2.position(), core1.position() - core2.position()));
        }
    }
    return result;
}

double GaussianSystem::particleDensity(const mat& C, double x, double y, double z) const {
    if(C.n_rows < m_basisFunctions.size() || C.n_cols < m_nParticles / 2) {
        cout << "C matrix has the wrong dimensions" << endl;
        throw exception();
    }
    double result = 0;

    int nk = m_nParticles / 2;
    nk = max(nk, 1);
    for(int i = 0; i < nk; i++) {
        double innerResult = 0;
        for(int j = 0; j < int(m_basisFunctions.size()); j++) {
            const GaussianContractedOrbital &bf = m_basisFunctions.at(j);
            double evaluation = bf.evaluated(x,y,z);
            innerResult += C(j,i) * C(j,i) * evaluation * evaluation;
        }
        result += innerResult * innerResult;
    }
    return result;
}

void GaussianSystem::addCore(const GaussianCore &core)
{
    m_nBasisFunctions += core.contractedOrbitals().size();
    m_nParticles += core.nElectrons();
    m_basisFunctions.insert(m_basisFunctions.end(), core.contractedOrbitals().begin(), core.contractedOrbitals().end());
    m_cores.push_back(core);
    int angularMomentumMax = 0;
    for(const GaussianContractedOrbital &contracted : m_basisFunctions) {
        for(const GaussianPrimitiveOrbital &primitive : contracted.primitiveBasisFunctions()) {
            angularMomentumMax = max(angularMomentumMax, primitive.xExponent());
            angularMomentumMax = max(angularMomentumMax, primitive.yExponent());
            angularMomentumMax = max(angularMomentumMax, primitive.zExponent());
        }
    }
    setAngularMomentumMax(angularMomentumMax);
}

void GaussianSystem::setAngularMomentumMax(int angularMomentumMax)
{
    m_angularMomentumMax = angularMomentumMax;
    electronInteractionIntegral = GaussianElectronInteractionIntegral(angularMomentumMax);
    kineticIntegral = GaussianKineticIntegral(angularMomentumMax);
    coloumbIntegral = GaussianColoumbAttractionIntegral(angularMomentumMax);
}
