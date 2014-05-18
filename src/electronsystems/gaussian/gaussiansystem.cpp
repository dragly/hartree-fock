#include "gaussiansystem.h"

#include "math/vector3.h"
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
    coulombIntegral(0)
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
            electronInteractionIntegral.setAB(pBF.corePosition(), qBF.corePosition(), pP, qP);
            for(const GaussianPrimitiveOrbital& rP : rBF.primitiveBasisFunctions()) {
                for(const GaussianPrimitiveOrbital& sP : sBF.primitiveBasisFunctions()) {
                    electronInteractionIntegral.setCD(rBF.corePosition(), sBF.corePosition(),
                                                      rP, sP);
                    result += pP.weight() * rP.weight() * qP.weight() * sP.weight()
                            * electronInteractionIntegral.electronInteractionIntegral(pP, qP, rP, sP);
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
            kineticIntegral.set(pBF.corePosition(), qBF.corePosition(), pP, qP);
            result += pP.weight() * qP.weight() * kineticIntegral.kineticIntegral(pP, qP);
            for(uint i = 0; i < m_cores.size(); i++) {
                const GaussianCore &core = m_cores.at(i);
                const Vector3 &corePositionC = core.position();
                coulombIntegral.set(pBF.corePosition(), qBF.corePosition(), corePositionC,
                                    pP, qP);
                result -= core.charge() * pP.weight() * qP.weight() * coulombIntegral.coloumbAttractionIntegral(pP, qP);
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
            GaussianOverlapIntegral integrator(pBF.corePosition(), qBF.corePosition(), pP, qP);
            result += pP.weight() * qP.weight() * integrator.overlapIntegral(pP, qP);
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

double GaussianSystem::orbitalDensity(uint orbital, const mat& C, double x, double y, double z) const
{
    double innerResult = 0;
    for(int j = 0; j < int(m_basisFunctions.size()); j++) {
        const GaussianContractedOrbital &bf = m_basisFunctions.at(j);
        double evaluation = bf.evaluated(x,y,z);
        innerResult += C(j,orbital) * C(j,orbital) * evaluation * evaluation;
    }
    return innerResult * innerResult;
}

double GaussianSystem::electronDensity(const mat& C, double x, double y, double z) const {
    if(C.n_rows < m_basisFunctions.size() || C.n_cols < m_nParticles / 2) {
        cout << "C matrix has the wrong dimensions" << endl;
        throw exception();
    }
    double result = 0;

    int nk = C.n_cols;
    for(int i = 0; i < nk; i++) {
        result += orbitalDensity(i, C, x, y, z);
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
    coulombIntegral = GaussianColoumbAttractionIntegral(angularMomentumMax);
}
