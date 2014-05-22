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
    m_overlapIntegral(0),
    m_kineticIntegral(0),
    m_coulombIntegral(0),
    m_electronInteractionIntegral(0)
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
            m_electronInteractionIntegral.setAB(pBF.corePosition(), qBF.corePosition(), pP, qP);
            for(const GaussianPrimitiveOrbital& rP : rBF.primitiveBasisFunctions()) {
                for(const GaussianPrimitiveOrbital& sP : sBF.primitiveBasisFunctions()) {
                    m_electronInteractionIntegral.setCD(rBF.corePosition(), sBF.corePosition(),
                                                      rP, sP);
                    result += pP.weight() * rP.weight() * qP.weight() * sP.weight()
                            * m_electronInteractionIntegral.electronInteractionIntegral(pP, qP, rP, sP);
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
            m_kineticIntegral.set(pBF.corePosition(), qBF.corePosition(), pP, qP);
            result += pP.weight() * qP.weight() * m_kineticIntegral.kineticIntegral(pP, qP);
            for(uint i = 0; i < m_cores.size(); i++) {
                const GaussianCore &core = m_cores.at(i);
                const Vector3 &corePositionC = core.position();
                m_coulombIntegral.set(pBF.corePosition(), qBF.corePosition(), corePositionC,
                                    pP, qP);
                result -= core.charge() * pP.weight() * qP.weight() * m_coulombIntegral.coloumbAttractionIntegral(pP, qP);
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
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) { // TODO Optimize for symmetry if pBF = qBF
            m_overlapIntegral.set(pBF.corePosition(), qBF.corePosition(), pP, qP);
            result += pP.weight() * qP.weight() * m_overlapIntegral.overlapIntegral(pP, qP);
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

double GaussianSystem::electrostaticPotential(const mat& C, const Vector3 &position)
{
    double result = 0;
    result += corePotential(position);
    result -= electronPotential(C, position);
    return result;
}

double GaussianSystem::corePotential(const Vector3 &position)
{
    double potential = 0;
    for(const GaussianCore &core : m_cores) {
        Vector3 diff = position - core.position();
        double distance = sqrt(dot(diff, diff));
        potential += core.charge() / distance;
    }
    return potential;
}

double GaussianSystem::electronPotential(const mat& C, const Vector3 position)
{
    double result = 0;
    for(uint p = 0; p < (m_basisFunctions.size()); p++) {
        double potential = electronPotential(p,p,position);
        for(uint orbital = 0; orbital < C.n_cols; orbital++) {
            result += C(p, orbital) * C(p,orbital) * potential;
        }
        for(uint q = p+1; q < (m_basisFunctions.size()); q++) {
            potential = electronPotential(p,q,position);
            for(uint orbital = 0; orbital < C.n_cols; orbital++) {
                result += 2.0 * C(p, orbital) * C(q, orbital) * potential;
            }
        }
    }
    return result;
}

double GaussianSystem::electronPotential(uint p, uint q, const Vector3 &position)
{
    double result = 0;
    const GaussianContractedOrbital &pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
            const Vector3 &corePositionC = position;
            m_coulombIntegral.set(pBF.corePosition(), qBF.corePosition(), corePositionC, pP, qP);
            result += pP.weight() * qP.weight() * m_coulombIntegral.coloumbAttractionIntegral(pP, qP);
        }
    }
    return result;
}

//double GaussianSystem::orbitalDensity(uint orbital, const mat& C, double x, double y, double z) const
//{
//    double innerProduct = 0;
//    for(uint p = 0; p < (m_basisFunctions.size()); p++) {
//        const GaussianContractedOrbital &bfp = m_basisFunctions.at(p);
//        double evaluationp = bfp.evaluated(x,y,z);
//        double Cp = C(p,orbital);
//        innerProduct += C(p,orbital) * C(p,orbital) * evaluationp * evaluationp;
//        for(uint q = p+1; q < (m_basisFunctions.size()); q++) {
//            const GaussianContractedOrbital &bfq = m_basisFunctions.at(q);
//            double evaluationq = bfq.evaluated(x,y,z);
//            innerProduct += 2.0 * Cp * C(q,orbital) * evaluationp * evaluationq;
//        }
//    }
//    return innerProduct * innerProduct;
//}

rowvec GaussianSystem::orbitalDensities(const mat& C, const Vector3 &position) const
{
    rowvec result = zeros(C.n_cols);
    for(uint p = 0; p < (m_basisFunctions.size()); p++) {
        const GaussianContractedOrbital &bfp = m_basisFunctions.at(p);
        double evaluationp = bfp.evaluated(position);
        for(uint orbital = 0; orbital < C.n_cols; orbital++) {
            result(orbital) += C(p,orbital) * C(p,orbital) * evaluationp * evaluationp;
        }
        for(uint q = p+1; q < (m_basisFunctions.size()); q++) {
            const GaussianContractedOrbital &bfq = m_basisFunctions.at(q);
            double evaluationq = bfq.evaluated(position);
            for(uint orbital = 0; orbital < C.n_cols; orbital++) {
                result(orbital) += 2.0 * C(p,orbital) * C(q,orbital) * evaluationp * evaluationq;
            }
        }
    }
    return result;
}

double GaussianSystem::electronDensity(const mat& C, const Vector3 &position) const
{
    double result = sum(orbitalDensities(C, position));
//    for(uint p = 0; p < (m_basisFunctions.size()); p++) {
//        const GaussianContractedOrbital &bfp = m_basisFunctions.at(p);
//        double evaluationp = bfp.evaluated(position);
//        for(uint orbital = 0; orbital < C.n_cols; orbital++) {
//            result += C(p,orbital) * C(p,orbital) * evaluationp * evaluationp;
//        }
//        for(uint q = p+1; q < (m_basisFunctions.size()); q++) {
//            const GaussianContractedOrbital &bfq = m_basisFunctions.at(q);
//            double evaluationq = bfq.evaluated(position);
//            for(uint orbital = 0; orbital < C.n_cols; orbital++) {
//                result += 2.0 * C(p,orbital) * C(q,orbital) * evaluationp * evaluationq;
//            }
//        }
//    }
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
    m_overlapIntegral = GaussianOverlapIntegral(angularMomentumMax);
    m_electronInteractionIntegral = GaussianElectronInteractionIntegral(angularMomentumMax);
    m_kineticIntegral = GaussianKineticIntegral(angularMomentumMax);
    m_coulombIntegral = GaussianColoumbAttractionIntegral(angularMomentumMax);
}
