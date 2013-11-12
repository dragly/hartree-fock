#include "gaussiansystem.h"

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

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

void GaussianSystem::addContractedOrbital(const GaussianContractedOrbital &contractedOrbital)
{
    m_basisFunctions.push_back(contractedOrbital);
}

void GaussianSystem::addContractedOrbitals(const vector<GaussianContractedOrbital> &contractedOrbitals)
{
    m_basisFunctions.insert(m_basisFunctions.end(), contractedOrbitals.begin(), contractedOrbitals.end());
}

double GaussianSystem::particleDensity(const mat& C, double x, double y, double z) {
    if(C.n_rows < m_basisFunctions.size() || C.n_cols < m_nParticles / 2) {
        cout << "C matrix has the wrong dimensions" << endl;
        throw exception();
    }
//    ivec bfIndices = linspace<ivec>(0, m_nParticles / 2 - 1, m_nParticles / 2);
//    ivec bfIndices2 = linspace<ivec>(0, m_nParticles / 2 - 1, m_nParticles / 2);

//    mat overlapIntegrals = zeros(m_basisFunctions.size(), m_basisFunctions.size());
//    for(uint p = 0; p < overlapIntegrals.n_rows; p++) {
//        for(uint q = 0; q < overlapIntegrals.n_cols; q++) {
//            overlapIntegrals(p,q) = overlapIntegral(p, q);
//        }
//    }

//    vec evaluations = zeros(m_basisFunctions.size());
//    for(uint i = 0; i < evaluations.n_elem; i++) {
//        const GaussianContractedOrbital &bf = m_basisFunctions.at(i);
//        evaluations(i) = bf.evaluated(x,y,z);
//    }

//    cout << "nParticles: " << m_nParticles << endl;
    double result = 0;

    // TODO This should include some integrals over a lot of overlap functions - if not, why?s
    for(int i = 0; i < m_nParticles / 2; i++) {
        double innerResult = 0;
        for(int j = 0; j < m_basisFunctions.size(); j++) {
            const GaussianContractedOrbital &bf = m_basisFunctions.at(j);
            double evaluation = bf.evaluated(x,y,z);
            innerResult += C(j,i) * C(j,i) * evaluation * evaluation;
        }
        result += innerResult * innerResult;
    }

////    uint p1 = 0;
////    do {
////        uint p2 = 0;
////        do {
//    int sign1 = 1 - 2 * (p1 % 2);
//    int sign2 = 1 - 2 * (p2 % 2);

//    // Multiply the gaussians for x_0
//    double productResult = 0;
//    for(int j = 0; j < m_basisFunctions.size(); j++) {
//        productResult += C(j, bfIndices(0)) * evaluations(j);
//    }
//    for(int j = 0; j < m_basisFunctions.size(); j++) {
//        productResult += C(j, bfIndices2(0)) * evaluations(j);
//    }

//    // Do the overlap integrals for the rest of the gaussians
//    double overlapResult = 0;
//    for(int i = 1; i < m_nParticles / 2; i++) {
//        for(int j = 0; j < m_basisFunctions.size(); j++) {
//            for(int k = 0; k < m_basisFunctions.size(); k++) {
//                overlapResult += C(j,bfIndices(i)) * C(k, bfIndices2(i)) * overlapIntegrals(j,k);
//            }
//        }
//    }
//    result += sign1 * sign2 * overlapResult * productResult;
//    p2++;
//        } while (next_permutation(bfIndices2.memptr(), bfIndices2.memptr() + bfIndices2.n_elem));
//        p1++;
//    } while (next_permutation(bfIndices.memptr(), bfIndices.memptr() + bfIndices.n_elem));
    // Constructing the total wave function
    //    uint nBasis = m_basisFunctions.size();
    //    for (uint p = 0; p < nBasis; ++p) {
    //        int sign = 1 - 2 * (p % 2);
    //        cout << sign << " ";
    //        for(uint i = 0; i < nBasis; i++) {
    //            uint index = (i + p) % nBasis;
    //            cout << index;
    //        }
    //        cout << endl;
    //    }
    return result;
}
