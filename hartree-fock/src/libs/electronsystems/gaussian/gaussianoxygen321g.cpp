#include "gaussianoxygen321g.h"

#include <basisfunctions/gaussian/integrals/gaussianoverlapintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiankineticintegral.h>
#include <basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h>
#include <basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h>

#include <iostream>

using namespace std;

GaussianOxygen321G::GaussianOxygen321G() :
    m_nParticles(2),
    m_nBasisFunctions(2 * 2),
    electronInteractionIntegral(0)
{
    rowvec corePositions = {-0.7, 0, 0,
                            0.7, 0, 0};
    m_corePositions = corePositions;
    m_corePositions.reshape(3,2);
    m_corePositions = m_corePositions.t();
    cout << m_corePositions << endl;
//    for(uint i = 0; i < m_corePositions.n_rows; i++) {
//        for(uint j = 0; j < m_nBasisFunctions; j++) {
//            GaussianContractedOrbital contracted(m_corePositions.row(i));
//            for(int k = 0; k < 4; k++) {
//                GaussianPrimitiveOrbital primitive(1,0,0,0,0.1 * (k + 1));
//                contracted.addPrimitiveBasisFunction(primitive);
//            }
//            m_basisFunctions.push_back(contracted);
//        }
//    }

    for(uint i = 0; i < m_corePositions.n_rows; i++) {
//        for(uint j = 0; j < m_nBasisFunctions; j++) {
//            GaussianContractedOrbital contracted(m_corePositions.row(i));
//            for(int k = 0; k < 4; k++) {
//                GaussianPrimitiveOrbital primitive(1,0,0,0,0.1 * (k + 1));
//                contracted.addPrimitiveBasisFunction(primitive);
//            }
//            m_basisFunctions.push_back(contracted);
//        }

//        double factorAlwaysNeeded = pow(2 * 5.4471780 / M_PI, 0.75);
        GaussianContractedOrbital contracted1(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive11(0.1562850 * pow(2 * 5.4471780 / M_PI, 0.75), 0, 0, 0, 5.4471780);
        GaussianPrimitiveOrbital primitive12(0.9046910 * pow(2 * 0.8245470 / M_PI, 0.75), 0, 0, 0, 0.8245470);
        contracted1.addPrimitiveBasisFunction(primitive11);
        contracted1.addPrimitiveBasisFunction(primitive12);

        GaussianContractedOrbital contracted2(m_corePositions.row(i));
        GaussianPrimitiveOrbital primitive21(1.0000000 * pow(2 * 0.1831920 / M_PI, 0.75), 0, 0, 0, 0.1831920);
        contracted2.addPrimitiveBasisFunction(primitive21);
        m_basisFunctions.push_back(contracted1);
        m_basisFunctions.push_back(contracted2);
    }
}


double GaussianOxygen321G::coupledIntegral(int p, int r, int q, int s)
{
    double result = 0;
    const GaussianContractedOrbital& pBasisFunction = m_basisFunctions.at(p);
    const GaussianContractedOrbital& rBasisFunction = m_basisFunctions.at(r);
    const GaussianContractedOrbital& qBasisFunction = m_basisFunctions.at(q);
    const GaussianContractedOrbital& sBasisFunction = m_basisFunctions.at(s);
    for(const GaussianPrimitiveOrbital& pPrimitive : pBasisFunction.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& rPrimitive : rBasisFunction.primitiveBasisFunctions()) {
            electronInteractionIntegral.setAB(pBasisFunction.corePosition(), rBasisFunction.corePosition(),
                                              pPrimitive.exponent(), rPrimitive.exponent());
            for(const GaussianPrimitiveOrbital& qPrimitive : qBasisFunction.primitiveBasisFunctions()) {
                for(const GaussianPrimitiveOrbital& sPrimitive : sBasisFunction.primitiveBasisFunctions()) {
                    electronInteractionIntegral.setCD(qBasisFunction.corePosition(), sBasisFunction.corePosition(),
                                                      qPrimitive.exponent(), sPrimitive.exponent());
                    result += electronInteractionIntegral.electronInteractionIntegral(pPrimitive.xExponent(), pPrimitive.yExponent(), pPrimitive.zExponent(),
                                                                                      rPrimitive.xExponent(), rPrimitive.yExponent(), rPrimitive.zExponent(),
                                                                                      qPrimitive.xExponent(), qPrimitive.yExponent(), qPrimitive.zExponent(),
                                                                                      sPrimitive.xExponent(), sPrimitive.yExponent(), sPrimitive.zExponent());
                }
            }
        }
    }

    return result;
}

double GaussianOxygen321G::uncoupledIntegral(int p, int q)
{
    double result = 0;
    const GaussianContractedOrbital& pBasisFunction = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBasisFunction = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pPrimitive : pBasisFunction.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qPrimitive : qBasisFunction.primitiveBasisFunctions()) {
            GaussianKineticIntegral kineticIntegral(pBasisFunction.corePosition(), qBasisFunction.corePosition(),
                                                    pPrimitive.exponent(), qPrimitive.exponent(),
                                                    0);
            result += kineticIntegral.kineticIntegral(pPrimitive.xExponent(), pPrimitive.yExponent(), pPrimitive.zExponent(),
                                                      qPrimitive.xExponent(), qPrimitive.yExponent(), qPrimitive.zExponent());
            for(uint i = 0; i < m_corePositions.n_rows; i++) {
                const rowvec &corePositionC = m_corePositions.row(i);
                GaussianColoumbAttractionIntegral coloumbIntegral(pBasisFunction.corePosition(), qBasisFunction.corePosition(), corePositionC,
                                                                  pPrimitive.exponent(), qPrimitive.exponent(),
                                                                  0);
                result += coloumbIntegral.coloumbAttractionIntegral(pPrimitive.xExponent(), pPrimitive.yExponent(), pPrimitive.zExponent(),
                                                                    qPrimitive.xExponent(), qPrimitive.yExponent(), qPrimitive.zExponent());
            }
        }
    }
    return result;
}

double GaussianOxygen321G::overlapIntegral(int p, int q)
{
    double result = 0;
    const GaussianContractedOrbital& pBasisFunction = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBasisFunction = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pPrimitive : pBasisFunction.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qPrimitive : qBasisFunction.primitiveBasisFunctions()) {
            GaussianOverlapIntegral primitiveOverlapIntegral(pBasisFunction.corePosition(), qBasisFunction.corePosition(),
                                                             pPrimitive.exponent(), qPrimitive.exponent(),
                                                             0);
            result += primitiveOverlapIntegral.overlapIntegral(pPrimitive.xExponent(), pPrimitive.yExponent(), pPrimitive.zExponent(),
                                                               qPrimitive.xExponent(), qPrimitive.yExponent(), qPrimitive.zExponent());
        }
    }
    return result;
}

uint GaussianOxygen321G::nBasisFunctions()
{
    return m_nBasisFunctions;
}

uint GaussianOxygen321G::nParticles()
{
    return m_nParticles;
}

double GaussianOxygen321G::additionalEnergyTerms()
{
    double result = 0;
    for(uint i = 0; i < m_corePositions.n_rows; i++) {
        for(uint j = i + 1; j < m_corePositions.n_rows; j++) {
            result += 1/sqrt(dot(m_corePositions.row(i) - m_corePositions.row(j), m_corePositions.row(i)- m_corePositions.row(j)));
        }
    }
    return result;
}
