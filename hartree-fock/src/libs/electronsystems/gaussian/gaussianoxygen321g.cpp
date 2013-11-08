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


double GaussianOxygen321G::coupledIntegral(int p, int q, int r, int s)
{
    double result = 0;
    const GaussianContractedOrbital& pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& rBF = m_basisFunctions.at(r);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    const GaussianContractedOrbital& sBF = m_basisFunctions.at(s);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& rP : rBF.primitiveBasisFunctions()) {
            electronInteractionIntegral.setAC(pBF.corePosition(), rBF.corePosition(), pP.exponent(), rP.exponent());
            for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
                for(const GaussianPrimitiveOrbital& sP : sBF.primitiveBasisFunctions()) {
                    electronInteractionIntegral.setBD(qBF.corePosition(), sBF.corePosition(), qP.exponent(), sP.exponent());
                    result += pP.weight() * rP.weight() * qP.weight() * sP.weight()
                            * electronInteractionIntegral.electronInteractionIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                                      rP.xExponent(), rP.yExponent(), rP.zExponent(),
                                                                                      qP.xExponent(), qP.yExponent(), qP.zExponent(),
                                                                                      sP.xExponent(), sP.yExponent(), sP.zExponent());
                }
            }
        }
    }

    return result;
}

double GaussianOxygen321G::uncoupledIntegral(int p, int q)
{
    double result = 0;
    const GaussianContractedOrbital& pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
            GaussianKineticIntegral kineticIntegral(pBF.corePosition(), qBF.corePosition(),
                                                    pP.exponent(), qP.exponent(),
                                                    0);
            result += pP.weight() * qP.weight() * kineticIntegral.kineticIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                      qP.xExponent(), qP.yExponent(), qP.zExponent());
            for(uint i = 0; i < m_corePositions.n_rows; i++) {
                const rowvec &corePositionC = m_corePositions.row(i);
                GaussianColoumbAttractionIntegral coloumbIntegral(pBF.corePosition(), qBF.corePosition(), corePositionC,
                                                                  pP.exponent(), qP.exponent(),
                                                                  0);
                result -= pP.weight() * qP.weight() * coloumbIntegral.coloumbAttractionIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                    qP.xExponent(), qP.yExponent(), qP.zExponent());
            }
        }
    }
    return result;
}

double GaussianOxygen321G::overlapIntegral(int p, int q)
{
    double result = 0;
    const GaussianContractedOrbital& pBF = m_basisFunctions.at(p);
    const GaussianContractedOrbital& qBF = m_basisFunctions.at(q);
    for(const GaussianPrimitiveOrbital& pP : pBF.primitiveBasisFunctions()) {
        for(const GaussianPrimitiveOrbital& qP : qBF.primitiveBasisFunctions()) {
            GaussianOverlapIntegral integrator(pBF.corePosition(), qBF.corePosition(),pP.exponent(), qP.exponent(), 0);
            result += pP.weight() * qP.weight() * integrator.overlapIntegral(pP.xExponent(), pP.yExponent(), pP.zExponent(),
                                                                             qP.xExponent(), qP.yExponent(), qP.zExponent());
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
