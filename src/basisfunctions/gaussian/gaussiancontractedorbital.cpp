#include "gaussiancontractedorbital.h"

#include "math/vector3.h"

GaussianContractedOrbital::GaussianContractedOrbital() :
    GaussianContractedOrbital(Vector3::createZeros())
{
}

GaussianContractedOrbital::GaussianContractedOrbital(const Vector3 &corePosition) :
    m_corePosition(corePosition)
{

}

void GaussianContractedOrbital::addPrimitiveBasisFunction(const GaussianPrimitiveOrbital &primitive)
{
    m_primitiveBasisFunctions.push_back(primitive);
}

const GaussianPrimitiveOrbital &GaussianContractedOrbital::primitiveBasisFunction(int index)
{
    return m_primitiveBasisFunctions.at(index);
}
Vector3 GaussianContractedOrbital::corePosition() const
{
    return m_corePosition;
}

void GaussianContractedOrbital::setCorePosition(const Vector3 &corePosition)
{
    m_corePosition = corePosition;
}
const vector<GaussianPrimitiveOrbital> &GaussianContractedOrbital::primitiveBasisFunctions() const
{
    return m_primitiveBasisFunctions;
}

void GaussianContractedOrbital::setPrimitiveBasisFunctions(const vector<GaussianPrimitiveOrbital> &primitiveBasisFunctions)
{
    m_primitiveBasisFunctions = primitiveBasisFunctions;
}

double GaussianContractedOrbital::evaluated(double x, double y, double z) const
{
    double xDiff = x - corePosition()(0);
    double yDiff = y - corePosition()(1);
    double zDiff = z - corePosition()(2);
    double rSquared = (xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);
    double result = 0;
    for(uint i = 0; i < m_primitiveBasisFunctions.size(); i++) {
        const GaussianPrimitiveOrbital &p = m_primitiveBasisFunctions.at(i);
        result += p.weight() * pow(xDiff, p.xExponent()) * pow (yDiff, p.yExponent()) * pow(zDiff, p.zExponent())
                * exp(-p.exponent() * rSquared);
    }
    return result;
}


