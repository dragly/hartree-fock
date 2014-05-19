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

double GaussianContractedOrbital::evaluated(const Vector3& position) const
{
    Vector3 diff = position - corePosition();
    double rSquared = dot(diff, diff);
    double result = 0;
    for(uint i = 0; i < m_primitiveBasisFunctions.size(); i++) {
        const GaussianPrimitiveOrbital &p = m_primitiveBasisFunctions.at(i);
        result += p.weight() * pow(diff.x(), p.xExponent()) * pow (diff.y(), p.yExponent()) * pow(diff.z(), p.zExponent())
                * exp(-p.exponent() * rSquared);
    }
    return result;
}


