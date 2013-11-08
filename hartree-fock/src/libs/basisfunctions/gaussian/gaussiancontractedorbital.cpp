#include "gaussiancontractedorbital.h"

//GaussianContractedOrbital::GaussianContractedOrbital()
//{
//}

GaussianContractedOrbital::GaussianContractedOrbital(const rowvec& corePosition) :
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
rowvec GaussianContractedOrbital::corePosition() const
{
    return m_corePosition;
}

void GaussianContractedOrbital::setCorePosition(const rowvec &corePosition)
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


