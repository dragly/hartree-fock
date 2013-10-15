#include "atomicbasisfunction.h"

AtomicBasisFunction::AtomicBasisFunction() :
    m_weight(0),
    m_nucleusPosition(zeros(3))
{
}
const vector<GaussianTypeOrbital> &AtomicBasisFunction::orbitals() const
{
    return m_orbitals;
}

void AtomicBasisFunction::setOrbitals(const vector<GaussianTypeOrbital> &orbitals)
{
    m_orbitals = orbitals;
}
double AtomicBasisFunction::weight() const
{
    return m_weight;
}

void AtomicBasisFunction::setWeight(double weight)
{
    m_weight = weight;
}
vec AtomicBasisFunction::nucleusPosition() const
{
    return m_nucleusPosition;
}

void AtomicBasisFunction::setNucleusPosition(const vec &nucleusPosition)
{
    m_nucleusPosition = nucleusPosition;
}



