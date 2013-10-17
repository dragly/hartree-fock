#include "atomicbasisfunction.h"

AtomicBasisFunction::AtomicBasisFunction() :
    m_weight(0),
    m_nucleusPosition(zeros(3)),
    m_xExponent(0),
    m_yExponent(0),
    m_zExponent(0)
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
int AtomicBasisFunction::xExponent() const
{
    return m_xExponent;
}

void AtomicBasisFunction::setXExponent(int xExponent)
{
    m_xExponent = xExponent;
}
int AtomicBasisFunction::yExponent() const
{
    return m_yExponent;
}

void AtomicBasisFunction::setYExponent(int yExponent)
{
    m_yExponent = yExponent;
}
int AtomicBasisFunction::zExponent() const
{
    return m_zExponent;
}

void AtomicBasisFunction::setZExponent(int zExponent)
{
    m_zExponent = zExponent;
}





