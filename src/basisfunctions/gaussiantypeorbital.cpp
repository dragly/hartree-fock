#include "gaussiantypeorbital.h"

GaussianTypeOrbital::GaussianTypeOrbital(double weight, double exponent)
{
}

double GaussianTypeOrbital::exponent() const
{
    return m_exponent;
}

void GaussianTypeOrbital::setExponent(double exponent)
{
    m_exponent = exponent;
}

double GaussianTypeOrbital::weight() const
{
    return m_weight;
}

void GaussianTypeOrbital::setWeight(double weight)
{
    m_weight = weight;
}

