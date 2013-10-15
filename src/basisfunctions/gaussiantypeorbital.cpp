#include "gaussiantypeorbital.h"

GaussianTypeOrbital::GaussianTypeOrbital(double weight, double exponent) :
    m_xExponent(0),
    m_yExponent(0),
    m_zExponent(0)
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
int GaussianTypeOrbital::xExponent() const
{
    return m_xExponent;
}

void GaussianTypeOrbital::setXExponent(int xExponent)
{
    m_xExponent = xExponent;
}
int GaussianTypeOrbital::yExponent() const
{
    return m_yExponent;
}

void GaussianTypeOrbital::setYExponent(int yExponent)
{
    m_yExponent = yExponent;
}
int GaussianTypeOrbital::zExponent() const
{
    return m_zExponent;
}

void GaussianTypeOrbital::setZExponent(int zExponent)
{
    m_zExponent = zExponent;
}



