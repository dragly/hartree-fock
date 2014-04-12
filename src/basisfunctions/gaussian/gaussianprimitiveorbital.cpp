#include "gaussianprimitiveorbital.h"

#include <algorithm>

using std::max;

GaussianPrimitiveOrbital::GaussianPrimitiveOrbital() :
    m_weight(0),
    m_xExponent(0),
    m_yExponent(0),
    m_zExponent(0),
    m_exponent(0)
{
}

GaussianPrimitiveOrbital::GaussianPrimitiveOrbital(double weight, int xExponent, int yExponent, int zExponent, double exponent) :
    m_weight(weight),
    m_xExponent(xExponent),
    m_yExponent(yExponent),
    m_zExponent(zExponent),
    m_exponent(exponent)
{

}
double GaussianPrimitiveOrbital::exponent() const
{
    return m_exponent;
}

void GaussianPrimitiveOrbital::setExponent(double exponent)
{
    m_exponent = exponent;
}
int GaussianPrimitiveOrbital::zExponent() const
{
    return m_zExponent;
}

void GaussianPrimitiveOrbital::setZExponent(int zExponent)
{
    m_zExponent = zExponent;
}
int GaussianPrimitiveOrbital::yExponent() const
{
    return m_yExponent;
}

void GaussianPrimitiveOrbital::setYExponent(int yExponent)
{
    m_yExponent = yExponent;
}
int GaussianPrimitiveOrbital::xExponent() const
{
    return m_xExponent;
}

void GaussianPrimitiveOrbital::setXExponent(int xExponent)
{
    m_xExponent = xExponent;
}
double GaussianPrimitiveOrbital::weight() const
{
    return m_weight;
}

void GaussianPrimitiveOrbital::setWeight(double weight)
{
    m_weight = weight;
}

int GaussianPrimitiveOrbital::angularMomentum() const
{
    return m_xExponent + m_yExponent + m_zExponent;
}

int GaussianPrimitiveOrbital::exponentMax() const
{
    return max(m_xExponent, max(m_yExponent, m_zExponent));
}

