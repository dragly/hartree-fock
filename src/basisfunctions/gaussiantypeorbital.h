#ifndef GAUSSIANTYPEORBITAL_H
#define GAUSSIANTYPEORBITAL_H

#include <armadillo>

using namespace arma;

class GaussianTypeOrbital
{
public:
    GaussianTypeOrbital(double weight = 0, double exponent = 0);

    double exponent() const;
    void setExponent(double exponent);

    double weight() const;
    void setWeight(double weight);

    int xExponent() const;
    void setXExponent(int xExponent);

    int yExponent() const;
    void setYExponent(int yExponent);

    int zExponent() const;
    void setZExponent(int zExponent);

private:
    double m_weight;
    double m_exponent;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
};

#endif // GAUSSIANTYPEORBITAL_H
