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

private:
    double m_weight;
    double m_exponent;
};

#endif // GAUSSIANTYPEORBITAL_H
