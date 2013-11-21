#ifndef GAUSSIANPRIMITIVEORBITAL_H
#define GAUSSIANPRIMITIVEORBITAL_H

class GaussianPrimitiveOrbital
{
public:
    explicit GaussianPrimitiveOrbital();
    explicit GaussianPrimitiveOrbital(double weight, int xExponent, int yExponent, int zExponent, double exponent);

    double exponent() const;
    void setExponent(double exponent);

    int zExponent() const;
    void setZExponent(int zExponent);

    int yExponent() const;
    void setYExponent(int yExponent);

    int xExponent() const;
    void setXExponent(int xExponent);

    double weight() const;
    void setWeight(double weight);

private:
    double m_weight;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
    double m_exponent;
};

#endif // GAUSSIANPRIMITIVEORBITAL_H
