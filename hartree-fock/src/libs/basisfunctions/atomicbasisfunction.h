#ifndef ATOMICBASISFUNCTION_H
#define ATOMICBASISFUNCTION_H

#include <vector>
#include <basisfunctions/gaussiantypeorbital.h>

using namespace std;

class AtomicBasisFunction
{
public:
    AtomicBasisFunction();

    const vector<GaussianTypeOrbital> &orbitals() const;
    void setOrbitals(const vector<GaussianTypeOrbital> &orbitals);

    double weight() const;
    void setWeight(double weight);

    vec nucleusPosition() const;
    void setNucleusPosition(const vec &nucleusPosition);

    int xExponent() const;
    void setXExponent(int xExponent);

    int yExponent() const;
    void setYExponent(int yExponent);

    int zExponent() const;
    void setZExponent(int zExponent);

private:
    vector<GaussianTypeOrbital> m_orbitals;
    double m_weight;
    vec m_nucleusPosition;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
};

#endif // ATOMICBASISFUNCTION_H
