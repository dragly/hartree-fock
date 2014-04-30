#ifndef GAUSSIANCORE_H
#define GAUSSIANCORE_H

#include <armadillo>
#include <string>
#include <hf.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>

class Vector3;

using namespace arma;
using namespace std;

class GaussianCore
{
public:
    GaussianCore(Vector3 position = Vector3::createZeros(), string fileName = "");
    GaussianCore(Vector3 position, int atomNumber, string basisName);
    GaussianCore(Vector3 position, string atomAbbreviation, string basisName);

    void load(string fileName);

    const Vector3 &position() const;
    void setPosition(const Vector3 &position);

    const vector<GaussianContractedOrbital> &contractedOrbitals() const;
    void setContractedOrbitals(const vector<GaussianContractedOrbital> &contractedOrbitals);

    int nElectrons() const;
    void setNElectrons(int nElectrons);

    HF::AtomType atomType() const;
    void setAtomType(const HF::AtomType &atomType);

    int charge() const;
    void setCharge(int charge);

private:
    Vector3 m_position;
    int m_charge;
    int m_nElectrons;
    HF::AtomType m_atomType;
    vector<GaussianContractedOrbital> m_contractedOrbitals;
};

#endif // GAUSSIANCORE_H
