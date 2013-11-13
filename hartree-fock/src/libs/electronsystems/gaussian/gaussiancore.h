#ifndef GAUSSIANCORE_H
#define GAUSSIANCORE_H

#include <armadillo>
#include <string>
#include <hf.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>

using namespace arma;
using namespace std;

class GaussianCore
{
public:
    GaussianCore(rowvec position = zeros(3), string fileName = "");

    void load(string fileName);

    rowvec position() const;
    void setPosition(const rowvec &position);

    const vector<GaussianContractedOrbital> &contractedOrbitals() const;
    void setContractedOrbitals(const vector<GaussianContractedOrbital> &contractedOrbitals);

    int nElectrons() const;
    void setNElectrons(int nElectrons);

    HF::AtomType atomType() const;
    void setAtomType(const HF::AtomType &atomType);

    int charge() const;
    void setCharge(int charge);

private:
    rowvec m_position;
    int m_charge;
    int m_nElectrons;
    HF::AtomType m_atomType;
    vector<GaussianContractedOrbital> m_contractedOrbitals;
};

#endif // GAUSSIANCORE_H
