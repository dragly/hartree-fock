#include "gaussiancore.h"

#include <parsers/turbomoleparser.h>

GaussianCore::GaussianCore(rowvec position, string fileName) :
    m_position(position)
{
    if(fileName != "") {
        load(fileName);
    }
}

void GaussianCore::load(string fileName)
{
    TurboMoleParser parser;
    parser.load(fileName);
    m_atomType = parser.atomType();
    m_nElectrons = int(m_atomType);
    m_charge = int(m_atomType);
    m_contractedOrbitals = parser.contractedBasisFunctions();
    for(GaussianContractedOrbital &contracted : m_contractedOrbitals) {
        contracted.setCorePosition(m_position);
    }
}

rowvec GaussianCore::position() const
{
    return m_position;
}

void GaussianCore::setPosition(const rowvec &position)
{
    m_position = position;
    for(GaussianContractedOrbital &contracted : m_contractedOrbitals) {
        contracted.setCorePosition(m_position);
    }
}

const vector<GaussianContractedOrbital> &GaussianCore::contractedOrbitals() const
{
    return m_contractedOrbitals;
}

void GaussianCore::setContractedOrbitals(const vector<GaussianContractedOrbital> &contractedOrbitals)
{
    m_contractedOrbitals = contractedOrbitals;
}

int GaussianCore::nElectrons() const
{
    return m_nElectrons;
}

void GaussianCore::setNElectrons(int nElectrons)
{
    m_nElectrons = nElectrons;
}
HF::AtomType GaussianCore::atomType() const
{
    return m_atomType;
}

void GaussianCore::setAtomType(const HF::AtomType &atomType)
{
    m_atomType = atomType;
}
int GaussianCore::charge() const
{
    return m_charge;
}

void GaussianCore::setCharge(int charge)
{
    m_charge = charge;
}





