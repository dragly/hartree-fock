#include "gaussiancore.h"

#include "math/vector3.h"
#include <parsers/turbomoleparser.h>

GaussianCore::GaussianCore(Vector3 position, string fileName) :
    m_position(position)
{
    if(fileName != "") {
        load(fileName);
    }
}

void GaussianCore::load(string fileName)
{
    TurboMoleParser parser;
    bool ok = parser.load(fileName);
    if(!ok) {
        cerr << "Could not load file " << fileName << endl;
        throw logic_error("Could not read basis file");
    }
    m_atomType = parser.atomType();
    m_nElectrons = int(m_atomType);
    m_charge = int(m_atomType);
    m_contractedOrbitals = parser.contractedBasisFunctions();
    for(GaussianContractedOrbital &contracted : m_contractedOrbitals) {
        contracted.setCorePosition(m_position);
    }
}

const Vector3& GaussianCore::position() const
{
    return m_position;
}

void GaussianCore::setPosition(const Vector3 &position)
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





