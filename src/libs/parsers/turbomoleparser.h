#ifndef TURBOMOLEPARSER_H
#define TURBOMOLEPARSER_H

#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>
#include <hf.h>

#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class TurboMoleParser
{
public:
    TurboMoleParser();

    bool load(string fileName);

    const vector<GaussianContractedOrbital> &contractedBasisFunctions() const;

    HF::AtomType atomType() const;

private:
    vector<GaussianContractedOrbital> m_contractedBasisFunctions;
    HF::AtomOrbitalType m_currentOrbitalType = HF::sOrbitalType;
    vector<GaussianPrimitiveOrbital> m_collectedPrimitiveBasisFunctions;
    void mergePrimitivesIntoContracted();
    HF::AtomType m_atomType;
};

#endif // TURBOMOLEPARSER_H
