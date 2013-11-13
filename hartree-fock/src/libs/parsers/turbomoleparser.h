#ifndef TURBOMOLEPARSER_H
#define TURBOMOLEPARSER_H

#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>
#include <basisfunctions/gaussian/gaussiancontractedorbital.h>

#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class TurboMoleParser
{
public:
    TurboMoleParser();

    enum AtomOrbitalType {
        sOrbitalType,
        pOrbitalType,
        dOrbitalType,
        fOrbitalType
    };

    enum AtomType {
        Hydrogen,
        Nitrogen,
        Oxygen,
        Sulfur
    };

    bool read(string fileName);
private:
    unordered_map<AtomType, int, hash<int>> nExpectedCoreSOrbitals;
    unordered_map<AtomType, int, hash<int>> nExpectedCorePOrbitals;
    unordered_map<AtomType, int, hash<int>> nExpectedValenceSOrbitals;
    unordered_map<AtomType, int, hash<int>> nExpectedValencePOrbitals;
    vector<GaussianContractedOrbital> contractedBasisFunctions;
    AtomOrbitalType currentOrbitalType = sOrbitalType;
    vector<GaussianPrimitiveOrbital> collectedPrimitiveBasisFunctions;
    void mergePrimitivesIntoContracted();
};

#endif // TURBOMOLEPARSER_H
