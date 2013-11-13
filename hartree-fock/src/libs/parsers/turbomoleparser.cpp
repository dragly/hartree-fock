#include "turbomoleparser.h"

#include <fstream>
#include <iostream>
#include <boost/regex.hpp>

using namespace std;
using namespace boost;

TurboMoleParser::TurboMoleParser()
{
    // Orbitals expected for hydrogen
    nExpectedValenceSOrbitals.insert(pair<AtomType, int>(Hydrogen, 1));
    // Orbitals expected for nitrogen
    nExpectedCoreSOrbitals.insert(pair<AtomType, int>(Nitrogen, 1));
    nExpectedValenceSOrbitals.insert(pair<AtomType, int>(Nitrogen, 1));
    nExpectedValencePOrbitals.insert(pair<AtomType, int>(Nitrogen, 1));
    // Orbitals expected for oxygen
    nExpectedCoreSOrbitals.insert(pair<AtomType, int>(Oxygen, 1));
    nExpectedValenceSOrbitals.insert(pair<AtomType, int>(Oxygen, 1));
    nExpectedValencePOrbitals.insert(pair<AtomType, int>(Oxygen, 1));
    // Orbitals expected for sulfur
    nExpectedCoreSOrbitals.insert(pair<AtomType, int>(Sulfur, 2));
    nExpectedCorePOrbitals.insert(pair<AtomType, int>(Sulfur, 1));
    nExpectedValenceSOrbitals.insert(pair<AtomType, int>(Sulfur, 1));
    nExpectedValencePOrbitals.insert(pair<AtomType, int>(Sulfur, 1));
}

bool TurboMoleParser::read(string fileName)
{
    ifstream dataFile(fileName);
    string line;

    vector<int> nValenceOrbitals;
    string atomType;
    while (getline(dataFile, line))
    {
        bool skip = false;
        skip |= regex_match(line, regex("#.*"));
        skip |= regex_match(line, regex("$basis.*"));
        skip |= regex_match(line, regex("$end.*"));
        skip |= regex_match(line, regex("\\*.*"));
        if(skip) {
            continue;
        }
        smatch what;
        regex basisRegex("\\s*([a-zA-Z])\\s*([0-9])-([0-9]+)([G])\\s*");
        while(regex_search(line, what, basisRegex)) { // n 4-31G
            atomType = string(what[1]);
            int nCoreOrbitals(stoi(string(what[2])));
            for(char& c : string(what[3])) {
                int nValenceOrbital = c - '0';
                nValenceOrbitals.push_back(nValenceOrbital);
            }
            string basisSetType(what[4]);

            cout << "Atom: " << atomType << " nCore: " << nCoreOrbitals << " type: " << basisSetType << " nValence: ";
            for(int nValenceOrbital : nValenceOrbitals) {
                cout << nValenceOrbital << " ";
            }
            cout << endl;

            break;
        }
        regex orbitalRegex("\\s*([0-9])\\s*([spdf])\\s*");
        if(regex_search(line, what, orbitalRegex)) { // 4 s
            mergePrimitivesIntoContracted();
            if(what[2] == "s") {
                currentOrbitalType = sOrbitalType;
            } else if(what[2] == "p") {
                currentOrbitalType = pOrbitalType;
            } else if(what[2] == "d") {
                currentOrbitalType = dOrbitalType;
                cout << "ERROR: d orbitals not yet supported" << endl;
                return false;
            } else if(what[2] == "f") {
                currentOrbitalType = fOrbitalType;
                cout << "ERROR: f orbitals not yet supported" << endl;
                return false;
            } else {
                cout << "Unknown orbital type " << string(what[1]) << endl;
                return false;
            }
        }
        regex exponentWeightRegex("\\s*(-?[0-9]+\\.?[0-9]+)\\s*(-?[0-9]+\\.?[0-9]+)\\s*"); // 2.1 4.9
        if(regex_search(line, what, exponentWeightRegex)) {
            double exponent = stod(string(what[1]));
            double weight = stod(string(what[2]));
            GaussianPrimitiveOrbital primitive;
            primitive.setWeight(weight);
            primitive.setExponent(exponent);
            collectedPrimitiveBasisFunctions.push_back(primitive);
        }
    }
    mergePrimitivesIntoContracted();
    for(GaussianContractedOrbital orbital : contractedBasisFunctions) {
        cout << "Contracted:" << endl;
        for(GaussianPrimitiveOrbital primitive : orbital.primitiveBasisFunctions()) {
            cout << primitive.weight() << " * " << primitive.exponent() << endl;
        }
    }
    cout << "Done!" << endl;
    return true;
}

void TurboMoleParser::mergePrimitivesIntoContracted()
{
    if(collectedPrimitiveBasisFunctions.size() > 0) {
        switch(currentOrbitalType) {
        case sOrbitalType: {
            GaussianContractedOrbital contractedOrbital;
            for(GaussianPrimitiveOrbital primitive : collectedPrimitiveBasisFunctions) {
                double weightAdjusted = primitive.weight() * pow(2 * primitive.exponent() / M_PI, 3.0/4.0);
                primitive.setWeight(weightAdjusted);
                contractedOrbital.addPrimitiveBasisFunction(primitive);
            }
            contractedBasisFunctions.push_back(contractedOrbital);
            break;
        }
        case pOrbitalType: {
            int dim = 3;
            for(int i = 0; i < dim; i++) {
                GaussianContractedOrbital contractedOrbital;
                for(GaussianPrimitiveOrbital primitive : collectedPrimitiveBasisFunctions) {
                    double weightAdjusted = primitive.weight() * pow(2 * primitive.exponent() / M_PI, 3.0/4.0) * 2 * sqrt(primitive.exponent());
                    primitive.setWeight(weightAdjusted);
                    primitive.setXExponent(i == 0);
                    primitive.setYExponent(i == 1);
                    primitive.setZExponent(i == 2);
                    contractedOrbital.addPrimitiveBasisFunction(primitive);
                }
                contractedBasisFunctions.push_back(contractedOrbital);
            }
            break;
            // TODO Support d and f type orbitals
        }
        case dOrbitalType: {
            cout << "ERROR: d orbitals not yet supported" << endl;
            break;
        }
        case fOrbitalType: {
            cout << "ERROR: f orbitals not yet supported" << endl;
            break;
        }
        default: {
            cout << "ERROR: Unknown orbital type." << endl;
            break;
        }
        }
        collectedPrimitiveBasisFunctions.clear();
    }
}
