#include "turbomoleparser.h"

#include <fstream>
#include <iostream>
#include <boost/regex.hpp>
#include <clocale>

using namespace std;
using namespace boost;

TurboMoleParser::TurboMoleParser()
{
}

bool TurboMoleParser::load(string fileName)
{
    setlocale(LC_ALL, "C");
    m_contractedBasisFunctions.clear();
    ifstream dataFile(fileName);
    string line;

    vector<int> nValenceOrbitals;
    string atomTypeAbbreviation;
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
        regex basisRegex("\\s*([a-zA-Z]+)\\s*([0-9])-([0-9]+)([G])\\s*");
        while(regex_search(line, what, basisRegex)) { // n 4-31G
            atomTypeAbbreviation = string(what[1]);
//            int nCoreOrbitals(stoi(string(what[2])));
            for(char& c : string(what[3])) {
                int nValenceOrbital = c - '0';
                nValenceOrbitals.push_back(nValenceOrbital);
            }
//            string basisSetType(what[4]);

//            cout << "Atom: " << atomTypeAbbreviation << " nCore: " << nCoreOrbitals << " type: " << basisSetType << " nValence: ";
//            for(int nValenceOrbital : nValenceOrbitals) {
//                cout << nValenceOrbital << " ";
//            }
//            cout << endl;
//            cout << "AtomType: " << HF::abbreviationToAtomType(atomTypeAbbreviation) << endl;
            m_atomType = HF::abbreviationToAtomType(atomTypeAbbreviation);

            break;
        }
        regex orbitalRegex("\\s*([0-9])\\s*([spdf])\\s*");
        if(regex_search(line, what, orbitalRegex)) { // 4 s
            mergePrimitivesIntoContracted();
            if(what[2] == "s") {
                m_currentOrbitalType = HF::sOrbitalType;
            } else if(what[2] == "p") {
                m_currentOrbitalType = HF::pOrbitalType;
            } else if(what[2] == "d") {
                m_currentOrbitalType = HF::dOrbitalType;
                cout << "ERROR: d orbitals not yet supported" << endl;
                return false;
            } else if(what[2] == "f") {
                m_currentOrbitalType = HF::fOrbitalType;
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
            m_collectedPrimitiveBasisFunctions.push_back(primitive);
        }
    }
    mergePrimitivesIntoContracted();
//    for(GaussianContractedOrbital orbital : m_contractedBasisFunctions) {
//        cout << "Contracted:" << endl;
//        for(GaussianPrimitiveOrbital primitive : orbital.primitiveBasisFunctions()) {
//            cout << primitive.weight() << " * " << primitive.exponent() << endl;
//        }
//    }
    cout << "Done!" << endl;
    return true;
}

const vector<GaussianContractedOrbital> &TurboMoleParser::contractedBasisFunctions() const
{
    return m_contractedBasisFunctions;
}

void TurboMoleParser::mergePrimitivesIntoContracted()
{
    if(m_collectedPrimitiveBasisFunctions.size() > 0) {
        switch(m_currentOrbitalType) {
        case HF::sOrbitalType: {
            GaussianContractedOrbital contractedOrbital;
            for(GaussianPrimitiveOrbital primitive : m_collectedPrimitiveBasisFunctions) {
                double weightAdjusted = primitive.weight() * pow(2 * primitive.exponent() / M_PI, 3.0/4.0);
                primitive.setWeight(weightAdjusted);
                contractedOrbital.addPrimitiveBasisFunction(primitive);
            }
            m_contractedBasisFunctions.push_back(contractedOrbital);
            break;
        }
        case HF::pOrbitalType: {
            int dim = 3;
            for(int i = 0; i < dim; i++) {
                GaussianContractedOrbital contractedOrbital;
                for(GaussianPrimitiveOrbital primitive : m_collectedPrimitiveBasisFunctions) {
                    double weightAdjusted = primitive.weight() * pow(2 * primitive.exponent() / M_PI, 3.0/4.0) * 2 * sqrt(primitive.exponent());
                    primitive.setWeight(weightAdjusted);
                    primitive.setXExponent(i == 0);
                    primitive.setYExponent(i == 1);
                    primitive.setZExponent(i == 2);
                    contractedOrbital.addPrimitiveBasisFunction(primitive);
                }
                m_contractedBasisFunctions.push_back(contractedOrbital);
            }
            break;
            // TODO Support d and f type orbitals
        }
        case HF::dOrbitalType: {
            cout << "ERROR: d orbitals not yet supported" << endl;
            break;
        }
        case HF::fOrbitalType: {
            cout << "ERROR: f orbitals not yet supported" << endl;
            break;
        }
        default: {
            cout << "ERROR: Unknown orbital type." << endl;
            break;
        }
        }
        m_collectedPrimitiveBasisFunctions.clear();
    }
}
HF::AtomType TurboMoleParser::atomType() const
{
    return m_atomType;
}
