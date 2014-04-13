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
        regex basisRegex("^\\s*([a-zA-Z]+)");
        while(regex_search(line, what, basisRegex)) { // n 4-31G
            atomTypeAbbreviation = string(what[1]);
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
            } else if(what[2] == "f") {
                m_currentOrbitalType = HF::fOrbitalType;
            } else {
                cerr << "Unknown orbital type " << string(what[1]) << endl;
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
    //    cout << "Done!" << endl;
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
        }
        case HF::dOrbitalType: {
            umat exponents;
            exponents << 2 << 0 << 0 << endr
                      << 0 << 2 << 0 << endr
                      << 0 << 0 << 2 << endr
                      << 1 << 1 << 0 << endr
                      << 1 << 0 << 1 << endr
                      << 0 << 1 << 1 << endr;
            for(int i = 0; i < int(exponents.n_rows); i++) {
                GaussianContractedOrbital contractedOrbital;
                for(GaussianPrimitiveOrbital primitive : m_collectedPrimitiveBasisFunctions) {
                    double weightAdjusted = primitive.weight() * normalizationFactor(primitive.exponent(), exponents.row(i));
                    primitive.setWeight(weightAdjusted);
                    primitive.setXExponent(exponents(i, 0));
                    primitive.setYExponent(exponents(i, 1));
                    primitive.setZExponent(exponents(i, 2));
                    contractedOrbital.addPrimitiveBasisFunction(primitive);
                }
                m_contractedBasisFunctions.push_back(contractedOrbital);
            }
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

//TODO: Change urowvec to std::array<int, 3>
double TurboMoleParser::normalizationFactor(double exp, urowvec pows)
{
    int i = pows.at(0);
    int j = pows.at(1);
    int k = pows.at(2);
    return pow((2*exp/M_PI),0.75)*sqrt(pow(8*exp,i+j+k)*factorial(i)*factorial(j)*factorial(k)/(factorial(2*i)*factorial(2*j)*factorial(2*k)));
}

int TurboMoleParser::factorial(int n)
{
    double value = 1;
    double i = 1;

    while(i < n){
        i += 1;
        value *= i;
    }
    return value;
}
