#ifndef HF_H
#define HF_H

#include <map>
#include <unordered_map>
#include <string>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>

struct iequal_to
        : std::binary_function<std::string, std::string, bool>
{
    bool operator()(std::string const& x,
                    std::string const& y) const
    {
        return boost::algorithm::iequals(x, y, std::locale());
    }
};

struct ihash
        : std::unary_function<std::string, std::size_t>
{
    std::size_t operator()(std::string const& x) const
    {
        std::size_t seed = 0;
        std::locale locale;

        for(std::string::const_iterator it = x.begin();
            it != x.end(); ++it)
        {
            boost::hash_combine(seed, std::toupper(*it, locale));
        }

        return seed;
    }
};

class HF
{
public:
    enum AtomOrbitalType {
        sOrbitalType,
        pOrbitalType,
        dOrbitalType,
        fOrbitalType
    };

    enum AtomType {
        Unknown,
        Hydrogen,
        Helium,
        Lithium,
        Beryllium,
        Boron,
        Carbon,
        Nitrogen,
        Oxygen,
        Fluorine,
        Neon,
        Sodium,
        Magnesium,
        Aluminum,
        Silicon,
        Phosphorus,
        Sulfur,
        Chlorine,
        Argon,
        Potassium,
        Calcium,
        Scandium,
        Titanium,
        Vanadium,
        Chromium,
        Manganese,
        Iron,
        Cobalt,
        Nickel,
        Copper,
        Zinc,
        Gallium,
        Germanium,
        Arsenic,
        Selenium,
        Bromine,
        Krypton,
        Rubidium,
        Strontium,
        Yttrium,
        Zirconium,
        Niobium,
        Molybdenum,
        Technetium,
        Ruthenium,
        Rhodium,
        Palladium,
        Silver,
        Cadmium,
        Indium,
        Tin,
        Antimony,
        Tellurium,
        Iodine,
        Xenon,
        Cesium,
        Barium,
        Lanthanum,
        Cerium,
        Praseodymium,
        Neodymium,
        Promethium,
        Samarium,
        Europium,
        Gadolinium,
        Terbium,
        Dysprosium,
        Holmium,
        Erbium,
        Thulium,
        Ytterbium,
        Lutetium,
        Hafnium,
        Tantalum,
        Tungsten,
        Rhenium,
        Osmium,
        Iridium,
        Platinum,
        Gold,
        Mercury,
        Thallium,
        Lead,
        Bismuth,
        Polonium,
        Astatine,
        Radon,
        Francium,
        Radium,
        Actinium,
        Thorium,
        Protactinium,
        Uranium,
        Neptunium,
        Plutonium,
        Americium,
        Curium,
        Berkelium,
        Californium,
        Einsteinium,
        Fermium,
        Mendelevium,
        Nobelium,
        Lawrencium,
        Rutherfordium,
        Dubnium,
        Seaborgium,
        Bohrium,
        Hassium,
        Meitnerium,
        Darmstadtium,
        Roentgenium,
        Copernicium,
        Ununtrium,
        Flerovium,
        Ununpentium,
        Livermorium,
        Ununseptium,
        Ununoctium
    };
    static AtomType abbreviationToAtomType(std::string abbreviation) {
        std::unordered_map<std::string, AtomType, ihash, iequal_to> abbreviationMap {
            {"H",Hydrogen},
            {"He",Helium},
            {"Li",Lithium},
            {"Be",Beryllium},
            {"B",Boron},
            {"C",Carbon},
            {"N",Nitrogen},
            {"O",Oxygen},
            {"F",Fluorine},
            {"Ne",Neon},
            {"Na",Sodium},
            {"Mg",Magnesium},
            {"Al",Aluminum},
            {"Si",Silicon},
            {"P",Phosphorus},
            {"S",Sulfur},
            {"Cl",Chlorine},
            {"Ar",Argon},
            {"K",Potassium},
            {"Ca",Calcium},
            {"Sc",Scandium},
            {"Ti",Titanium},
            {"V",Vanadium},
            {"Cr",Chromium},
            {"Mn",Manganese},
            {"Fe",Iron},
            {"Co",Cobalt},
            {"Ni",Nickel},
            {"Cu",Copper},
            {"Zn",Zinc},
            {"Ga",Gallium},
            {"Ge",Germanium},
            {"As",Arsenic},
            {"Se",Selenium},
            {"Br",Bromine},
            {"Kr",Krypton},
            {"Rb",Rubidium},
            {"Sr",Strontium},
            {"Y",Yttrium},
            {"Zr",Zirconium},
            {"Nb",Niobium},
            {"Mo",Molybdenum},
            {"Tc",Technetium},
            {"Ru",Ruthenium},
            {"Rh",Rhodium},
            {"Pd",Palladium},
            {"Ag",Silver},
            {"Cd",Cadmium},
            {"In",Indium},
            {"Sn",Tin},
            {"Sb",Antimony},
            {"Te",Tellurium},
            {"I",Iodine},
            {"Xe",Xenon},
            {"Cs",Cesium},
            {"Ba",Barium},
            {"La",Lanthanum},
            {"Ce",Cerium},
            {"Pr",Praseodymium},
            {"Nd",Neodymium},
            {"Pm",Promethium},
            {"Sm",Samarium},
            {"Eu",Europium},
            {"Gd",Gadolinium},
            {"Tb",Terbium},
            {"Dy",Dysprosium},
            {"Ho",Holmium},
            {"Er",Erbium},
            {"Tm",Thulium},
            {"Yb",Ytterbium},
            {"Lu",Lutetium},
            {"Hf",Hafnium},
            {"Ta",Tantalum},
            {"W",Tungsten},
            {"Re",Rhenium},
            {"Os",Osmium},
            {"Ir",Iridium},
            {"Pt",Platinum},
            {"Au",Gold},
            {"Hg",Mercury},
            {"Tl",Thallium},
            {"Pb",Lead},
            {"Bi",Bismuth},
            {"Po",Polonium},
            {"At",Astatine},
            {"Rn",Radon},
            {"Fr",Francium},
            {"Ra",Radium},
            {"Ac",Actinium},
            {"Th",Thorium},
            {"Pa",Protactinium},
            {"U",Uranium},
            {"Np",Neptunium},
            {"Pu",Plutonium},
            {"Am",Americium},
            {"Cm",Curium},
            {"Bk",Berkelium},
            {"Cf",Californium},
            {"Es",Einsteinium},
            {"Fm",Fermium},
            {"Md",Mendelevium},
            {"No",Nobelium},
            {"Lr",Lawrencium},
            {"Rf",Rutherfordium},
            {"Db",Dubnium},
            {"Sg",Seaborgium},
            {"Bh",Bohrium},
            {"Hs",Hassium},
            {"Mt",Meitnerium},
            {"Ds",Darmstadtium},
            {"Rg",Roentgenium},
            {"Cn",Copernicium},
            {"Uut",Ununtrium},
            {"Fl",Flerovium},
            {"Uup",Ununpentium},
            {"Lv",Livermorium},
            {"Uus",Ununseptium},
            {"Uuo",Ununoctium}
        };
        if(abbreviationMap.find(abbreviation) != abbreviationMap.end()) {
            return abbreviationMap.find(abbreviation)->second;
        } else {
            return Unknown;
        }
    }
};

#endif // HF_H
