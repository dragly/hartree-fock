#include <solvers/restrictedhartreefocksolver.h>
#include <solvers/unrestrictedhartreefocksolver.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <math/vector3.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <H5Cpp.h>

using namespace std;

void operator >> (const YAML::Node& node, Vector3& v)
{
    double x;
    double y;
    double z;
    node[0] >> x;
    node[1] >> y;
    node[2] >> z;
    v = Vector3(x,y,z);
}

class Method {
public:
    enum MethodName {
        Unknown,
        Restricted,
        Unrestricted
    };
};


class Output {
public:
    enum OutputName {
        Unknown,
        Energy,
        Density,
        ElectrostaticPotential
    };
};

int main(int argc, char* argv[])
{
    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        cerr << "Usage: hartree-fock <config.yaml>" << endl;
        return 1;
    }
    string outputPath = "";
    if(argc >= 3) {
        outputPath = argv[2];
        outputPath += "/";
        cout << "Output will be written to " << outputPath << endl;
    }
    ifstream fin(argv[1]);
    if(fin.fail()) {
        cerr << "Could not open the configuration file " << argv[1] << endl;
        return 1;
    }
    YAML::Parser parser(fin);

    YAML::Node rootNode;
    parser.GetNextDocument(rootNode);
    unsigned int nAtoms = rootNode["atoms"].size();
    string defaultBasis = "3-21G";
    try {
        rootNode["basis"] >> defaultBasis;
        cout << "Using basis " << defaultBasis << endl;
    } catch( YAML::TypedKeyNotFound<std::string> ) {
        cout << "Default basis not found, expecting atom to have basis name or assuming " << defaultBasis << endl;
    }

    struct AtomData {
        double x;
        double y;
        double z;
        double partialCharge;
    };

    H5::CompType atomCompound( sizeof(AtomData) );
    atomCompound.insertMember("x", HOFFSET(AtomData, x), H5::PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("y", HOFFSET(AtomData, y), H5::PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("z", HOFFSET(AtomData, z), H5::PredType::NATIVE_DOUBLE);

    struct AtomMetaData {
        int    type;
        char basisName[64];
    };

    H5::CompType atomMetaCompound( sizeof(AtomMetaData) );
    atomMetaCompound.insertMember("type", HOFFSET(AtomMetaData, type), H5::PredType::NATIVE_INT);
    atomMetaCompound.insertMember("basisName", HOFFSET(AtomMetaData, basisName), H5::StrType(H5::PredType::C_S1, 64));

    H5::H5File outFile(outputPath + "results.h5", H5F_ACC_TRUNC);
    hsize_t dim[] = {nAtoms};
    H5::DataSet stateDataSet = outFile.createDataSet("state", atomCompound, H5::DataSpace(1, dim));
    H5::DataSet atomMetaDataSet(outFile.createDataSet("atomMeta", atomMetaCompound, H5::DataSpace(1, dim)));

    AtomData atoms[nAtoms];
    AtomMetaData atomsMeta[nAtoms];

    GaussianSystem system;
    Method::MethodName method = Method::Unknown;
    vector<Output::OutputName> outputs;
    for(YAML::Iterator it=rootNode.begin();it!=rootNode.end();++it) {
        string rootKey;
        it.first() >> rootKey;
        if(rootKey == "atoms") {
            int atomCounter = 0;
            const YAML::Node &atomsNode = it.second();
            for(YAML::Iterator it2=atomsNode.begin();it2!=atomsNode.end();++it2) {
                const YAML::Node &atomNode = *it2;
                string typeAbbreviation;
                atomNode["type"] >> typeAbbreviation;
                Vector3 position;
                atomNode["position"] >> position;
                string basis;
                try {
                    atomNode["basis"] >> basis;
                } catch( YAML::TypedKeyNotFound<std::string> ) {
                    basis = defaultBasis;
                }

                basis = HF::escapeBasis(basis);

                strncpy(atomsMeta[atomCounter].basisName, basis.c_str(), 63);
                atomsMeta[atomCounter].type = HF::abbreviationToAtomType(typeAbbreviation);
                system.addCore(GaussianCore(position, typeAbbreviation, basis));
                atoms[atomCounter].x = position.x();
                atoms[atomCounter].y = position.y();
                atoms[atomCounter].z = position.z();
                atomCounter++;
            }
        } else if(rootKey == "method") {
            string methodName;
            it.second() >> methodName;
            if(methodName == "unrestricted") {
                method = Method::Unrestricted;
            } else if(methodName == "restricted") {
                method = Method::Restricted;
            } else {
                cerr << "Error: Unknown method name " << methodName << endl;
                return 1;
            }
        } else if(rootKey == "electronsDown" || rootKey == "electronsUp") {
            int electronsDown;
            it.second() >> electronsDown;
            system.setNParticlesDown(electronsDown);
        } else if(rootKey == "output") {
            const YAML::Node &outputNode = it.second();
            for(YAML::Iterator it2=outputNode.begin();it2!=outputNode.end();++it2) {
                string outputName;
                *it2 >> outputName;
                if(outputName == "energy") {
                    outputs.push_back(Output::Energy);
                } else if(outputName == "density") {
                    outputs.push_back(Output::Density);
                } else if(outputName == "electrostatic_potential") {
                    outputs.push_back(Output::ElectrostaticPotential);
                }
            }
        }
    }
    stateDataSet.write(atoms, atomCompound);
    atomMetaDataSet.write(atomsMeta, atomMetaCompound);

    double energy = 0.0;
    mat coefficientsUp;
    mat coefficientsDown;
    if(method == Method::Unrestricted) {
        cout << "Setting up the unrestricted Hartree-Fock solver..." << endl;
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-10);
        solver.setDensityMixFactor(0.95);
        solver.setNIterationsMax(1e4);
        cout << "Solving..." << endl;
        solver.solve();
        energy = solver.energy();
        cout << "Energy: " << setprecision(16) << solver.energy() << endl;

        coefficientsUp = solver.coeffcientMatrixUp();
        coefficientsDown = solver.coeffcientMatrixDown();

    } else {
        cout << "Setting up the restricted Hartree-Fock solver..." << endl;
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setDensityMixFactor(0.95);
        solver.setNIterationsMax(1e4);
        cout << "Solving..." << endl;
        solver.solve();
        energy = solver.energy();
        cout << "Energy: " << setprecision(16) << solver.energy() << endl;

        coefficientsUp = solver.coefficientMatrix();
    }
    cout << "Writing requested output to file:" << endl;
    for(Output::OutputName output : outputs) {
        if(output == Output::Energy) {
            cout << "Writing energy..." << endl;
            H5::Attribute energyAttribute(stateDataSet.createAttribute("energy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR));
            energyAttribute.write(H5::PredType::NATIVE_DOUBLE, &energy);
        } else if(output == Output::Density) {
            cout << "Calculating density..." << endl;

            mat coefficientsAll;
            if(method == Method::Unrestricted) {
                coefficientsAll = join_rows(coefficientsUp, coefficientsDown);
            } else if(method == Method::Restricted) {
                coefficientsAll = coefficientsUp;
            }
            vec x = linspace(-3, 3, 100);
            vec y = linspace(-3, 3, 100);
            vec z = linspace(-3, 3, 100);
            cube totalDensity = zeros(x.n_elem, y.n_elem, z.n_elem);
            field<cube> orbitalDensities;
            orbitalDensities.set_size(system.nParticles());
            for(cube& orbitalDensity : orbitalDensities) {
                orbitalDensity = zeros(x.n_elem, y.n_elem, z.n_elem);
            }
            for(uint i = 0; i < x.n_elem; i++) {
                for(uint j = 0; j < y.n_elem; j++) {
                    for(uint k = 0; k < z.n_elem; k++) {
                        Vector3 position(x(i), y(j), z(k));
                        rowvec density = system.orbitalDensities(coefficientsAll, position);
                        for(uint orbital = 0; orbital < orbitalDensities.size(); orbital++) {
                            cube &orbitalDensity = orbitalDensities(orbital);
                            orbitalDensity(k,j,i) = density(orbital);
                        }
                        totalDensity(k,j,i) = sum(density);
                    }
                }
            }
            for(uint orbital = 0; orbital < orbitalDensities.size(); orbital++) {
                cube &orbitalDensity = orbitalDensities(orbital);
                orbitalDensity.save(outputPath + "orbital_density_" + to_string(orbital) + ".h5", hdf5_binary);
            }
            totalDensity.save(outputPath + "density.h5", hdf5_binary);
        } else if(output == Output::ElectrostaticPotential) {
            cout << "Calculating electrostatic potential..." << endl;
            if(method == Method::Unrestricted) {

                vec x = linspace(-5, 5, 50);
                vec y = linspace(-5, 5, 50);
                vec z = linspace(-5, 5, 50);
                cube electrostaticPotential = zeros(x.n_elem, y.n_elem, z.n_elem);
                mat coefficientsAll = join_rows(coefficientsUp, coefficientsDown);
                for(uint i = 0; i < x.n_elem; i++) {
                    for(uint j = 0; j < y.n_elem; j++) {
                        for(uint k = 0; k < z.n_elem; k++) {
                            Vector3 position(x(i), y(j), z(k));
                            electrostaticPotential(k,j,i) = system.electrostaticPotential(coefficientsAll, position);
                        }
                    }
                }
                electrostaticPotential.save(outputPath + "electrostatic_potential.h5", hdf5_binary);
            }
        }
    }
    outFile.close();
    cout << "Done." << endl;
    return 0;
}

