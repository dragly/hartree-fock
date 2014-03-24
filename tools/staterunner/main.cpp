#include <iostream>
#include <iomanip>
#include <fstream>
#include <libconfig.h++>

#include <hartreefocksolver.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>

#include <H5Cpp.h>

using namespace std;
using namespace H5;

int main(int argc, char* argv[])
{
    if(argc < 2) {
        cout << "Error: No file provided." << endl;
        cout << "Usage: programname <filename>" << endl;
    }
    /* First structure  and dataset*/
    struct AtomData {
        double x;
        double y;
        double z;
        double partialCharge;
    };

    CompType atomCompound( sizeof(AtomData) );
    atomCompound.insertMember("x", HOFFSET(AtomData, x), PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("y", HOFFSET(AtomData, y), PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("z", HOFFSET(AtomData, z), PredType::NATIVE_DOUBLE);


    struct AtomMetaData {
        int    type;
        char basisName[64];
    };
    H5::StrType string_type(H5::PredType::C_S1, 64);
    CompType atomMetaCompound( sizeof(AtomMetaData) );
    atomMetaCompound.insertMember( "type", HOFFSET(AtomMetaData, type), PredType::NATIVE_INT);
    atomMetaCompound.insertMember( "basisName", HOFFSET(AtomMetaData, basisName), string_type);

    cout << "Opening file " << argv[1] << endl;
    H5File inFile(argv[1], H5F_ACC_RDONLY );
    cout << inFile.getNumObjs() << " states found" << endl;

    cout << "Getting atom metadata" << endl;
    Group rootGroup(inFile.openGroup("/"));
    string metaInformationName = "atomMetaInformation";
    Attribute atomMetaAttribute(rootGroup.openAttribute(metaInformationName));
    hsize_t dims[1];
    atomMetaAttribute.getSpace().getSimpleExtentDims(dims);
    int nAtoms = dims[0];
    cout << nAtoms << endl;

    AtomMetaData atomMetaData[nAtoms];
    atomMetaAttribute.read(atomMetaCompound, atomMetaData);

    for(int i = 0; i < nAtoms; i++) {
        cout << atomMetaData[i].type << endl;
        cout << atomMetaData[i].basisName << endl;
    }

    cout << "Opening new file for writing..." << endl;
    H5File outFile("results.h5", H5F_ACC_TRUNC );
    Group rootGroupOut(outFile.openGroup("/"));

    cout << "Writing atom metadata..." << endl;
    Attribute atomMetaAttributeOut(rootGroupOut.createAttribute(metaInformationName,
                                                     atomMetaCompound,
                                                     DataSpace(atomMetaAttribute.getSpace())));

    atomMetaAttributeOut.write(atomMetaCompound, atomMetaData);

    for(int stateIndex = 0; stateIndex < int(inFile.getNumObjs()); stateIndex++) {
        string objectName = inFile.getObjnameByIdx(stateIndex);
        if(objectName.find("state") == string::npos) {
            continue;
        }
        cout << "Opening " << objectName << endl;
        DataSet atomDataSet(inFile.openDataSet(objectName));
        hsize_t dims2[1];
        atomDataSet.getSpace().getSimpleExtentDims(dims2);
        int nAtoms2 = dims2[0];
        cout << nAtoms2 << endl;

        if(nAtoms != nAtoms2) {
            cerr << "Error! The number of atoms in " << objectName << " (nAtoms = " << nAtoms2 << ") does not match "
                 << "the number of atoms in the metadata (= " << nAtoms << ")" << endl;
            cerr << "Cannot continue" << endl;
            exit(0);
        }

        AtomData *atoms = new AtomData[nAtoms2];
        atomDataSet.read(atoms, atomCompound);

        cout << "Writing atoms to new file..." << endl;

        DataSet atomDataSetOut(outFile.createDataSet(objectName, atomCompound, DataSpace(atomDataSet.getSpace())));
        atomDataSetOut.write(atoms, atomCompound);

        GaussianSystem system;
        for(int i = 0; i < nAtoms2; i++) {
            string fileName = "";
            if(atomMetaData[i].type == 8 && !strcmp(atomMetaData[i].basisName, "3-21G")) {
                fileName = "oxygen321g.tm";
            } else if(atomMetaData[i].type == 1 && !strcmp(atomMetaData[i].basisName, "3-21G")) {
                fileName = "hydrogen321g.tm";
            }
            cout << fileName << endl;
            system.addCore(GaussianCore({ atoms[i].x, atoms[i].y, atoms[i].z}, fileName));
            cout << atomMetaData[i].type << endl;
            cout << atoms[i].x << endl;
            cout << atoms[i].y << endl;
            cout << atoms[i].z << endl;
        }
        HartreeFockSolver solver(&system);
        solver.setNIterationsMax(1e3);
        cout << "Solving..." << endl;
        solver.solve();
        cout << "Energy: " << solver.energy() << endl;

        double energy = solver.energy();

        Attribute energyAttribute(atomDataSetOut.createAttribute("energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
        energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);
        outFile.flush(H5F_SCOPE_GLOBAL);

        delete atoms;
    }
    return 0;
}

