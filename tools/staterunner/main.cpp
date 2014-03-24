#include <iostream>
#include <iomanip>
#include <fstream>
#include <libconfig.h++>

#include <hartreefocksolver.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <boost/mpi.hpp>

#include <H5Cpp.h>

using namespace std;
using namespace H5;
using namespace boost;

int blockLow(int id, int np, int n) {
    return (id * n) / np;
}

int blockHigh(int id, int np, int n) {
    return blockLow(id + 1, np, n) - 1;
}

int blockSize(int id, int p, int n) {
    return blockLow(id + 1,p,n) - blockLow(id,p,n);
}

int main(int argc, char* argv[])
{
    mpi::environment env;
    mpi::communicator world;
    mpi::timer timer;
    timer.restart();
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

    H5File inFile(argv[1], H5F_ACC_RDONLY );

    Group rootGroup(inFile.openGroup("/"));
    string metaInformationName = "atomMetaInformation";
    Attribute atomMetaAttribute(rootGroup.openAttribute(metaInformationName));
    hsize_t dims[1];
    atomMetaAttribute.getSpace().getSimpleExtentDims(dims);
    int nAtoms = dims[0];
    cout << nAtoms << endl;

    AtomMetaData atomMetaData[nAtoms];
    atomMetaAttribute.read(atomMetaCompound, atomMetaData);

    stringstream outFileName;
    outFileName << "results.h5." << setfill('0') << setw(4) << world.rank();
    H5File outFile(outFileName.str(), H5F_ACC_TRUNC );
    Group rootGroupOut(outFile.openGroup("/"));

    cout << "Writing atom metadata..." << endl;
    Attribute atomMetaAttributeOut(rootGroupOut.createAttribute(metaInformationName,
                                                                atomMetaCompound,
                                                                DataSpace(atomMetaAttribute.getSpace())));

    atomMetaAttributeOut.write(atomMetaCompound, atomMetaData);

    vector<string> allStates;
    for(int stateIndex = 0; stateIndex < int(inFile.getNumObjs()); stateIndex++) {
        string objectName = inFile.getObjnameByIdx(stateIndex);
        if(objectName.find("state") == string::npos) {
            continue;
        }
        allStates.push_back(objectName);
    }
    vector<string> states;
    for(int i = blockLow(world.rank(), world.size(), allStates.size());
        i <= blockHigh(world.rank(), world.size(), allStates.size());
        i++) {
        states.push_back(allStates.at(i));
    }
    cout << "I got " << states.size() << " states to handle." << endl;
    int nTotal = states.size();
    int currentState = 0;
    for(const string& stateName : states) {
        if(world.rank() == 0) {
            cout << "Rank " << world.rank()
                 << ", progress: " << fixed << setprecision(2) << double(currentState) / nTotal * 100 << " %"
                 << ", time: " << timer.elapsed() << " s"
                 << endl;
        }
        DataSet atomDataSet(inFile.openDataSet(stateName));
        hsize_t dims2[1];
        atomDataSet.getSpace().getSimpleExtentDims(dims2);
        int nAtoms2 = dims2[0];

        if(nAtoms != nAtoms2) {
            cerr << "Error! The number of atoms in " << stateName << " (nAtoms = " << nAtoms2 << ") does not match "
                 << "the number of atoms in the metadata (= " << nAtoms << ")" << endl;
            cerr << "Cannot continue" << endl;
            exit(0);
        }

        AtomData *atoms = new AtomData[nAtoms2];
        atomDataSet.read(atoms, atomCompound);

        DataSet atomDataSetOut(outFile.createDataSet(stateName, atomCompound, DataSpace(atomDataSet.getSpace())));
        atomDataSetOut.write(atoms, atomCompound);

        GaussianSystem system;
        for(int i = 0; i < nAtoms2; i++) {
            string fileName = "";
            if(atomMetaData[i].type == 8 && !strcmp(atomMetaData[i].basisName, "3-21G")) {
                fileName = "oxygen321g.tm";
            } else if(atomMetaData[i].type == 1 && !strcmp(atomMetaData[i].basisName, "3-21G")) {
                fileName = "hydrogen321g.tm";
            } else if(atomMetaData[i].type == 1 && !strcmp(atomMetaData[i].basisName, "6-311G")) {
                fileName = "hydrogen6311g.tm";
            }
            system.addCore(GaussianCore({ atoms[i].x, atoms[i].y, atoms[i].z}, fileName));
        }
        HartreeFockSolver solver(&system);
        solver.setNIterationsMax(1e3);
        solver.solve();

        double energy = solver.energy();

        Attribute energyAttribute(atomDataSetOut.createAttribute("energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
        energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);
        outFile.flush(H5F_SCOPE_GLOBAL);

        currentState++;

        delete atoms;
    }
    world.barrier();
    if(world.rank() == 0) {
        cout << "Total time "  << fixed << setprecision(2) << timer.elapsed() << " s" << endl;
    }
    return 0;
}

