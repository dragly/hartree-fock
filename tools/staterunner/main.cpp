#include <iostream>
#include <iomanip>
#include <fstream>
#include <libconfig.h++>

#include <solvers/unrestrictedhartreefocksolver.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <boost/mpi.hpp>

#include <H5Cpp.h>

using namespace std;
using namespace H5;
namespace mpi = boost::mpi;

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

    string inFileName = argv[1];

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

    vector<int> stateIDs;
    if(world.rank() == 0) {
        H5File inFile(inFileName, H5F_ACC_RDONLY );
        H5::Group statesGroup(inFile.openGroup("/states"));
        int nStates = statesGroup.getNumObjs();
        vector<int> allStates;
        allStates.reserve(nStates);
        for(int i = 0; i < nStates; i++) {
            allStates.push_back(i);
        }
        random_shuffle(allStates.begin(), allStates.end());
        for(int p = 0; p < world.size(); p++) {
            vector<int> states;
            for(int i = blockLow(p, world.size(), allStates.size());
                i <= blockHigh(p, world.size(), allStates.size());
                i++) {

                states.push_back(allStates.at(i));

            }
            if(p == 0) {
                stateIDs = states;
            } else {
                world.send(p, 0, states);
            }
        }
        inFile.close();
    } else {
        world.recv(0, 0, stateIDs);
    }
    cout << "Rank " << world.rank() << ": I have " << stateIDs.size() << " states to dig!" << endl;
    world.barrier();

    // Copy HDF5 file
    stringstream outFileName;
    outFileName << inFileName << ".results." << setfill('0') << setw(4) << world.rank();

    // Open outFile for writing
    H5File outFile(outFileName.str(), H5F_ACC_TRUNC);

    for(int p = 0; p < world.size(); p++) {
        world.barrier();
        if(p != world.rank()) {
            continue;
        }
        // Open inFile for reading
        H5File inFile(inFileName, H5F_ACC_RDONLY);

        // Copy states to outFile

        for(int i = 0; i < int(inFile.getNumObjs()); i++) {
            string objectName = inFile.getObjnameByIdx(i);
            if(objectName == string("states")) {
                H5::Group statesGroup(inFile.openGroup("/states"));
                H5::Group statesGroupOut(outFile.createGroup("/states"));
                for(int stateID : stateIDs) {
                    string stateName = statesGroup.getObjnameByIdx(stateID);
                    H5Ocopy(statesGroup.getId(), stateName.c_str(), statesGroupOut.getId(), stateName.c_str(), H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
                }
            } else {
                H5Ocopy(inFile.getId(), objectName.c_str(), outFile.getId(), objectName.c_str(), H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
            }
        }
        inFile.close();
    }

    H5::Group rootGroup(outFile.openGroup("/"));
    string metaInformationName = "atomMeta";
    DataSet atomMeta(rootGroup.openDataSet(metaInformationName));
    hsize_t dims[1];
    atomMeta.getSpace().getSimpleExtentDims(dims);
    int nAtoms = dims[0];

    AtomMetaData atomMetaData[nAtoms];
    atomMeta.read(atomMetaData, atomMetaCompound);

    H5::Group statesGroup(outFile.openGroup("/states"));

    int nTotal = statesGroup.getNumObjs();
    int currentState = 0;
    for(int stateID = 0; stateID < nTotal; stateID++) {
        if(world.rank() == 0) {
            cout << "Rank " << world.rank()
                 << ", progress: " << fixed << setprecision(2) << double(currentState) / nTotal * 100 << " %"
                 << ", time: " << timer.elapsed() << " s"
                 << endl;
        }
        string stateName = statesGroup.getObjnameByIdx(stateID);
        DataSet stateDataSet(statesGroup.openDataSet(stateName));
        hsize_t dims2[1];
        stateDataSet.getSpace().getSimpleExtentDims(dims2);
        int nAtoms2 = dims2[0];

        if(nAtoms != nAtoms2) {
            cerr << "Error! The number of atoms in " << stateName << " (nAtoms = " << nAtoms2 << ") does not match "
                 << "the number of atoms in the metadata (= " << nAtoms << ")" << endl;
            cerr << "Cannot continue" << endl;
            throw std::logic_error("Mismatching number of atoms");
        }

        AtomData *atoms = new AtomData[nAtoms2];
        stateDataSet.read(atoms, atomCompound);

        GaussianSystem system;
        for(int i = 0; i < nAtoms2; i++) {
            stringstream basisFile;
            basisFile << "atom_" << atomMetaData[i].type << "_basis_" << atomMetaData[i].basisName << ".tm";
            string fileName = basisFile.str();
            system.addCore(GaussianCore({ atoms[i].x, atoms[i].y, atoms[i].z}, fileName));
        }
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.95);
        solver.setConvergenceTreshold(1e-9);
        solver.solve();

        double energy = solver.energy();

        Attribute energyAttribute(stateDataSet.createAttribute("energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
        energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);
        outFile.flush(H5F_SCOPE_GLOBAL);

        currentState++;

        delete atoms;
    }
    world.barrier();
    if(world.rank() == 0) {
        cout << "Total time "  << fixed << setprecision(2) << timer.elapsed() << " s" << endl;
    }
    outFile.close();
    return 0;
}

