#include <iostream>
#include <iomanip>
#include <fstream>
#include <libconfig.h++>
#include <string.h>

#include <solvers/unrestrictedhartreefocksolver.h>
#include <solvers/restrictedhartreefocksolver.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <boost/mpi.hpp>
#include <yaml-cpp/yaml.h>

#include <H5Cpp.h>

using namespace std;
using namespace H5;
namespace mpi = boost::mpi;

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
        cout << "Usage: programname <infile> [outfile]" << endl;
        exit(0);
    }

    string inFileName = argv[1];

    stringstream outFileName;
    if(argc < 3) {
        outFileName << "states";
    } else {
        outFileName << argv[2];
    }

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
    atomMetaCompound.insertMember("type", HOFFSET(AtomMetaData, type), PredType::NATIVE_INT);
    atomMetaCompound.insertMember("basisName", HOFFSET(AtomMetaData, basisName), string_type);

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

    outFileName << "." << setfill('0') << setw(4) << world.rank() << ".h5";

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

    string method("");
    Attribute hartreeFockMethodAttribute = atomMeta.openAttribute("method");
    hartreeFockMethodAttribute.read(hartreeFockMethodAttribute.getDataType(), method);

    mat coefficientMatrixUp;
    mat coefficientMatrixDown;

    // Precalculate energy offset only if using unrestricted.
    cout << "Calculating energy offset..." << endl;
    Attribute energyOffsetAttribute = atomMeta.createAttribute("energyOffset", PredType::NATIVE_DOUBLE, H5S_SCALAR);
    double energyOffset = 0.0;
    for(int i = 0; i < nAtoms; i++) {
        GaussianSystem system;
        system.addCore(GaussianCore({0,0,0}, atomMetaData[i].type, atomMetaData[i].basisName));
        // Restricted HF cannot be done on one atom currently.
        // Restricted tends to have convergence problems with high distances anyway,
        // so just use the energy offset from UHF in case we want to plot the two together with offset
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setDiisEnabled(false);
        solver.setNIterationsMax(1e3);
        solver.setDensityMixFactor(0.95);
        solver.setConvergenceTreshold(1e-9);
        solver.solve();
        energyOffset += solver.energy();
    }
    energyOffsetAttribute.write(PredType::NATIVE_DOUBLE, &energyOffset);
    cout << "Energy offset: " << energyOffset << endl;
    // Done precalculating energy offset

    // Precalculate ground state
    cout << "Calculating ground state to get coefficients..." << endl;

    bool foundGroundState = false;
    DataSet groundStateDataSet;
    try {
        groundStateDataSet = rootGroup.openDataSet("groundState");
        foundGroundState = true;
    } catch (GroupIException exception) {
        cout << "Did not find ground state data set, skipping initialization of coefficient matrices." << endl;
    }

    if(foundGroundState) {
        hsize_t dims2[1];
        groundStateDataSet.getSpace().getSimpleExtentDims(dims2);
        int groundStateNAtoms2 = dims2[0];

        if(nAtoms != groundStateNAtoms2) {
            cerr << "Error! The number of atoms in the ground state (nAtoms = " << groundStateNAtoms2 << ") does not match "
                 << "the number of atoms in the metadata (= " << nAtoms << ")" << endl;
            cerr << "Cannot continue" << endl;
            throw std::logic_error("Mismatching number of atoms");
        }

        AtomData *groundStateAtoms = new AtomData[groundStateNAtoms2];
        groundStateDataSet.read(groundStateAtoms, atomCompound);

        GaussianSystem groundStateSystem;
        for(int i = 0; i < groundStateNAtoms2; i++) {
            stringstream basisFile;
            basisFile << "atom_" << atomMetaData[i].type << "_basis_" << atomMetaData[i].basisName << ".tm";
            string fileName = basisFile.str();
            groundStateSystem.addCore(GaussianCore({ groundStateAtoms[i].x, groundStateAtoms[i].y, groundStateAtoms[i].z}, fileName));
        }
        if(method == "unrestricted") {
            UnrestrictedHartreeFockSolver groundStateSolver(&groundStateSystem);
            groundStateSolver.setDiisEnabled(false);
            groundStateSolver.setNIterationsMax(1e3);
            groundStateSolver.setDensityMixFactor(0.95);
            groundStateSolver.setConvergenceTreshold(1e-9);
            groundStateSolver.solve();

            coefficientMatrixUp = groundStateSolver.coeffcientMatrixUp();
            coefficientMatrixDown = groundStateSolver.coeffcientMatrixDown();
            cout << "Ground state energy: " << groundStateSolver.energy() << endl;
        } else {
            RestrictedHartreeFockSolver groundStateSolver(&groundStateSystem);
            groundStateSolver.setNIterationsMax(1e3);
            groundStateSolver.setDensityMixFactor(0.95);
            groundStateSolver.setConvergenceTreshold(1e-9);
            groundStateSolver.solve();

            coefficientMatrixUp = groundStateSolver.coefficientMatrix();
            cout << "Ground state energy: " << groundStateSolver.energy() << endl;
        }

        delete groundStateAtoms;
    }
    // Done precalculate ground state

    H5::Group statesGroup(outFile.openGroup("/states"));

    int nTotal = statesGroup.getNumObjs();
    int currentState = 0;
    for(int stateID = 0; stateID < nTotal; stateID++) {
        if(world.rank() == 0) {
            cout << "Progress: " << fixed << setprecision(2) << double(currentState) / nTotal * 100 << " %"
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
        double energy;
        if(method == "unrestricted") {
            UnrestrictedHartreeFockSolver solver(&system);
            solver.setInitialCoefficientMatrices(coefficientMatrixUp, coefficientMatrixDown);
            solver.setNIterationsMax(1e3);
            solver.setDensityMixFactor(0.95);
            solver.setConvergenceTreshold(1e-9);
            solver.solve();
            energy = solver.energy();
        } else {
            RestrictedHartreeFockSolver solver(&system);
            solver.setInitialCoefficientMatrix(coefficientMatrixUp);
            solver.setNIterationsMax(1e3);
            solver.setDensityMixFactor(0.95);
            solver.setConvergenceTreshold(1e-9);
            solver.solve();
            energy = solver.energy();
        }

        Attribute energyAttribute(stateDataSet.createAttribute("energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
        energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);
        outFile.flush(H5F_SCOPE_GLOBAL);

        currentState++;

        delete atoms;
    }
    world.barrier();
    if(world.rank() == 0) {
        cout << "Progress: 100 %, time "  << fixed << setprecision(2) << timer.elapsed() << " s" << endl;
        cout << "Done!" << endl;
    }
    outFile.close();
    return 0;
}

