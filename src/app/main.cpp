#include <iostream>
#include <iomanip>
#include <fstream>
#include <libconfig.h++>

#include <hartreefocksolver.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>

#include <H5Cpp.h>

using namespace std;
using namespace libconfig;
using namespace H5;

#define H5FILE_NAME          "SDScompound.h5"
#define DATASETNAME   "atoms"
#define LENGTH        10
#define RANK          1

int main(int argc, char* argv[])
{
    /* First structure  and dataset*/
    typedef struct AtomData {
        int    type;
        double posx;
        double posy;
        double posz;
        double partialCharge;
    } AtomData;
    AtomData       s1[LENGTH];

//    hid_t      file, dataset, space; /* Handles */
    H5File* file = new H5File( H5FILE_NAME, H5F_ACC_TRUNC );
    CompType mtype1( sizeof(AtomData) );
    mtype1.insertMember( "type", HOFFSET(AtomData, type), PredType::NATIVE_INT);
    mtype1.insertMember( "posx", HOFFSET(AtomData, posx), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "posy", HOFFSET(AtomData, posy), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "posz", HOFFSET(AtomData, posz), PredType::NATIVE_DOUBLE);

//        ofstream outFile("energies.dat");
    vec y1Range = linspace( 0.5, 4.0, 30);
    vec x2Range = linspace( 0.5, 4.0, 30);
    vec y2Range = linspace(-4.0, 4.0, 30);
    int configCounter = 0;
    for(int i = 0; i < y1Range.n_elem; i++) {
        for(int j = 0; j < x2Range.n_elem; j++) {
            for(int k = 0; k < y2Range.n_elem; k++) {
                double y1 = y1Range(i);
                double x2 = x2Range(j);
                double y2 = y2Range(k);
                cout << 0.0 << " " << 0.0 << " " << 0.0 << " "
                     << 0.0 << " " << y1 << " " << 0.0 << " "
                     << x2 << " " << y2 << " " << 0.0 << " " << endl;
                stringstream configurationName;
                configurationName << "configuration" << std::setw(9) << std::setfill('0') << configCounter;
                cout << "Creating dataset " << configurationName.str() << endl;
                try {
                    vector<GaussianCore> cores;
                    cores.push_back(GaussianCore({ 0.000, 0.000, 0.000}, "oxygen321g.tm"));
                    cores.push_back(GaussianCore({ 0.000, y1, 0.000}, "hydrogen321g.tm"));
                    cores.push_back(GaussianCore({ x2, y2, 0.000}, "hydrogen321g.tm"));
                    GaussianSystem system;
                    for(const GaussianCore &core : cores) {
                        system.addCore(core);
                    }
                    HartreeFockSolver solver(&system);
                    solver.setNIterationsMax(1e3);
                    solver.solve();
//                    outFile << solver.energy() << endl;
                    cout << solver.energy() << endl;

                    hsize_t    dim[] = {cores.size()};   /* Dataspace dimensions */
                    DataSpace space( RANK, dim );
                    DataSet dataset(file->createDataSet(configurationName.str().c_str(), mtype1, space));

                    IntType int_type(PredType::NATIVE_DOUBLE);
                    DataSpace att_space(H5S_SCALAR);
                    Attribute att = dataset.createAttribute("energy", int_type, att_space );
                    double energy = solver.energy();
                    att.write( int_type, &energy );

                    int atomIndex = 0;
                    for(const GaussianCore &core : cores) {
                        s1[atomIndex].type = int(core.atomType());
                        s1[atomIndex].posx = core.position()(0);
                        s1[atomIndex].posy = core.position()(1);
                        s1[atomIndex].posz = core.position()(2);
                        atomIndex++;
                    }
                    dataset.write( s1, mtype1 );
                } catch(std::logic_error) {
//                    outFile << "error" << endl;
                }
                configCounter++;
                file->flush(H5F_SCOPE_GLOBAL);
            }
        }
    }
    delete file;
//        outFile.close();
    return 0;
}

