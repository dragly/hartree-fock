#include <electronsystems/gaussian/simplegaussian.h>
#include <electronsystems/gaussian/gaussianoxygen431g.h>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <hartreefocksolver.h>
#include <unittest++/UnitTest++.h>

using namespace std;
using namespace arma;

SUITE(GaussianDensity) {
//    TEST(SimpleDensity) {
//        SimpleGaussian system;
//        mat C = {0.5,0.0,
//                 0.0,0.5,
//                 0.5,0.0,
//                 0.0,0.5};
//        C.reshape(2,4);
//        C = C.t();
//        cout << system.particleDensity(C, 0.0, 0.0, 0.0) << endl;
////        vec x = linspace(-5, 5, 100);
////        vec y = linspace(-5, 5, 100);
////        ofstream dataFile("density.dat");
////        for(uint i = 0; i < x.n_elem; i++) {
////            for(uint j = 0; j < y.n_elem; j++) {
////                dataFile << system.particleDensity(C, x(i), y(j), 0.0) << " ";
////            }
////            dataFile << endl;
////        }
////        dataFile.close();
//    }
    TEST(SimpleDensity2) {
        SimpleGaussian system;
        mat C = {0.25,0.0,
                 0.25,0.5,
                 0.25,0.0,
                 0.25,0.5};
        C.reshape(2,4);
        C = C.t();
        vec x = linspace(-5, 5, 100);
        vec y = linspace(-5, 5, 100);
        ofstream dataFile("density.dat");
        for(uint i = 0; i < x.n_elem; i++) {
            for(uint j = 0; j < y.n_elem; j++) {
                dataFile << system.particleDensity(C, x(i), y(j), 0.0) << " ";
            }
            dataFile << endl;
        }
        dataFile.close();
    }
    TEST(OxygenDensity) {
        GaussianOxygen431G system;
        HartreeFockSolver solver(&system);
//        for(int i = 0; i < 100; i++) {
//            solver.advance();
//        }
//        cout << "Energy: " << solver.energy() << endl;
//        mat C = solver.coefficientMatrix();
//        C.save("coefficients.dat");
        mat C;
        C.load("coefficients.dat");
        vec x = linspace(-2, 2, 100);
        vec y = linspace(-2, 2, 100);
        double densitySum = 0;
        ofstream dataFile("density-oxygen.dat");
        for(uint i = 0; i < x.n_elem; i++) {
            for(uint j = 0; j < y.n_elem; j++) {
                double density = system.particleDensity(C, x(i), y(j), 0.0);
                dataFile << density << " ";
                densitySum += density;
            }
            dataFile << endl;
        }
        dataFile.close();
        cout << densitySum << endl;
    }
}
