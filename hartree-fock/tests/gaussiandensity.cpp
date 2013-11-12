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
//        HartreeFockSolver solver(&system);
//        for(int i = 0; i < 100; i++) {
//            solver.advance();
//        }
//        cout << "Energy: " << solver.energy() << endl;
//        mat C = solver.coefficientMatrix();
//        C.save("coefficients.dat");
//        mat S = solver.overlapMatrix();
//        S.save("overlap.dat");
        mat C;
        mat S;
        C.load("coefficients.dat");
        S.load("overlap.dat");
        int nCores = 2;
        int nBF = system.nBasisFunctions();
        int nP = system.nParticles();
        mat bC = zeros(nBF,nP / 2);
        bC.submat(span(0,nBF - 1),span(0,nP / 2 / nCores - 1)) = eye(nBF, nP / 2 / nCores);
        bC.submat(span(nBF / 2,nBF - 1), span(nP / 2 / nCores, nP / 2 - 1)) = eye(nBF / 2, nP / 2 / nCores);
        for(uint k = 0; k < nP / 2; k++) {
            double factor = 0.0;
            for(uint p = 0; p < nBF; p++){
                for(uint q = 0; q < nBF; q++){
                    factor += bC(p,k) * S(p,q) * bC(q,k);
                }
            }
            bC.col(k) = bC.col(k) / sqrt(factor);
        }
        cout << "C = \n" << C << endl;
        cout << "bC = \n" << bC << endl;
//        mat SC = S * C;
        vec x = linspace(-2, 2, 50);
        vec y = linspace(-2, 2, 50);
        vec z = linspace(-2, 2, 50);
        double dx = x(1) - x(0);
        double dy = y(1) - y(0);
        double dz = z(1) - z(0);
        double basisDensitySum = 0;
        double densitySum = 0;
        ofstream dataFile("density-oxygen.dat");
        ofstream basisDataFile("density-oxygen-basis.dat");
        ofstream diffDataFile("density-oxygen-diff.dat");
        for(uint i = 0; i < x.n_elem; i++) {
            for(uint j = 0; j < y.n_elem; j++) {
                for(uint k = 0; k < z.n_elem; k++) {
                    double density = system.particleDensity(C, x(i), y(j), z(k));
                    double basisDensity = system.particleDensity(bC, x(i), y(j), z(k)) / nCores;
                    double densityDiff = density - basisDensity;
                    densitySum += density * dx * dy * dz;
                    basisDensitySum += basisDensity * dx * dy * dz;
                    if(k == z.n_elem / 2) {
                        dataFile << density << " ";
                        basisDataFile << basisDensity << " ";
                        diffDataFile << densityDiff << " ";
                    }
                }
            }
            dataFile << endl;
            basisDataFile << endl;
            diffDataFile << endl;
        }
        dataFile.close();
        diffDataFile.close();
        basisDataFile.close();
        cout << densitySum << endl;
        cout << basisDensitySum << endl;
    }
//    TEST(OrthogonalityCheck) {
//        mat C;
//        mat S;
//        C.load("coefficients.dat");
//        S.load("overlap.dat");
//        mat SC = S * C;
//        mat CSC = C.t() * SC;
//        mat Ct = C.t();

//        mat CSC2 = zeros(8,8);
//        for(int i = 0; i < 8; i++) {
//            for(int j = 0; j < 8; j++) {
//                double result = 0;
//                for(int k = 0; k < 18; k++) {
//                    for(int l = 0; l < 18; l++) {
//                        result += Ct(i,k) * C(l,j) * S(k,l);
//                    }
//                }
//                CSC2(i,j) = result;
//            }
//        }
////        for(int i = 0; i < C.n_cols; i++) {
////            for(int j = 0; j < C.n_cols; j++) {
////                cout << "C(" << i << ") * C(" << j << ") = " << dot(C.col(i), C.col(j)) << endl;
////            }
////        }
//        cout << CSC << endl;
//        cout << CSC2 << endl;
//    }
}
