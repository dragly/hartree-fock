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
    //    TEST(OxygenDensity) {
    //        GaussianOxygen431G system;
    ////        HartreeFockSolver solver(&system);
    ////        for(int i = 0; i < 100; i++) {
    ////            solver.advance();
    ////        }
    ////        cout << "Energy: " << solver.energy() << endl;
    ////        mat C = solver.coefficientMatrix();
    ////        C.save("coefficients.dat");
    ////        mat S = solver.overlapMatrix();
    ////        S.save("overlap.dat");
    //        mat C;
    //        mat S;
    //        C.load("coefficients.dat");
    //        S.load("overlap.dat");
    //        int nCores = 2;
    //        int nBF = system.nBasisFunctions();
    //        int nP = system.nParticles();
    //        mat bC = zeros(nBF,nP / 2);
    //        bC.submat(span(0,nBF - 1),span(0,nP / 2 / nCores - 1)) = eye(nBF, nP / 2 / nCores);
    //        bC.submat(span(nBF / 2,nBF - 1), span(nP / 2 / nCores, nP / 2 - 1)) = eye(nBF / 2, nP / 2 / nCores);
    //        for(uint k = 0; k < nP / 2; k++) {
    //            double factor = 0.0;
    //            for(uint p = 0; p < nBF; p++){
    //                for(uint q = 0; q < nBF; q++){
    //                    factor += bC(p,k) * S(p,q) * bC(q,k);
    //                }
    //            }
    //            bC.col(k) = bC.col(k) / sqrt(factor);
    //        }
    //        cout << "C = \n" << C << endl;
    //        cout << "bC = \n" << bC << endl;
    ////        mat SC = S * C;
    //        vec x = linspace(-2, 2, 50);
    //        vec y = linspace(-2, 2, 50);
    //        vec z = linspace(-2, 2, 50);
    //        double dx = x(1) - x(0);
    //        double dy = y(1) - y(0);
    //        double dz = z(1) - z(0);
    //        double basisDensitySum = 0;
    //        double densitySum = 0;
    //        ofstream dataFile("density-oxygen.dat");
    //        ofstream basisDataFile("density-oxygen-basis.dat");
    //        ofstream diffDataFile("density-oxygen-diff.dat");
    //        for(uint i = 0; i < x.n_elem; i++) {
    //            for(uint j = 0; j < y.n_elem; j++) {
    //                for(uint k = 0; k < z.n_elem; k++) {
    //                    double density = system.particleDensity(C, x(i), y(j), z(k));
    //                    double basisDensity = system.particleDensity(bC, x(i), y(j), z(k)) / nCores;
    //                    double densityDiff = density - basisDensity;
    //                    densitySum += density * dx * dy * dz;
    //                    basisDensitySum += basisDensity * dx * dy * dz;
    //                    if(k == z.n_elem / 2) {
    //                        dataFile << density << " ";
    //                        basisDataFile << basisDensity << " ";
    //                        diffDataFile << densityDiff << " ";
    //                    }
    //                }
    //            }
    //            dataFile << endl;
    //            basisDataFile << endl;
    //            diffDataFile << endl;
    //        }
    //        dataFile.close();
    //        diffDataFile.close();
    //        basisDataFile.close();
    //        cout << densitySum << endl;
    //        cout << basisDensitySum << endl;
    //    }

    void densityToFile(string fileName, const GaussianSystem& system, const mat& C) {
        vec x = linspace(-3, 3, 200);
        vec y = linspace(-3, 3, 200);
        //        vec z = linspace(-3, 3, 3);
//        double dx = x(1) - x(0);
//        double dy = y(1) - y(0);
        //        double dz = z(1) - z(0);
        //        double basisDensitySum = 0;
//        double densitySum = 0;
        ofstream dataFile(fileName);
        //        ofstream basisDataFile("density-oxygen-basis.dat");
        //        ofstream diffDataFile("density-oxygen-diff.dat");
        for(uint i = 0; i < x.n_elem; i++) {
//            cout << "Calculating density for x = " << x(i) << endl;
            for(uint j = 0; j < y.n_elem; j++) {
                //                for(uint k = 0; k < z.n_elem; k++) {
                double density = system.particleDensity(C, x(i), y(j), 0.0);
                //                    double basisDensity = system.particleDensity(bC, x(i), y(j), z(k)) / nCores;
                //                    double densityDiff = density - basisDensity;
                //                    densitySum += density * dx * dy * dz;
                //                    basisDensitySum += basisDensity * dx * dy * dz;
                //                    if(k == z.n_elem / 2) {
                dataFile << density << " ";
                //                        basisDataFile << basisDensity << " ";
                //                        diffDataFile << densityDiff << " ";
                //                    }
                //                }
            }
            dataFile << endl;
            //            basisDataFile << endl;
            //            diffDataFile << endl;
        }
        dataFile.close();
        //        diffDataFile.close();
        //        basisDataFile.close();
//        cout << "Density sum: " << densitySum << endl;
    }

    TEST(GaussianCore) {
        vector<GaussianCore> cores;
        cores.push_back(GaussianCore({0,0,0}, "oxygen431g.tm"));
        cores.push_back(GaussianCore({-1.43,1.108,0}, "hydrogen431g.tm"));
        cores.push_back(GaussianCore({1.43,1.108,0}, "hydrogen431g.tm"));
        GaussianSystem system;
        for(const GaussianCore &core : cores) {
            system.addCore(core);
        }
        mat C;
        HartreeFockSolver solver(&system);
        for(int i = 0; i < 100; i++) {
            solver.advance();
        }
        cout << "Energy: " << solver.energy() << endl;
        C = solver.coefficientMatrix();
        C.save("coefficients-water.dat");
        densityToFile("density-water.dat", system, C);

        int counter = 0;
        for(const GaussianCore &core : cores) {
//            cout << "Subsystem: " << counter << endl;
            GaussianSystem system2;
            system2.addCore(core);
            mat C2;
            HartreeFockSolver solver2(&system2);
            for(int i = 0; i < 100; i++) {
                solver2.advance();
            }
//            cout << "Energy: " << solver2.energy() << endl;
            C2 = solver2.coefficientMatrix();
            stringstream coefficientsFileName;
            coefficientsFileName << "coefficients-water" << counter << ".dat";
            C2.save(coefficientsFileName.str(), raw_ascii);
            stringstream densityFileName;
            densityFileName << "density-water" << counter << ".dat";
            densityToFile(densityFileName.str(), system2, C2);
            counter++;
        }
    }

    //    TEST(GaussianCore2) {
    //        GaussianCore core1({0,0,0}, "oxygen431g.tm");
    //        GaussianCore core2({-1.43,1.108,0}, "hydrogen431g.tm");
    //        GaussianCore core3({1.43,1.108,0}, "hydrogen431g.tm");
    //        //        GaussianCore core4({2.86,2.216,0}, "oxygen431g.tm");
    //        //        GaussianCore core5({2.86,4.026,0}, "hydrogen431g.tm");
    //        //        GaussianCore core6({4.29,1.108,0}, "hydrogen431g.tm");
    //        GaussianCore core4({0,-1.0,0}, "oxygen431g.tm");
    //        GaussianCore core5({-1.43,-2.108,0}, "hydrogen431g.tm");
    //        GaussianCore core6({1.43,-2.108,0}, "hydrogen431g.tm");
    //        GaussianSystem system;
    //        system.addCore(core1);
    //        system.addCore(core2);
    //        system.addCore(core3);
    //        system.addCore(core4);
    //        system.addCore(core5);
    //        system.addCore(core6);
    //        bool recalc = false;
    //        mat C;
    //        if(recalc) {
    //            HartreeFockSolver solver(&system);
    //            for(int i = 0; i < 100; i++) {
    //                solver.advance();
    //            }
    //            cout << "Energy: " << solver.energy() << endl;
    //            mat C = solver.coefficientMatrix();
    //            C.save("coefficients-water2.dat");
    //        } else {
    //            C.load("coefficients-water2.dat");
    //        }
    //        vec x = linspace(-5, 5, 75);
    //        vec y = linspace(-5, 5, 75);
    //        vec z = linspace(-5, 5, 75);
    //        double dx = x(1) - x(0);
    //        double dy = y(1) - y(0);
    //        double dz = z(1) - z(0);
    //        //        double basisDensitySum = 0;
    //        double densitySum = 0;
    //        ofstream dataFile("density-water2.dat");
    //        //        ofstream basisDataFile("density-oxygen-basis.dat");
    //        //        ofstream diffDataFile("density-oxygen-diff.dat");
    //        for(uint i = 0; i < x.n_elem; i++) {
    //            cout << "Density calc for x = " << x(i) << endl;
    //            for(uint j = 0; j < y.n_elem; j++) {
    //                for(uint k = 0; k < z.n_elem; k++) {
    //                    double density = system.particleDensity(C, x(i), y(j), z(k));
    //                    //                    double basisDensity = system.particleDensity(bC, x(i), y(j), z(k)) / nCores;
    //                    //                    double densityDiff = density - basisDensity;
    //                    densitySum += density * dx * dy * dz;
    //                    //                    basisDensitySum += basisDensity * dx * dy * dz;
    //                    if(k == z.n_elem / 2) {
    //                        dataFile << density << " ";
    //                        //                        basisDataFile << basisDensity << " ";
    //                        //                        diffDataFile << densityDiff << " ";
    //                    }
    //                }
    //            }
    //            dataFile << endl;
    //            //            basisDataFile << endl;
    //            //            diffDataFile << endl;
    //        }
    //        dataFile.close();
    //        //        diffDataFile.close();
    //        //        basisDataFile.close();
    //        cout << "density sum: " << densitySum << endl;
    //        //        cout << basisDensitySum << endl;
    //    }
}
