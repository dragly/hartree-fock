#include <iostream>
#include <unittest++/UnitTest++.h>
#include <src/hartreesolver.h>
#include <src/hartreefocksolver.h>
#include <src/basisfunctions/helium/heliumhartree.h>
#include <src/basisfunctions/hydrogen/hydrogenmolecule.h>

#include <fstream>

using namespace std;

//TEST(MatrixElement) {
//    Helium basisFunction;
//    int p = 0;
//    int q = 1;
//    int r = 2;
//    int s = 3;

//    double value = basisFunction.electronInteractionIntegral(p,r,q,s);

//    CHECK_CLOSE(0.075820496873881579, value, 1e-9);
//}

//TEST(KineticIntegral) {
//    Helium basisFunction;
//    int p = 0;
//    int q = 3;

//    double value = basisFunction.kineticIntegral(p,q);

//    CHECK_CLOSE(0.0204654960291562, value, 1e-9);
//}

//TEST(NuclearAttractionIntegral) {
//    Helium basisFunction;
//    int p = 1;
//    int q = 2;

//    double value = basisFunction.nuclearAttractionIntegral(p,q);

//    CHECK_CLOSE(-1.78867607774792, value, 1e-9);
//}

//TEST(OverlapIntegral) {
//    Helium basisFunction;
//    int p = 3;
//    int q = 2;

//    double value = basisFunction.overlapIntegral(p, q);

//    CHECK_CLOSE(0.01891203832678, value, 1e-9);
//}

//TEST(HeliumAdvanceMany) {
//    Helium basisFunction;
//    HartreeSolver solver(&basisFunction);
//    solver.reset();
//    for(int i = 0; i < 500; i++) {
//        solver.advance();
//    }
//    double value = solver.energy();
//    CHECK_CLOSE(-2.855160382370257377, value, 1e-7);
//}

//TEST(HydrogenOverlapIntegral)  {
//    HydrogenMolecule basisFunction;
//    double value = basisFunction.overlapIntegral(6,2);
//    cout << "value: " << value << endl;
//}

//TEST(HydrogenElectronInteractionIntegral)  {
//    HydrogenMolecule basisFunction;
//    double value = basisFunction.electronInteractionIntegral(1,5,7,2);
//    cout << "value: " << value << endl;
//}

//TEST(HydrogenNuclearAttractionIntegral)  {
//    HydrogenMolecule basisFunction;
//    double value = basisFunction.nuclearAttractionIntegral(3,7);
//    cout << "value: " << value << endl;
//}

TEST(HydrogenAdvanceMany) {
    HydrogenMolecule basisFunction(1.0);
    HartreeSolver solver(&basisFunction);
    for(int i = 0; i < 100; i++) {
        solver.advance();
    }
    cout << "Energy: " << solver.energy() << endl;
    cout << "Energy with repulsion: " << solver.energy() + basisFunction.nuclearRepulsion() << endl;
}

TEST(HydrogenAdvanceManyHF) {
    cout << "Hydrogen with HF:" << endl;
    HydrogenMolecule basisFunction(1.0);
    HartreeFockSolver solver(&basisFunction);
    for(int i = 0; i < 100; i++) {
        solver.advance();
    }
    cout << "Energy: " << solver.energy() << endl;
    cout << "Energy with repulsion: " << solver.energy() + basisFunction.nuclearRepulsion() << endl;
}

//TEST(HydrogenPlot) {
//    vec distances = linspace(1.0, 4.0, 100);
//    std::ofstream energyFile;
//    energyFile.open("energies.dat");
//    for(double distance : distances) {
//        cout << distance << endl;
//        HydrogenMolecule basisFunction(distance);
//        HartreeSolver solver(&basisFunction);
//        solver.reset();
//        for(int i = 0; i < 100; i++) {
//            solver.advance();
//        }
//        double energy = solver.energy() + basisFunction.nuclearRepulsion();
//        energyFile << distance << " " << energy << "\n";
//    }
//    energyFile.close();
//}

int main()
{
    UnitTest::RunAllTests();
    return 0;
}

