#include "hartreefocksolver.h"

#include <electronsystems/electronsystem.h>

#include <armadillo>
#include <iomanip>
#include <cmath>

using namespace arma;
using namespace std;

HartreeFockSolver::HartreeFockSolver(ElectronSystem *basisFunction) :
    m_electronSystem(basisFunction)
{
    cout << setprecision(20);
    allocateQMemory();
    reset();
}

HartreeFockSolver::~HartreeFockSolver()
{
    cleanUpQMemory();
}

void HartreeFockSolver::reset() {
    setuph();
    setupS();
    setupQ();
    resetC();
//    normalizeCwithRegardsToS();
    setupP();
}

void HartreeFockSolver::setuph() {
//    cout << "Setting up h" << endl;
    ElectronSystem* f = m_electronSystem;
    uint nOrbitals = f->nBasisFunctions();
    h.reset();
    h = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            //            h(p,q) = f->kineticIntegral(p,q) + f->nuclearAttractionIntegral(p,q);
            h(p,q) = f->uncoupledIntegral(p,q);
        }
    }
//    cout << h << endl;
}

void HartreeFockSolver::setupS() {
//    cout << "Setting up S" << endl;
    ElectronSystem* f = m_electronSystem;
    uint nOrbitals = f->nBasisFunctions();
    S.reset();
    S = zeros(nOrbitals,nOrbitals);
    for(uint p = 0; p < nOrbitals; p++) {
        for(uint q = 0; q < nOrbitals; q++) {
            S(p,q) = f->overlapIntegral(p, q);
        }
    }
//    cout << S << endl;
}

void HartreeFockSolver::allocateQMemory() {
    if(!isQAllocated) {
//        cout << "Allocating memory for Q matrix" << endl;
        ElectronSystem* f = m_electronSystem;
        uint n = f->nBasisFunctions();
        Q.set_size(n,n);
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                Q(i,j) = zeros(n,n);
            }
        }
        //        QData = new double[n*n*n*n];
        //        Q = new double***[n];
        //        for(uint p = 0; p < n; p++) {
        //            Q[p] = new double**[n];
        //            for(uint r = 0; r < n; r++) {
        //                Q[p][r] = new double *[n];
        //                for(uint q = 0; q < n; q++) {
        //                    Q[p][r][q] = &QData[n*n*n*p + n*n*r + n*q];
        //                    for(uint s = 0; s < n; s++) {
        //                        Q[p][r][q][s] = 0;
        //                    }
        //                }
        //            }
        //        }
        isQAllocated = true;
//        cout << "Done allocating memory" << endl;
    }
}

void HartreeFockSolver::cleanUpQMemory() {
    if(isQAllocated) {
        Q.reset();
        //        ElectronSystem* f = m_electronSystem;
        //        uint n = f->nOrbitals();
        //        for (uint i = 0; i < n; ++i) {
        //            for (uint j = 0; j < n; ++j){
        //                delete [] Q[i][j];
        //            }
        //            delete [] Q[i];
        //        }
        //        delete [] Q;
        //        delete []QData;
    }
}

void HartreeFockSolver::setupQ() {
//    cout << "Setting up Q" << endl;
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
//    cout << "n=" << n <<endl;
//    fstream qFile;
//    qFile.open("Q2.dat", ios::out);
//    qFile << setprecision(9);
    for(uint p = 0; p < n; p++) {
        cout << "Q calculating for p " << p << " of " << n << endl;
        for(uint r = 0; r < n; r++) {
            for(uint q = 0; q < n; q++) {
                for(uint s = 0; s < n; s++) {
                    Q(p,r)(q,s) = f->coupledIntegral(p, q, r, s); // NOTE: Indexes changed on purpose in element and call
//                    qFile << Q(p,r)(q,s) << "\n";
                }
            }
        }
    }
//    qFile.close();
}

void HartreeFockSolver::resetC() {
//    cout << "Resetting C" << endl;
    ElectronSystem* f = m_electronSystem;
    C.reset();
    C = zeros(f->nBasisFunctions(), f->nParticles() / 2);
}

void HartreeFockSolver::advance() {
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;
    nk = max(int(nk), 2); // TODO Fix this hack for the Hydrogen atom
    setupF();

    vec s;
    mat U;
    eig_sym(s, U, S);

    mat V = U*diagmat(1.0/sqrt(s));

    F = V.t() * F * V;

    vec eps;
    mat Cmat;
    eig_sym(eps, Cmat, F);


    C = V*Cmat.submat(0, 0, no - 1, nk - 1);
    normalizeCwithRegardsToS();

    setupP();

    double energy = 0;

    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            energy += P(p,q) * h(p,q);
        }
    }

    for(uint p = 0; p < no; p++) {
        for(uint q = 0; q < no; q++) {
            for(uint r = 0; r < no; r++) {
                for(uint s = 0; s < no; s++) {
                    energy += 0.25 * Qtilde(p, q, r, s) * P(p,q) * P(s,r);
                }
            }
        }
    }
    energy += m_electronSystem->additionalEnergyTerms();
    m_energy = energy;
//    cout << m_energy << endl;
}

void HartreeFockSolver::normalizeCwithRegardsToS(){
//    cout << "Normalizing C" << endl;
    ElectronSystem* f = m_electronSystem;
    uint no = f->nBasisFunctions();
    uint nk = f->nParticles() / 2;

    for(uint k = 0; k < nk; k++) {
        double factor = 0.0;
        for(uint p = 0; p < no; p++){
            for(uint q = 0; q < no; q++){
                factor += C(p,k) * S(p,q) * C(q,k);
            }
        }
        C.col(k) = C.col(k) / sqrt(factor);
    }
}

double HartreeFockSolver::Qtilde(int p, int q, int r, int s) {
    return 2 * Q(p,r)(q,s) - Q(p,r)(s,q);
//    return 2 * Q[p][r][q][s] - Q[p][r][s][q];
}

void HartreeFockSolver::setupF() {
    ElectronSystem* f = m_electronSystem;
    uint n = f->nBasisFunctions();
    F = zeros(n,n);
    for(uint p = 0; p < n; p++) {
        for(uint q = 0; q < n; q++) {
            F(p,q) = h(p,q);
            for(uint r = 0; r < n; r++) {
                for(uint s = 0; s < n; s++) {
                    F(p,q) += 0.5 * Qtilde(p, q, r, s) * P(s,r);
                }
            }
        }
    }
}

void HartreeFockSolver::setupP() {
//    cout << "Setting up P" << endl;
    P = 2 * C * C.t();
}
