#include "hartreesolver.h"

#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;

HartreeSolver::HartreeSolver()
{
    cout << setprecision(20);
    allocateQMemory();
    reset();
}

HartreeSolver::~HartreeSolver()
{
    cleanUpQMemory();
}

double HartreeSolver::basisFunction(int index, vec position) {
    double r2 = dot(position, position);
    return exp(-alpha[index] * r2);
}

double HartreeSolver::matrixElement(int p, int r, int q, int s) {
    double denominator = (alpha[p] + alpha[q])*(alpha[r] + alpha[s])*sqrt(alpha[p] + alpha[q] + alpha[r] + alpha[s]);
    return 2 * powPi5over2 / denominator;
}

/*!
 * \brief HartreeSolver::kineticIntegral is the solution of the the integral
 * $$\langle \chi_p | -\frac{1}{2} \nabla^2 | \chi_q \rangle$$
 * \param p
 * \param q
 * \return
 */
double HartreeSolver::kineticIntegral(int p, int q) {
    double alpha_p = alpha[p];
    double alpha_q = alpha[q];
    return 3*pow(M_PI, 3.0/2.0)*alpha_q/(pow(alpha_p, 3.0/2.0)*pow(1 + alpha_q/alpha_p, 5.0/2.0));
}

/*!
 * \brief HartreeSolver::kineticIntegral is the solution of the the integral
 * $$\langle \chi_p | -\frac{1}{2} \nabla^2 | \chi_q \rangle$$
 * \param p
 * \param q
 * \return
 */
double HartreeSolver::nuclearAttractionIntegral(int p, int q) {
    double alpha_p = alpha[p];
    double alpha_q = alpha[q];
    return -4*M_PI/(alpha_p*(1 + alpha_q/alpha_p));
}

double HartreeSolver::overlapIntegral(int p, int q) {
    double alpha_p = alpha[p];
    double alpha_q = alpha[q];
    return pow(M_PI, 3.0/2.0)/(pow(alpha_p, 3.0/2.0)*pow(1 + alpha_q/alpha_p, 3.0/2.0));
}

void HartreeSolver::reset() {
    setupAlpha();
    setuph();
    setupS();
    setupQ();
    resetC();
}

void HartreeSolver::setupAlpha() {
    alpha[0] = 0.298073;
    alpha[1] = 1.242567;
    alpha[2] = 5.782948;
    alpha[3] = 38.474970;
}

void HartreeSolver::setuph() {
    h.reset();
    h = zeros(4,4);
    for(int p = 0; p < 4; p++) {
        for(int q = 0; q < 4; q++) {
            h(p,q) = kineticIntegral(p,q) + nuclearAttractionIntegral(p,q);
        }
    }
}

void HartreeSolver::setupS() {
    S.reset();
    S = zeros(4,4);
    for(int p = 0; p < 4; p++) {
        for(int q = 0; q < 4; q++) {
            S(p,q) = overlapIntegral(p, q);
        }
    }
}

void HartreeSolver::allocateQMemory() {
    if(!isQAllocated) {
        QData = new double[4*4*4*4];
        Q = new double***[4];
        for(int p = 0; p < 4; p++) {
            Q[p] = new double**[4];
            for(int r = 0; r < 4; r++) {
                Q[p][r] = new double *[4];
                for(int q = 0; q < 4; q++) {
                    Q[p][r][q] = &QData[4*4*4*p + 4*4*r + 4*q];
                    for(int s = 0; s < 4; s++) {
                        Q[p][r][q][s] = 0;
                    }
                }
            }
        }
        isQAllocated = true;
    }
}

void HartreeSolver::cleanUpQMemory() {
    if(isQAllocated) {
        for (uint i = 0; i < 4; ++i) {
            for (uint j = 0; j < 4; ++j){
                delete [] Q[i][j];
            }
            delete [] Q[i];
        }
        delete [] Q;
        delete []QData;
    }
}

void HartreeSolver::setupQ() {
    for(int p = 0; p < 4; p++) {
        for(int r = 0; r < 4; r++) {
            for(int q = 0; q < 4; q++) {
                for(int s = 0; s < 4; s++) {
                    Q[p][r][q][s] = matrixElement(p, r, q, s);
                }
            }
        }
    }
}

void HartreeSolver::resetC() {
    C.reset();
    C = ones(4);
}

void HartreeSolver::advance() {
    normalizeCwithRegardsToS();
    setupF();

    vec s;
    mat U;
    eig_sym(s, U, S);

    mat V = U*diagmat(1.0/sqrt(s));

    F = V.t() * F * V;

    vec eps;
    mat Cmat;
    eig_sym(eps, Cmat, F);

    C = V*Cmat.col(0);
    normalizeCwithRegardsToS();

    double energy = 0;

    for(int p = 0; p < 4; p++) {
        for(int q = 0; q < 4; q++) {
            energy += 2 * C(p) * C(q) * h(p,q);
        }
    }

    for(int p = 0; p < 4; p++) {
        for(int q = 0; q < 4; q++) {
            for(int r = 0; r < 4; r++) {
                for(int s = 0; s < 4; s++) {
                    energy += Q[p][r][s][q] * C(p) * C(q) * C(r) * C(s);
                }
            }
        }
    }
    cout << energy << endl;
}

void HartreeSolver::normalizeCwithRegardsToS(){
    double factor = 0.0;

    for(uint i= 0; i < C.n_elem; i++){
        for(uint j= 0; j < C.n_elem; j++){
            factor += C(i)*S(i,j)*C(j);
        }
    }

    C = C/sqrt(factor);
}

void HartreeSolver::setupF() {
    F = zeros(4,4);
    for(int p = 0; p < 4; p++) {
        for(int q = 0; q < 4; q++) {
            F(p,q) = h(p,q);
            for(int r = 0; r < 4; r++) {
                for(int s = 0; s < 4; s++) {
                    F(p,q) += Q[p][r][q][s] * C(r) * C(s);
                }
            }
        }
    }
}
