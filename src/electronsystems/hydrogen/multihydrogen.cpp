#include "multihydrogen.h"

#include <fstream>

MultiHydrogen::MultiHydrogen(mat nucleiPositions) :
    R(nucleiPositions),
    m_nNuclei(nucleiPositions.n_rows),
    m_nOrbitalsPerNuclei(4)
{
    alpha = zeros(m_nOrbitalsPerNuclei);
    alpha(0) = 13.00773;
    alpha(1) = 1.962079;
    alpha(2) = 0.444529;
    alpha(3) = 0.1219492;
}

MultiHydrogen::~MultiHydrogen()
{
}

double MultiHydrogen::electronInteractionIntegral(int p, int r, int q, int s) {
    int pIndex = p % m_nOrbitalsPerNuclei;
    int qIndex = q % m_nOrbitalsPerNuclei;
    int rIndex = r % m_nOrbitalsPerNuclei;
    int sIndex = s % m_nOrbitalsPerNuclei;

    int Rp = p / m_nOrbitalsPerNuclei;
    int Rq = q / m_nOrbitalsPerNuclei;
    int Rr = r / m_nOrbitalsPerNuclei;
    int Rs = s / m_nOrbitalsPerNuclei;

    double A = alpha[pIndex] + alpha[qIndex];
    double B = alpha[rIndex] + alpha[sIndex];

    rowvec Ra = (alpha[pIndex]*R.row(Rp) + alpha[qIndex]*R.row(Rq))/A;
    rowvec Rb = (alpha[rIndex]*R.row(Rr) + alpha[sIndex]*R.row(Rs))/B;


    double t = (A*B/(A + B))*dot(Ra-Rb,Ra-Rb);

    double arg = 2*sqrt(A*B/(acos(-1)*(A+B)))*errorFunction(t)*overlapIntegral(p,q)*overlapIntegral(s,r);
    return arg;
}

double MultiHydrogen::additionalEnergyTerms() {
    double repulsionEnergy = 0;
    for(uint i = 0; i < m_nNuclei; i++) {
        for(uint j = i + 1; j < m_nNuclei; j++) {
            rowvec Rij = R.row(i) - R.row(j);
            repulsionEnergy += 1/sqrt(dot(Rij, Rij));
        }
    }
    return repulsionEnergy;
}

double MultiHydrogen::errorFunction(double arg){

    if (arg < 1.0E-6){
        return 1.0;
    }

    else{
        arg = sqrt(arg);
        double f = 1.0/arg * erf(arg) *sqrt(acos(-1))/2.0;
        return f;
    }

}

/*!
 * \brief HartreeSolver::kineticIntegral is the solution of the the integral
 * $$\langle \chi_p | -\frac{1}{2} \nabla^2 | \chi_q \rangle$$
 * \param p
 * \param q
 * \return
 */
double MultiHydrogen::kineticIntegral(int p, int q) {
    int pIndex = p % m_nOrbitalsPerNuclei;
    int qIndex = q % m_nOrbitalsPerNuclei;

    double ap = alpha(pIndex);
    double aq = alpha(qIndex);

    double factor = ap*aq/(ap+aq);

    int RpIndex = p / m_nOrbitalsPerNuclei;
    int RqIndex = q / m_nOrbitalsPerNuclei;

    rowvec Rp = R.row(RpIndex);
    rowvec Rq = R.row(RqIndex);

    double Rpq = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-factor*Rpq);
    double kin = 0.5*factor*(6-4*factor*Rpq)*pow(acos(-1)/(ap+aq),3.0/2.0)*expTerm;

    return kin;
}

/*!
 * \brief HartreeSolver::kineticIntegral is the solution of the the integral
 * $$\langle \chi_p | -\frac{1}{2} \nabla^2 | \chi_q \rangle$$
 * \param p
 * \param q
 * \return
 */
double MultiHydrogen::nuclearAttractionIntegral(int p, int q) {
    int pIndex = p % m_nOrbitalsPerNuclei;
    int qIndex = q % m_nOrbitalsPerNuclei;

    double ap = alpha(pIndex);
    double aq = alpha(qIndex);

    double apPlusAq = ap + aq;
    double oneOverApPlusAq = 1.0/(apPlusAq);

    int RpIndex = p / m_nOrbitalsPerNuclei;
    int RqIndex = q / m_nOrbitalsPerNuclei;

    rowvec Rp = R.row(RpIndex);
    rowvec Rq = R.row(RqIndex);

    int Z = 1;

    double Rpq =dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-ap*aq*oneOverApPlusAq*Rpq);

    rowvec RP = (ap*Rp + aq*Rq)*oneOverApPlusAq;

    double errorFunctionSum = 0;

    for(int i = 0; i < m_nNuclei; i++) {
        errorFunctionSum += errorFunction(apPlusAq*dot(RP-R.row(i),RP-R.row(i)));
    }

    double nucAtt = -2*Z*oneOverApPlusAq*acos(-1)*expTerm*errorFunctionSum;

    return nucAtt;
}

double MultiHydrogen::overlapIntegral(int p, int q) {
    int pIndex = p % m_nOrbitalsPerNuclei;
    int qIndex = q % m_nOrbitalsPerNuclei;

    double ap = alpha(pIndex);
    double aq = alpha(qIndex);

    double factor = 1.0/(ap+aq);

    int RpIndex = p / m_nOrbitalsPerNuclei;
    int RqIndex = q / m_nOrbitalsPerNuclei;

    rowvec Rp = R.row(RpIndex);
    rowvec Rq = R.row(RqIndex);

    double Rpq = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-ap*aq*factor*Rpq);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;
}
