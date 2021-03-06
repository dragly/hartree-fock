#include "hydrogenmolecule.h"

#include <fstream>

HydrogenMolecule::HydrogenMolecule(double distance)
{

    alpha = zeros<rowvec>(4);
    alpha(0) = 13.00773;
    alpha(1) = 1.962079;
    alpha(2) = 0.444529;
    alpha(3) = 0.1219492;

    R = zeros(2,3);

    R(1,0) = distance;
}

HydrogenMolecule::~HydrogenMolecule()
{
}

double HydrogenMolecule::coupledIntegral(int p, int q, int r, int s) {
    int pIndex = p % nOrbitalsPerNuclei;
    int qIndex = q % nOrbitalsPerNuclei;
    int rIndex = r % nOrbitalsPerNuclei;
    int sIndex = s % nOrbitalsPerNuclei;

    int Rp = p / nOrbitalsPerNuclei;
    int Rq = q / nOrbitalsPerNuclei;
    int Rr = r / nOrbitalsPerNuclei;
    int Rs = s / nOrbitalsPerNuclei;

    double A = alpha[pIndex] + alpha[qIndex];
    double B = alpha[rIndex] + alpha[sIndex];

    rowvec Ra = (alpha[pIndex]*R.row(Rp) + alpha[qIndex]*R.row(Rq))/A;
    rowvec Rb = (alpha[rIndex]*R.row(Rr) + alpha[sIndex]*R.row(Rs))/B;


    double t = (A*B/(A + B))*dot(Ra-Rb,Ra-Rb);

    double arg = 2*sqrt(A*B/(acos(-1)*(A+B)))*boysFunction(t)*overlapIntegral(p,q)*overlapIntegral(s,r);
    return arg;
}

double HydrogenMolecule::additionalEnergyTerms() {
    return 1/sqrt(dot(R.row(0) - R.row(1),R.row(0)- R.row(1)));
}

double HydrogenMolecule::boysFunction(double arg){

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
double HydrogenMolecule::kineticIntegral(int p, int q) {
    int pIndex = p % nOrbitalsPerNuclei;
    int qIndex = q % nOrbitalsPerNuclei;

    double ap = alpha(pIndex);
    double aq = alpha(qIndex);

    double factor = ap*aq/(ap+aq);

    int RpIndex = p / nOrbitalsPerNuclei;
    int RqIndex = q / nOrbitalsPerNuclei;

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
double HydrogenMolecule::nuclearAttractionIntegral(int p, int q) {
    int pIndex = p % nOrbitalsPerNuclei;
    int qIndex = q % nOrbitalsPerNuclei;

    double ap = alpha(pIndex);
    double aq = alpha(qIndex);

    double factor = 1.0/(ap+aq);

    int Rp = p / nOrbitalsPerNuclei;
    int Rq = q / nOrbitalsPerNuclei;
    int Z = 1;

    double Rpq =dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-ap*aq*factor*Rpq);

    rowvec Rmc = (ap*R.row(Rp) + aq*R.row(Rq))*factor;


    double F0p = boysFunction(1.0/factor*dot(Rmc-R.row(0),Rmc-R.row(0)));
    double F0q = boysFunction(1.0/factor*dot(Rmc-R.row(1),Rmc-R.row(1)));

    double nucAtt = -2*Z*factor*acos(-1)*expTerm*(F0p+F0q);

    return nucAtt;
}

double HydrogenMolecule::overlapIntegral(int p, int q) {
    int pIndex = p % nOrbitalsPerNuclei;
    int qIndex = q % nOrbitalsPerNuclei;

    double ap = alpha(pIndex);
    double aq = alpha(qIndex);

    double factor = 1.0/(ap+aq);

    int RpIndex = p / nOrbitalsPerNuclei;
    int RqIndex = q / nOrbitalsPerNuclei;

    rowvec Rp = R.row(RpIndex);
    rowvec Rq = R.row(RqIndex);

    double Rpq = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-ap*aq*factor*Rpq);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;
}


double HydrogenMolecule::uncoupledIntegral(int p, int q)
{
    return kineticIntegral(p,q) + nuclearAttractionIntegral(p,q);
}
