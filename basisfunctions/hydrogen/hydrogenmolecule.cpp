#include "hydrogenmolecule.h"

HydrogenMolecule::HydrogenMolecule()
{
    alpha = zeros(4);
    alpha(0) = 13.00773;
    alpha(1) = 1.962079;
    alpha(2) = 0.444529;
    alpha(3) = 0.1219492;

    R = zeros(2,3);

    R(1,0) = 1.0;
}

double HydrogenMolecule::electronInteractionIntegral(int p, int r, int q, int s) {
    double A = alpha[p] + alpha[q];
    double B = alpha[r] + alpha[s];


    // TODO: FIX THESE!
    int Rp = 0;
    int Rr = 0;
    int Rq = 0;
    int Rs = 0;

    rowvec Ra = (alpha[p]*R.row(Rp) + alpha[q]*R.row(Rq))/A;
    rowvec Rb = (alpha[r]*R.row(Rr) + alpha[s]*R.row(Rs))/B;


    double t = (A*B/(A + B))*dot(Ra-Rb,Ra-Rb);

    double arg = 2*sqrt(A*B/(acos(-1)*(A+B)))*errorFunction(t)*overlapIntegral(p+Rp*4,q+Rq*4)*overlapIntegral(s+Rs*4,r+Rr*4);

    return arg;
}

double HydrogenMolecule::errorFunction(double arg){

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
    double factor = p*q/(p+q);

    // TODO: FIX THESE!
    rowvec Rp = zeros(3);
    rowvec Rq = zeros(3);

    double Rpq = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-factor*Rpq);
    double kin = 0.5*factor*(6-4*factor*Rpq)*pow(acos(-1)/(p+q),3.0/2.0)*expTerm;

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
    double factor = 1.0/(p+q);

    // TODO: FIX THESE!
    int Rp = 0;
    int Rr = 0;
    int Rq = 0;
    int Rs = 0;
    int Z = 1;

    double Rpq =dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-p*q*factor*Rpq);

    rowvec Rmc = (p*R.row(Rp) + q*R.row(Rq))*factor;


    double F0p = errorFunction(1.0/factor*dot(Rmc-R.row(0),Rmc-R.row(0)));
    double F0q = errorFunction(1.0/factor*dot(Rmc-R.row(1),Rmc-R.row(1)));

    double nucAtt = -2*Z*factor*acos(-1)*expTerm*(F0p+F0q);

    return nucAtt;
}

double HydrogenMolecule::overlapIntegral(int p, int q) {
    double factor = 1.0/(p+q);

    // TODO: FIX THESE!
    rowvec Rp = zeros(3);
    rowvec Rq = zeros(3);

    double Rpq = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-p*q*factor*Rpq);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;
}

uint HydrogenMolecule::nOrbitals()
{
    return 8;
}
