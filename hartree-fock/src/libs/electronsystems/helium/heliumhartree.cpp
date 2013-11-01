#include "heliumhartree.h"

Helium::Helium()
{
    alpha = zeros(4);
    alpha(0) = 0.298073;
    alpha(1) = 1.242567;
    alpha(2) = 5.782948;
    alpha(3) = 38.474970;
}

double Helium::electronInteractionIntegral(int p, int r, int q, int s) {
    double denominator = (alpha(p) + alpha(q))*(alpha(r) + alpha(s))*sqrt(alpha(p) + alpha(q) + alpha(r) + alpha(s));
    return 2 * powPi5over2 / denominator;
}

/*!
 * \brief HartreeSolver::kineticIntegral is the solution of the the integral
 * $$\langle \chi_p | -\frac{1}{2} \nabla^2 | \chi_q \rangle$$
 * \param p
 * \param q
 * \return
 */
double Helium::kineticIntegral(int p, int q) {
    double alpha_p = alpha(p);
    double alpha_q = alpha(q);
    return 3*pow(M_PI, 3.0/2.0)*alpha_q/(pow(alpha_p, 3.0/2.0)*pow(1 + alpha_q/alpha_p, 5.0/2.0));
}

/*!
 * \brief HartreeSolver::kineticIntegral is the solution of the the integral
 * $$\langle \chi_p | -\frac{1}{2} \nabla^2 | \chi_q \rangle$$
 * \param p
 * \param q
 * \return
 */
double Helium::nuclearAttractionIntegral(int p, int q) {
    double alpha_p = alpha(p);
    double alpha_q = alpha(q);
    return -4*M_PI/(alpha_p*(1 + alpha_q/alpha_p));
}

double Helium::overlapIntegral(int p, int q) {
    double alpha_p = alpha(p);
    double alpha_q = alpha(q);
    return pow(M_PI, 3.0/2.0)/(pow(alpha_p, 3.0/2.0)*pow(1 + alpha_q/alpha_p, 3.0/2.0));
}

double Helium::additionalEnergyTerms()
{
    return 0;
}

uint Helium::nOrbitals()
{
    return 4;
}
