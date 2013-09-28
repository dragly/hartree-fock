#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

#include <sys/types.h>

class BasisFunction
{
public:
    BasisFunction();

    virtual double electronInteractionIntegral(int p, int r, int q, int s) = 0;
    virtual double kineticIntegral(int p, int q) = 0;
    virtual double nuclearAttractionIntegral(int p, int q) = 0;
    virtual double overlapIntegral(int p, int q) = 0;

    virtual uint nOrbitals() = 0;
};

#endif // BASISFUNCTION_H
