#ifndef GAUSSIANCONTRACTEDORBITAL_H
#define GAUSSIANCONTRACTEDORBITAL_H

#include "math/vector3.h"
#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class GaussianContractedOrbital
{
public:
    GaussianContractedOrbital();
    GaussianContractedOrbital(const Vector3 &corePosition);

    void addPrimitiveBasisFunction(const GaussianPrimitiveOrbital& primitive);
    const GaussianPrimitiveOrbital& primitiveBasisFunction(int index);

    Vector3 corePosition() const;
    void setCorePosition(const Vector3 &corePosition);

    const vector<GaussianPrimitiveOrbital>& primitiveBasisFunctions() const;
    void setPrimitiveBasisFunctions(const vector<GaussianPrimitiveOrbital> &primitiveBasisFunctions);

    double evaluated(double x, double y, double z) const;

private:
    Vector3 m_corePosition;
    vector<GaussianPrimitiveOrbital> m_primitiveBasisFunctions;
};

#endif // GAUSSIANCONTRACTEDORBITAL_H
