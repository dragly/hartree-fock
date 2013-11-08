#ifndef GAUSSIANCONTRACTEDORBITAL_H
#define GAUSSIANCONTRACTEDORBITAL_H

#include <basisfunctions/gaussian/gaussianprimitiveorbital.h>

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class GaussianContractedOrbital
{
public:
//    GaussianContractedOrbital();
    GaussianContractedOrbital(const rowvec &corePosition);

    void addPrimitiveBasisFunction(const GaussianPrimitiveOrbital& primitive);
    const GaussianPrimitiveOrbital& primitiveBasisFunction(int index);

    rowvec corePosition() const;
    void setCorePosition(const rowvec &corePosition);

    const vector<GaussianPrimitiveOrbital>& primitiveBasisFunctions() const;
    void setPrimitiveBasisFunctions(const vector<GaussianPrimitiveOrbital> &primitiveBasisFunctions);

private:
    rowvec m_corePosition;
    vector<GaussianPrimitiveOrbital> m_primitiveBasisFunctions;
};

#endif // GAUSSIANCONTRACTEDORBITAL_H
