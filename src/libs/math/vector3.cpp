//#include <boost/mpi.hpp>

#include "vector3.h"

std::ostream& operator<< (std::ostream &out, const Vector3 &vector)
{
    out << vector.mem_local[0] << ", " << vector.mem_local[1] << ", " << vector.mem_local[2];
    return out;
}

Vector3 Vector3::ones()
{
    return Vector3(1,1,1);
}


Vector3 Vector3::createZeros()
{
    return Vector3(0,0,0);
}
