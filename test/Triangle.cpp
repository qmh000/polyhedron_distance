#include <Eigen/Eigen>
#include "Triangle.hpp"

Triangle::Triangle()
{
    v[0] << 0, 0, 0;
    v[1] << 0, 0, 0;
    v[2] << 0, 0, 0;
}

void Triangle::setVertex(int ind, Eigen::Vector3f ver) 
{ 
    v[ind] = ver; 
}



