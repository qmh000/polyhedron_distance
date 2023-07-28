#pragma once
#include <vector>
#include <Eigen/Eigen>
#include "Triangle.hpp"

class Polyhedron {
public:
	std::vector<Triangle*> TriangleList;
	float PolyDist(Polyhedron B);
};