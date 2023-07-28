#pragma once
#include <Eigen/Eigen>

class Triangle
{
public:
	Eigen::Vector3f v[3];
	Triangle();

	Eigen::Vector3f a() const { return v[0]; }
	Eigen::Vector3f b() const { return v[1]; }
	Eigen::Vector3f c() const { return v[2]; }

	void setVertex(int ind, Eigen::Vector3f ver);
};