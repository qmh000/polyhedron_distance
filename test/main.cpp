#include <iostream>
#include "Read_off.hpp"
#include "Triangle.hpp"
#include "Polyhedron.hpp"


int main()
{
	Polyhedron obj1, obj2;
	read_off("nuclei.off", obj1);
	read_off("vessel.off", obj2);
	std::cout << obj1.PolyDist(obj2) << std::endl;
}
