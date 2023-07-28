#include <fstream>
#include "Polyhedron.hpp"
#include "Triangle.hpp"
#include "Read_off.hpp"


void read_off(const std::string filename, Polyhedron& A)
{
	if (filename.empty()) {
		return;
	}
	std::ifstream fin;
	fin.open(filename);

	if (!fin)
	{
		printf("文件有误\n");
		return;
	}
	else
	{
		printf("文件打开成功\n");

		// 读取OFF字符串
		int nVertices = 0;
		int nFaces = 0;
		int nEdges = 0;
		std::vector<Eigen::Vector3f> vertices;
		std::string str;
		fin >> str;
		// 读取文件中顶点数、面片数、边数
		fin >> nVertices >> nFaces >> nEdges;

		// 根据顶点数，循环读取每个顶点坐标，将其保存到vertices
		float x, y, z;
		for (int i = 0; i < nVertices; i++) {
			fin >> x >> y >> z;
			vertices.push_back({ x, y, z });
		}

		// 根据面片数，循环读取每个面片信息，并用构建的vec3i结构体保存到faces
		unsigned int n, a, b, c;
		for (int i = 0; i < nFaces; i++) {
			fin >> n >> a >> b >> c;
			Triangle* face = new Triangle();
			face->setVertex(0, vertices[a]);
			face->setVertex(1, vertices[b]);
			face->setVertex(2, vertices[c]);
			A.TriangleList.push_back(face);
		}
	}
	fin.close();
}
