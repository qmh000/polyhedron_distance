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
		printf("�ļ�����\n");
		return;
	}
	else
	{
		printf("�ļ��򿪳ɹ�\n");

		// ��ȡOFF�ַ���
		int nVertices = 0;
		int nFaces = 0;
		int nEdges = 0;
		std::vector<Eigen::Vector3f> vertices;
		std::string str;
		fin >> str;
		// ��ȡ�ļ��ж���������Ƭ��������
		fin >> nVertices >> nFaces >> nEdges;

		// ���ݶ�������ѭ����ȡÿ���������꣬���䱣�浽vertices
		float x, y, z;
		for (int i = 0; i < nVertices; i++) {
			fin >> x >> y >> z;
			vertices.push_back({ x, y, z });
		}

		// ������Ƭ����ѭ����ȡÿ����Ƭ��Ϣ�����ù�����vec3i�ṹ�屣�浽faces
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
