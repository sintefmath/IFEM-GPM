#ifndef _PRIMITIVES_H
#define _PRIMITIVES_H

#include <GoTools/utils/Point.h>
#include <GoTools/geometry/Streamable.h>
#include <vector>
#include <set>
#include <iostream>


class Volume;
class Face;
class Line;
class Vertex;

class Volume {
	public:
		Volume(int id) ;
		std::vector<int> getVertexEnumeration(Vertex *v);
		void getEdgeEnumeration(Line *l, std::vector<int> &numb, std::vector<int> &parDir, std::vector<int> &parStep);
		std::vector<int> getSurfaceEnumeration(Face *f);

		Vertex *corner[8];
		Line   *line[12];
		Face   *face[6];
		int    id;
};

class Face : public Go::Streamable {
	public:
		Face(int n1, int n2);
		bool equals(Face *f, double tol);
		void write(std::ostream &os) const;
		void read(std::istream &is);
		std::vector<std::vector<Go::Point> > cp;

		Volume *v1;
		Volume *v2;
		int    face1;  // face index on volume 1
		int    face2;  // face index on volume 2
		bool   uv_flip;
		bool   u_reverse;
		bool   v_reverse;
};

class Line : public Go::Streamable {
	public:
		bool equals(Line *l, double tol);
		void write(std::ostream &os) const;
		void read(std::istream &is);
		static int getLineEnumeration(int vert1, int vert2);

		std::vector<Go::Point> cp;
		std::set<Volume*> volume;
		Vertex *v1;
		Vertex *v2;
};

class Vertex {
	public:
		Go::Point cp;
		std::set<Volume*> volume;
};

#endif
