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

/**********************************************************************************//**
 * \brief Class describing topological volume primitives
 *
 * A volume is defined by 8 corner vertices, 12 edge lines and 6 faces (collectively referred 
 * to as primitives). The primitives are stored as pointers and carefully not duplicated such that
 * equality testing can be done by address testing. Moreover, the class gives support for degenerate
 * volumes, making it possible that one unique vertex is stored in several of the 8 corner positions,
 * or likewise for the lines and faces.
 *
 *************************************************************************************/
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
		const char* material_code; //!< Material properties
};

/**********************************************************************************//**
 * \brief Class describing topological face primitives
 *
 * A face is is defined by its internal controlpoints, which are also used for checking
 * uniqueness. Once two (or more) faces are found equal, the form a twin-pair with access 
 * to both parent volumes. There are eight different ways that these faces can be oriented 
 * with respect to each other, and these are stored in the boolean variables uv_flip, 
 * u_reverse and v_reverse. Note that in the case of degenerate surfaces, it is possible
 * to neighbour an arbitrary number of volumes.
 *
 *************************************************************************************/
class Face : public Go::Streamable {
	public:
		Face(int id);
		Face(int n1, int n2);
		bool equals(Face *f, double tol);
		bool isDegen();
		std::vector<int> getVertexEnumeration(Vertex *v);
		void getEdgeEnumeration(Line *l, std::vector<int> &numb, std::vector<int> &parDir, std::vector<int> &parStep);
		void write(std::ostream &os) const;
		void read(std::istream &is);

		std::vector<std::vector<Go::Point> > cp; //!< Defining control points
		Vertex *corner[4];           //!< The four corners (only used in surface models)
		Line   *line[4];             //!< The four edges (only used in surface models)
		std::vector<Volume*> volume;     //!< Neighbouring volumes
		std::vector<int >    face;       //!< face index on volume 1 (volumetric models only)
		std::vector<bool>    uv_flip;    //!< is the first/second parameter direction swapped wrt the first volume (volumetric models only)
		std::vector<bool>    u_reverse;  //!< is the first parameter direction reversed (volumetric models only)
		std::vector<bool>    v_reverse;  //!< is the second parameter direction reversed (volumetric models only)

		bool   degen1;            //!< is the first parameter direction degenerated
		bool   degen2;            //!< is the second parameter direction degenerated

		int    id;
		const char* bc_code;      //!< Boundary condition code
};

/**********************************************************************************//**
 * \brief Class describing topological line primitives
 *
 * A line is defined by its internal controlpoints and contains access to all neighbouring volumes
 * as well as its start- and end-vertex
 *
 *************************************************************************************/
class Line : public Go::Streamable {
	public:
		Line();
		~Line();
		bool equals(Line *l, double tol);
		void write(std::ostream &os) const;
		void read(std::istream &is);
		static std::vector<int> getLineEnumeration(int face);
		static int              getLineEnumeration(int vert1, int vert2);

		std::vector<Go::Point> cp; //!< Defining control points
		std::set<Volume*> volume;  //!< Neighbouring volumes
		std::set<Face*> face;      //!< Neighbouring faces
		Vertex *v1;                //!< Start or stop vertex.
		Vertex *v2;                //!< Start or stop vertex.
		const char* bc_code;       //!< Boundary condition code
		bool degen;                //!< Is the line degenerated to a point
};

/**********************************************************************************//**
 * \brief Class describing topological vertex primitives
 *
 * A vertex is a single topological corner point, defined by it's corresponding controlpoint
 * and has access to all neighbouring volumes
 *
 *************************************************************************************/
class Vertex {
	public:
		Go::Point cp;              //!< Defining control point
		std::set<Volume*> volume;  //!< Neighbouring volumes
		std::set<Face*> face;      //!< Neighbouring face
		const char* bc_code;       //!< Boundary condition code

		Vertex();
		static void             getVertexEnumerationOnFace(int line, int &v1, int &v2);
		static void             getVertexEnumerationOnVolume(int line, int &v1, int &v2);
		static std::vector<int> getVertexEnumeration(int face);
};

#endif
