/**********************************************************************************//**
 * \file main_convert.cpp
 *
 * \author Kjetil A. Johannessen
 * \date November 2013
 *
 *************************************************************************************/

#include <iostream>
#include <sstream>
#include <set>
#include <fstream>
#include <string.h>

#include "SplineModel.h"
#include "TopologySet.h"
#include "primitives.h"

#include <GoTools/trivariate/SplineVolume.h>

/*
 *        Nodal enumeration of the resulting hex
 *
 *             5              6
 *              .+----------+  
 *            .' |        .'| 
 *       4  .'   |   7  .'  | 
 *         +-----+----+'    | 
 *         |     |    |     |          
 *         |   1 |    |     |           ^
 *         |    ,+----+-----+ 2         |   _, y
 *         |  .'      |   .'          z |  ,'|
 *         |.'        | .'              |.' 
 *         +----------+'                +--------> 
 *       0             3                     x
 */


typedef unsigned int uint;

using namespace std;
using namespace Go;

bool stl               = false;
bool onlyTets          = false;
bool oneIndexed        = false;
bool storeFile         = false;
bool uniform           = false;
int evaluationPoints   = 3;
const char *outputFile = new char[1024];
string fileUsage       = "\
File usage: convert [-out=<prefix>] [-tetr] [-1] [-ev=<n>] [-unif] <inputFile> \n\
  \n\
  Arguments\n\
    <inputFile> : one or more .g2-files describing the spline volumes \n\
  FLAGS\n\
    -out        : stores the resulting FEM model in the files starting with <prefix>\n\
    -tetr       : stores element list tetrahedrals and not hexahedrals \n\
    -stl        : create an stl-file rather than elements\n\
    -1          : stores element list as 1-indexed\n\
    -ev         : evaluates geometry with <n> points per knot span [default=3]\n\
    -unif       : evaluates geometry using <n> points per patch, ignoring knot spans\n\
    -help       : display this help screen";
SplineModel model;

void processParameters(int argc, char** argv) {
	bool fileRead = false;
	outputFile    = "out"; // default file prefix
	for(int argi=1; argi<argc; argi++) {
		const char *arg = argv[argi];
		if(strcmp(arg, "-help") == 0) {
			cout << fileUsage << endl;
			exit(0);
		} else if(strncmp(arg, "-out=", 5) == 0) {
			outputFile = arg+5;
			storeFile  = true;
		} else if(strncmp(arg, "-ev=", 4) == 0) {
			evaluationPoints = atoi(arg+4);
		} else if(strncmp(arg, "-unif", 5) == 0) {
			uniform = true;
		} else if(strncmp(arg, "-1", 2) == 0) {
			oneIndexed = true;
		} else if(strncmp(arg, "-stl", 4) == 0) {
			stl = true;
		} else if(strncmp(arg, "-tetr", 5) == 0) {
			onlyTets = true;
		} else {
			ifstream inFile;
			inFile.open(arg);
			if(!inFile.good()) {
				cerr << "Error reading input file \"" << arg << "\"" << endl;
				exit(1);
			}
			model.readSplines(inFile, false);
			inFile.close();
			fileRead = true;
		}
	}
	if(!fileRead) {
		cout << fileUsage << endl;
		exit(1);
	}
	cout << boolalpha;
	cout << "Converting spline model to a classical finite element model" << endl;
	cout << "Parameters: " << endl;
	cout << "  Model type                    : " << ((model.isVolumetricModel())?"Volumetric":"Surface") << endl;
	cout << "  Number of spline patches      : " << model.getNumbPatches() << endl;
	cout << "  Number of evaluation points   : " << evaluationPoints << endl;
	cout << "  Evaluation type               : " << ((uniform)?"Per patch":"Per knot-span") << endl;
	cout << "  Element type                  : " << ((stl) ? ("stl"):((onlyTets)?("Tetrahedral"):("Hexahedral"))) << endl;
	if(!stl)
	cout << "  Store results using 1-indexing: " << oneIndexed << endl;
	cout << "  Output file prefix            : \"" << outputFile << "\"" << endl;
	cout << endl;
}

/**********************************************************************************//**
 * \brief clears up any degenerate tetrahedral elements 
 * \param elms list of elements
 * \details Clears out any elements whos nodal points appear in multiple corners.
 *          For tetrahedral elements it makes sense to remove these elements since they will
 *          always have zero volume, but for more general elements (hex or prisms), you might
 *          want to change this to alter the element type from hex to prism or from prism to
 *          tetrahedral.
 *************************************************************************************/
void clearDegenerices (std::vector<std::vector<int> >& elms) {
	vector<vector<int> > result;
	for(vector<int> el : elms) {

		// add up all unique node numbers
		std::set<int> nodes; 
		for(int index : el) 
			nodes.insert(index); 

		// if any duplicate nodal points, don't include in results
		if(nodes.size() == el.size())
			result.push_back(el);
	}

	// overwrite the input element list with the non-degenerate one
	elms.clear();
	elms.resize(result.size());
	copy(result.begin(), result.end(), elms.begin());
}

/**********************************************************************************//**
 * \brief converts quadrilateral elements to triangle ones
 * \param quadrilaterals[in] list of elements (4 node elements only)
 * \param triangles[out]     list of elements (3 node elements)
 *************************************************************************************/
void quad2tri (const std::vector<std::vector<int> >& quadrilaterals, std::vector<std::vector<int> >& triangles) {
	int q2t[][3] = {{0,1,2},
	                {2,3,0}};  // quad-to-tri conversion matrix
	
	triangles.clear();
	for(vector<int> elm : quadrilaterals) {
		for(int i=0; i<2; i++) {
			vector<int> tri(3);
			for(int j=0; j<3; j++)
				tri[j] = elm[q2t[i][j]];
			triangles.push_back(tri);
		}
	}
}

/**********************************************************************************//**
 * \brief converts hexahedral elements to tetrahedral ones
 * \param hexahedrals[in]   list of elements (8 node elements only)
 * \param tetrahedrals[out] list of elements (4 node elements)
 *************************************************************************************/
void hex2tetr (const std::vector<std::vector<int> >& hexahedrals, std::vector<std::vector<int> >& tetrahedrals) {
	int h2t[][4] = {{0,1,5,6},
	                {0,4,5,6},
	                {0,1,2,6},
	                {0,2,3,6},
	                {0,3,6,7},
	                {0,4,6,7}}; // hex-to-tet conversion matrix
	
	tetrahedrals.clear();
	for(vector<int> elm : hexahedrals) {
		for(int i=0; i<6; i++) {
			vector<int> tet(4);
			for(int j=0; j<4; j++)
				tet[j] = elm[h2t[i][j]];
			tetrahedrals.push_back(tet);
		}
	}
}

/**********************************************************************************//**
 * \brief converts hexahedral elements to tetrahedral ones
 * \param outFile output stream 
 * \param tri     list of triangles
 * \param pts     list of points
 *************************************************************************************/
void writeSTL (std::ostream& outFile, const std::vector<std::vector<int> >& tri, const std::vector<Go::Point>& pts) {
	outFile << setprecision(16);
	outFile << "solid SplineModel" << endl;
	for(vector<int> el : tri) {
		outFile << "facet normal 0 0 0" << endl;
		outFile << "  outer loop" << endl;
		for(int i : el)
			outFile << "    vertex " << pts[i] << endl;
		outFile << "  endloop" << endl;
		outFile << "endfacet"  << endl;
	}
	outFile << "endsolid SplineModel" << endl;
}

/**********************************************************************************//**
 * \brief Main program for converting spline geometry to classical FEM geometry
 *
 * This samples the multipatch spline model (with possible degeneracies) and creates a
 * traditional finite element representation of this. The output format can be given
 * as hexahedral 8-node elements or tetrahedral 4-node elements in case of volume models.
 * Support for surface models is not yet implemented.
 *
 *************************************************************************************/

int main(int argc, char **argv) {
	processParameters(argc, argv);
	model.enforceRightHandSystem();
	model.buildTopology();
	model.generateGlobalNumbers();

	vector<Go::Point>    pts;
	vector<vector<int> > hex;
	vector<vector<int> > tets;
	vector<vector<int> > bndryHex;
	vector<vector<int> > bndryTet;

	cout << "Evaluating mesh..." << endl;
	model.getTesselation(pts, hex, bndryHex, evaluationPoints, uniform);
	
	if(onlyTets || stl) {
		cout << "Converting to tetrahedral elements..." << endl;
		if(!stl)
			hex2tetr(hex,tets);
		quad2tri(bndryHex,bndryTet);
		cout << "Removing degenerate tetrahedrals..." << endl;
		if(!stl)
			clearDegenerices(tets);
		clearDegenerices(bndryTet);
	}

/*
	for(int patch=0; patch<model.getNumbPatches(); patch++) {
		int n1 = model.getNumbPts(patch, 0);
		int n2 = model.getNumbPts(patch, 1);
		int n3 = model.getNumbPts(patch, 2);
		for(int w=0; w<n3; w++) {
			for(int v=0; v<n2; v++) {
				for(int u=0; u<n1; u++)
					cout << model.getGlobalNumber(patch, u,v,w) << " " ;
				cout << "\n";
			}
			cout << "\n";
		}
		cout << "-------------\n\n";
	}
	for(vector<int> el : ((onlyTets)?bndryTet:bndryHex))  {
		for(int d : el)
			cout << d << " ";
		cout << endl;
	}
*/
	

	cout << "Writing results to file..." << endl;
	if(stl) {
		string stlFile    = string(outputFile) + string(".stl");
		ofstream outFile(stlFile);
		writeSTL(outFile, bndryTet, pts);
		outFile.close();
		cout << "Results written to \"" << stlFile << "\"" << endl;
	} else {
		string nodeFile    = string(outputFile) + string("_nodes.m");
		string elementFile = string(outputFile) + string("_element.m");
		string bndryFile   = string(outputFile) + string("_bndry.m");
		ofstream outFile1(nodeFile);
		ofstream outFile2(elementFile);
		ofstream outFile3(bndryFile);
		outFile1 << setprecision(16);
	
		for(Go::Point p : pts)
			outFile1 << p << endl;
		outFile1.close();
		cout << "Nodal points successfully written to \"" << nodeFile << "\"" << endl;
	
		for(vector<int> el : ((onlyTets)?tets:hex))  {
			for(int i : el)
				outFile2 << (i+oneIndexed) << " ";
			outFile2  << endl;
		}
		outFile2.close();
		cout << "Elements successfully written to \"" << elementFile << "\"" << endl;
	
		for(vector<int> el : ((onlyTets)?bndryTet:bndryHex))  {
			for(int i : el)
				outFile3 << (i+oneIndexed) << " ";
			outFile3  << endl;
		}
		outFile3.close();
		cout << "Boundary successfully written to \"" << bndryFile << "\"" << endl;
	}
	
	stringstream ss;
	ss.imbue(std::locale(""));

	ss << endl;
	ss << "Summary: " << endl;
	if(stl) {
		ss << "  Total numb. of triangles: " << bndryTet.size() << endl;
	} else {
		ss << "  Total nodes written     : " << pts.size() << endl;
		ss << "  Total elements written  : " << ((onlyTets)?tets.size():hex.size()) << endl;
		ss << "  Total boundary elements : " << ((onlyTets)?bndryTet.size():bndryHex.size()) << endl;
		ss << "  Element type            : " << ((onlyTets)?"4-node tetrahedral":"8-node hexahedral") << endl;
	}
	cout << ss.str() << endl;


	return 0;
}

