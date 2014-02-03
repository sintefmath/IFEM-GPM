#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <GoTools/geometry/ObjectHeader.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/trivariate/SplineVolume.h>

using namespace Go;

/*!
  \brief Program for raising order of NURBS finite element.

  This program increases the order of a spline finite element
  mesh in the form of a g2-file by one. 

  The input to the program is specified through the following 
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file: Input file in g2-format
  \arg -2D:           Input file with spline surfaces
  \arg -3D:           Input file with spline volumes
*/


int main(int argc, char** argv)
{
  bool surface_model = true;
  char* infile = 0;

  size_t i = 0;
  for (i = 1; i < argc; i++)
    if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;
  
  i=0;
  while (i < strlen(infile) && isspace(infile[i])) i++;
  std::ifstream isp(infile+i);

  // For spline surface models
  std::vector<SplineSurface*> surf;

  // For spline volume models
  std::vector<SplineVolume*> vol;

  ObjectHeader head;
  int n = 0;
  while (!isp.eof()) {
    head.read(isp);
    if (head.classType() == Class_SplineVolume) {
      SplineVolume* v(new SplineVolume());
      v->read(isp);
      vol.push_back(v);
      surface_model = false;
    }
    else if (head.classType() == Class_SplineSurface) {
      SplineSurface* s(new SplineSurface());
      s->read(isp);
      surf.push_back(s);
      surface_model = true;
    }
    else
      std::cerr << "Unknown spline model" << std::endl;
    
    // Ignore blanks
    ws(isp); 
  }

  if (surface_model) {
    std::vector<SplineSurface*>::iterator s_it;
    for (s_it = surf.begin();s_it != surf.end();s_it++) {
      (*s_it)->raiseOrder(1,1);
      (*s_it)->writeStandardHeader(std::cout);
      (*s_it)->write(std::cout);
    }
  
    for (i = 0;i < surf.size();i++) 
      delete surf[i];
  }
  else {
    std::vector<SplineVolume*>::iterator v_it;
    for (v_it = vol.begin();v_it != vol.end();v_it++) {
      (*v_it)->raiseOrder(1,1,1);
      (*v_it)->writeStandardHeader(std::cout);
      (*v_it)->write(std::cout);
    }

    for (i = 0;i < vol.size();i++) 
      delete vol[i];
  }

  return 0;
}  
