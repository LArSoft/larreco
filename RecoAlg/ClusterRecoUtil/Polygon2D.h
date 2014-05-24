#ifndef POLYGON2D_H
#define POLYGON2D_H

//
// Author: David Caratelli @ Columbia U.
//

#include <vector>
#include <utility>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>

//a polygon is a vector of std::pairs with first = x coordinate
//and second = y coordinate of that vertex
//access vertices with Point function. Points are:
//0   = first ordered vertex
//n-1 = last ordered vertex (n=size of polygon)
//n   = first vertex again
//>n  = invalid...return error message

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//                BEGIN POLYGON CLASS               //
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
class Polygon2D{
  
private:
  std::vector< std::pair<float,float> > vertices;
  
 public:

  Polygon2D() { }
  Polygon2D(const std::vector< std::pair<float,float> > &points) { vertices = points; }
  unsigned int Size() const { return vertices.size(); } 
  const std::pair<float,float>& Point(unsigned int p) const; 
  std::pair<float,float> Project(const std::pair<float,float>&,float) const;
  float Area() const;
  bool Overlap(float slope, const Polygon2D &poly2, const std::pair<float,float> &origin) const;
  bool PolyOverlap(const Polygon2D &poly2) const;
  bool PolyOverlapSegments(const Polygon2D &poly2) const;
  bool PointInside(const std::pair<float,float> &point) const;
  bool Contained(const Polygon2D &poly2) const;
};


#endif
