//
//  VertexFitMinuitStruct.h
//  work
//
//  Created by Bruce Baller on 6/18/15.
//  Copyright (c) 2015 Bruce Baller. All rights reserved.
//

#ifndef VertexFitMinuitStruct_h

#include "TVector3.h"

#include <array>
#include <vector>

struct VertexFitMinuitStruct {

  unsigned short TPC;
  unsigned short Cstat;
  unsigned short NPlanes;
  double WirePitch;
  std::array<double, 3> XFactor;               // The denominator factor in ConvertXToTicks
  std::array<double, 3> TickOff;               // The tick offset in ConvertXToTicks
  std::array<double, 3> OrthY;
  std::array<double, 3> OrthZ;
  std::array<double, 3> FirstWire;             // the FirstWireProj in WireCoordinate
  TVector3 VtxPos;           // Vertex position (detector units)

  std::vector<std::vector<double>> HitX; // hit X
  //   track      X
  std::vector< std::vector<double>> HitXErr; // hit X errors
  std::vector<std::vector<unsigned short>> Plane;
  std::vector<std::vector<unsigned short>> Wire;
  //   track
  std::vector<TVector3> Dir;
  std::vector<TVector3> DirErr;
  double DoF;
  float ChiDoF;                             // fit Chisq/DOF

};

#define VertexFitMinuitStruct_h


#endif
