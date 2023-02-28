/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
// ----------------------------------------------------------
// Please see GFWirepointHitPolicy.h  before using this class.
// ----------------------------------------------------------

#include "larreco/Genfit/GFWirepointHitPolicy.h"

#include <cassert>
#include <cmath>

#include "TMath.h"
#include "TVector3.h"

#include "larreco/Genfit/GFAbsRecoHit.h"
#include "larreco/Genfit/GFAbsTrackRep.h"
#include "larreco/Genfit/GFException.h"

const std::string genf::GFWirepointHitPolicy::fPolicyName = "GFWirepointHitPolicy";

genf::GFWirepointHitPolicy::GFWirepointHitPolicy() : fMaxdistance(1.E50)
{
  ;
}

TMatrixT<Double_t> genf::GFWirepointHitPolicy::hitCoord(GFAbsRecoHit* hit, const GFDetPlane& plane)
{
  TMatrixT<Double_t> returnMat(2, 1);

  checkPlane(hit, plane);

  // raw x1, y1, z1, x2, y2, z2, rdrift, zreco
  TMatrixT<Double_t> rC = hit->getRawHitCoord();

  returnMat[0][0] = rC[6][0];
  returnMat[1][0] = rC[7][0];
  return returnMat;
}

TMatrixT<Double_t> genf::GFWirepointHitPolicy::hitCov(GFAbsRecoHit* hit, const GFDetPlane& plane)
{
  checkPlane(hit, plane);

  TMatrixT<Double_t> returnCov(2, 2);
  TMatrixT<Double_t> rawCov = hit->getRawHitCov();

  returnCov[0][0] = rawCov[6][6];
  returnCov[1][0] = rawCov[7][6];
  returnCov[0][1] = rawCov[6][7];
  returnCov[1][1] = rawCov[7][7];

  return returnCov;
}

void genf::GFWirepointHitPolicy::checkPlane(GFAbsRecoHit* hit, const GFDetPlane& plane)
{
  // raw x1, y1, z1, x2, y2, z2, rdrift, zreco
  TMatrixT<Double_t> rC = hit->getRawHitCoord();

  assert(rC.GetNrows() == 8);

  TVector3 wire1(rC[0][0], rC[1][0], rC[2][0]);
  TVector3 wire2(rC[3][0], rC[4][0], rC[5][0]);
  TVector3 wiredirection = wire1 - wire2;

  TVector3 vaxis = plane.getV();
  wiredirection.SetMag(1.);
  vaxis.SetMag(1.);

  if (fabs(TMath::Abs(wiredirection.Dot(vaxis)) - 1) > 1e-3) {

    std::cout << "GFWirepointHitPolicy: plane not valid!!" << std::endl;
  }
}

const genf::GFDetPlane& genf::GFWirepointHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{

  TMatrixT<Double_t> x = hit->getRawHitCoord();
  assert(x.GetNrows() == 8);
  TVector3 wire1(x[0][0], x[1][0], x[2][0]);
  TVector3 wire2(x[3][0], x[4][0], x[5][0]);

  //  distance of one (the first) of the wire extremities from the plane
  Double_t d_from_refplane = fDetPlane.dist(wire1).Mag();
  if (d_from_refplane < 1e-5) return fDetPlane;

  // point of closest approach
  TVector3 poca, poca_onwire, dirInPoca;

  rep->extrapolateToLine(wire1, wire2, poca, dirInPoca, poca_onwire);

  Double_t distance;
  distance = TMath::Sqrt(fabs(
    ((wire1 - poca).Mag2() * (wire2 - wire1).Mag2() - pow((wire1 - poca).Dot(wire2 - wire1), 2)) /
    (wire2 - wire1).Mag2()));

  // check poca inside tube
  if (distance > fMaxdistance) {
    throw GFException("distance poca-wire > maxdistance", __LINE__, __FILE__) /* .setFatal() */;
  }

  // find plane
  // unitary vector along distance
  // poca (on track), poca_onwire (on wire)
  TVector3 fromwiretoextr = poca - poca_onwire;
  fromwiretoextr.SetMag(1.);
  // unitary vector along the wire
  TVector3 wiredirection = wire2 - wire1;
  wiredirection.SetMag(1.);

  // check orthogonality
  if (fabs(fromwiretoextr * wiredirection) > 1e-3) {
    throw GFException("fromwiretoextr*wiredirection > 1e-3", __LINE__, __FILE__) /* .setFatal() */;
  }

  TVector3 U;
  U = fromwiretoextr;
  TVector3 V;
  V = wiredirection;
  U.SetMag(1.);
  V.SetMag(1.);

  TVector3 O = (wire1 + wire2) * 0.5;

  fDetPlane = GFDetPlane(O, U, V);

  return fDetPlane;
}

//ClassImp(GFWirepointHitPolicy)
