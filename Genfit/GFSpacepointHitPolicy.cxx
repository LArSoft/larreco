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
#include "Genfit/GFSpacepointHitPolicy.h"

#include "assert.h"

#include "TMath.h"

#include "Genfit/GFAbsRecoHit.h"

const std::string genf::GFSpacepointHitPolicy::fPolicyName = "GFSpacepointHitPolicy";


TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TMatrixT<Double_t> returnMat(2,1);

  TMatrixT<Double_t> _D(3,1);
  TVector3 _U;
  TVector3 _V;

  _D[0][0] = (plane.getO())[0];
  _D[1][0] = (plane.getO())[1];
  _D[2][0] = (plane.getO())[2];

  _D *= -1.; 
  _D += hit->getRawHitCoord();
  //now the vector _D points from the origin of the plane to the hit point


  _U = plane.getU();
  _V = plane.getV();


  returnMat[0][0] = _D[0][0] * _U[0] + _D[1][0] * _U[1] + _D[2][0] * _U[2];
  returnMat[1][0] = _D[0][0] * _V[0] + _D[1][0] * _V[1] + _D[2][0] * _V[2];
  //std::cout << "hitCoord="<<std::endl;
  //returnMat.Print();
  return returnMat;
}


TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane,const GFDetPlane& planePrev)
{
  TMatrixT<Double_t> returnMat(5,1); // Just return last 2 elements. Will calculate the rest in GFKalman.cxx.

  TMatrixT<Double_t> _D(3,1);
  TVector3 _U;
  TVector3 _V;

  _D[0][0] = (plane.getO())[0];
  _D[1][0] = (plane.getO())[1];
  _D[2][0] = (plane.getO())[2];

  _D *= -1.; 
  _D += hit->getRawHitCoord();
  //now the vector _D points from the origin of the plane to the hit point


  _U = plane.getU();
  _V = plane.getV();


  returnMat[0][0] = _D[0][0] * _U[0] + _D[1][0] * _U[1] + _D[2][0] * _U[2];
  returnMat[1][0] = _D[0][0] * _V[0] + _D[1][0] * _V[1] + _D[2][0] * _V[2];
  //std::cout << "hitCoord="<<std::endl;
  //returnMat.Print();
  return returnMat;
}


TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TVector3 _U;
  TVector3 _V;

  _U = plane.getU();
  _V = plane.getV();

  TMatrixT<Double_t> rawCov = hit->getRawHitCov();
  //rawCov[0][0] = 12.e12; // Force this. EC, 3-Apr-2012.

  TMatrixT<Double_t> jac(3,2);
  
  // jac = dF_i/dx_j = s_unitvec * t_untivec, with s=u,v and t=x,y,z
  jac[0][0] = _U[0];
  jac[1][0] = _U[1];
  jac[2][0] = _U[2];
  jac[0][1] = _V[0];
  jac[1][1] = _V[1];
  jac[2][1] = _V[2];

  TMatrixT<Double_t> jac_t = jac.T();
  TMatrixT<Double_t> jac_orig = jac;

  TMatrixT<Double_t> result=jac_t * (rawCov * jac_orig);
  //std::cout << "hitCov="<<std::endl;
  //result.Print();
  return  result;
}

TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane, const GFDetPlane& planePrev, const TMatrixT<Double_t>& state, const Double_t& mass)
{
  TVector3 _U;
  TVector3 _V;

  _U = plane.getU();
  _V = plane.getV();
  
//  Double_t eps(1.0e-6);
  TMatrixT<Double_t> rawCov = hit->getRawHitCov();

  //rawCov[0][0] = 12.e12; // Force this. EC, 3-Apr-2012.
  std::vector <double> tmpRawCov;
  tmpRawCov.push_back(rawCov[0][0]);
  tmpRawCov.push_back(rawCov[1][1]);
  tmpRawCov.push_back(rawCov[2][2]);
  tmpRawCov.push_back(rawCov[2][1]); // y-z correlated cov element.

  static std::vector <double> oldRawCov(tmpRawCov);
  static std::vector <double> oldOldRawCov(tmpRawCov);
  static GFDetPlane planePrevPrev(planePrev);
  rawCov.ResizeTo(7,7); // this becomes used now for something else entirely:
  // the 7x7 raw errors on x,y,z,px,py,pz,th, which sandwiched between 5x7
  // Jacobian converts it to the cov matrix for the 5x5 state space.
  rawCov[0][0] = tmpRawCov[0]; // x
  rawCov[1][1] = tmpRawCov[1]; // y
  rawCov[2][2] = tmpRawCov[2]; // z
  rawCov[1][2] = tmpRawCov[3]; // yz
  rawCov[2][1] = tmpRawCov[3]; // yz
  
  // This Jacobian concerns the transfrom from planar to detector coords.
  //   TMatrixT<Double_t> jac(3,2);
  TMatrixT<Double_t> jac(7,5); // X,Y,Z,UX,UY,UZ,Theta in detector coords
  
  // jac = dF_i/dx_j = s_unitvec * t_unitvec, with s=|q|/p,du/dw, dv/dw, u,v and t=th, x,y,z

  // More correctly :
  Double_t dist(0.3); // place holder!!
  Double_t C;
  //Double_t p (2.5); // place holder!!!
  //p = abs(1/state[0][0]);

  Double_t mom = fabs(1.0/state[0][0]);
  Double_t beta = mom/sqrt(mass*mass+mom*mom);
  dist = (plane.getO()-planePrev.getO()).Mag();
  C = 0.0136/beta*sqrt(dist/14.0)*(1+0.038*log(dist/14.0)); 
  TVector3 _D (plane.getNormal());

  // px^ = (x1-x2)/d; 
  // (sigpx)^2 = (sig_x1^2+sig_x2^2)*[(y1-y2)^2 + (z1-z2)^2]/d^4 +
  //                (x1-x2)^2*(y1-y2)^2*sigma_y1^2/d^6 + ...
  //                (x1-x2)^2*(z1-z2)^2*sigma_z2^2/d^6 + 

  rawCov[3][3] = (oldRawCov[0]+tmpRawCov[0])/pow(dist,4.) *
    ( 
     pow((plane.getO()-planePrev.getO()).Y(),2.) +
     pow((plane.getO()-planePrev.getO()).Z(),2.)
    )
    +
     pow((plane.getO()-planePrev.getO()).X(),2.) *
    (
     pow((plane.getO()-planePrev.getO()).Y(),2.) * (oldRawCov[1]+tmpRawCov[1]) + 
     pow((plane.getO()-planePrev.getO()).Z(),2.) * (oldRawCov[2]+tmpRawCov[2]) 
     ) / pow(dist,6.)
    // y-z cross-term
    +
    3.0* (plane.getO()-planePrev.getO()).X() * (plane.getO()-planePrev.getO()).Y() * (plane.getO()-planePrev.getO()).Z() / pow(dist,5.)  * (tmpRawCov[3] + oldRawCov[3])
    ;

  rawCov[4][4] = (oldRawCov[1]+tmpRawCov[1])/pow(dist,4.) *
    ( 
     pow((plane.getO()-planePrev.getO()).X(),2.) +
     pow((plane.getO()-planePrev.getO()).Z(),2.)
    )
    +
     pow((plane.getO()-planePrev.getO()).Y(),2.) *
    (
     pow((plane.getO()-planePrev.getO()).X(),2.) * (oldRawCov[0]+tmpRawCov[0]) + 
     pow((plane.getO()-planePrev.getO()).Z(),2.) * (oldRawCov[2]+tmpRawCov[2]) 
     ) / pow(dist,6.)
    // y-z cross-term
    +
    3.0* ( (plane.getO()-planePrev.getO()).X() * (plane.getO()-planePrev.getO()).Y() * (plane.getO()-planePrev.getO()).Z() / pow(dist,5.) + (plane.getO()-planePrev.getO()).Z() / pow(dist,3.) ) * ( tmpRawCov[3] + oldRawCov[3] )
    ;

  rawCov[5][5] = (oldRawCov[2]+tmpRawCov[2])/pow(dist,4.) *
    ( 
     pow((plane.getO()-planePrev.getO()).X(),2.) +
     pow((plane.getO()-planePrev.getO()).Y(),2.)
    )
    +
     pow((plane.getO()-planePrev.getO()).Z(),2.) *
    (
     pow((plane.getO()-planePrev.getO()).Y(),2.) * (oldRawCov[0]+tmpRawCov[0])+
     pow((plane.getO()-planePrev.getO()).X(),2.) * (oldRawCov[1]+tmpRawCov[1])
     ) / pow(dist,6.)
    +
    // y-z cross-term
    3.0* ( (plane.getO()-planePrev.getO()).X() * (plane.getO()-planePrev.getO()).Y() * (plane.getO()-planePrev.getO()).Z() / pow(dist,5.)  + (plane.getO()-planePrev.getO()).Y() / pow(dist,3.) ) * ( tmpRawCov[3] + oldRawCov[3] )
    ;
    
  
  Double_t num = (plane.getO()-planePrev.getO())*(planePrev.getO()-planePrevPrev.getO());
  Double_t d1 = (plane.getO()-planePrev.getO()).Mag(); // same as dist above
  Double_t d2 = (planePrev.getO()-planePrevPrev.getO()).Mag();

  // This is the error on cosTh^2. There are 9 terms for diagonal errors, one for each
  // x,y,z coordinate at 3 points. Below the first 3 terms are at pt 1. 
  // Second 3 are at 3. Third 3 are at 2, which is slightly more complicated.
  // There are three more terms for non-diagonal erros.
  double dcosTh = 
    pow(((planePrev.getO()-planePrevPrev.getO()).X()*d1*d2 -
		      num * (plane.getO()-planePrev.getO()).X() *d2/d1) /
		     pow(d1,2.0)/pow(d2,2.0),2.0) * tmpRawCov[0] +
    pow(((planePrev.getO()-planePrevPrev.getO()).Y()*d1*d2 + 
	 num * (plane.getO()-planePrev.getO()).Y() *d2/d1) /
	pow(d1,2.0)/pow(d2,2.0),2.0) * tmpRawCov[1] +
    pow(((planePrev.getO()-planePrevPrev.getO()).Z()*d1*d2 + 
	 num * (plane.getO()-planePrev.getO()).Z() *d2/d1) /
	pow(d1,2.0)/pow(d2,2.0),2.0) * tmpRawCov[2] +
    
    pow(((plane.getO()-planePrev.getO()).X()*d1*d2 - 
		      num * (planePrev.getO()-planePrevPrev.getO()).X() *d1/d2) /
		     pow(d1,2.0)/pow(d2,2.0),2.0) * oldOldRawCov[0] +
    pow(((plane.getO()-planePrev.getO()).Y()*d1*d2 - 
	 num * (planePrev.getO()-planePrevPrev.getO()).Y() *d1/d2) /
	pow(d1,2.0)/pow(d2,2.0),2.0) * oldOldRawCov[1] +
    pow(((plane.getO()-planePrev.getO()).Z()*d1*d2 - 
	 num * (planePrev.getO()-planePrevPrev.getO()).Z() *d1/d2) /
	pow(d1,2.0)/pow(d2,2.0),2.0) * oldOldRawCov[2] +
    
    pow(((plane.getO()-planePrevPrev.getO()).X()*d1*d2 - 
		      num * (plane.getO()-planePrev.getO()).X() *d2/d1 +
		      num * (planePrev.getO()-planePrevPrev.getO()).X() *d1/d2
	 ) /
		     pow(d1,2.0)/pow(d2,2.0),2.0) * oldRawCov[0] +
    pow(((plane.getO()-planePrevPrev.getO()).Y()*d1*d2 - 
		      num * (plane.getO()-planePrev.getO()).Y() *d2/d1 +
		      num * (planePrev.getO()-planePrevPrev.getO()).Y() *d1/d2
	 ) /
		     pow(d1,2.0)/pow(d2,2.0),2.0) * oldRawCov[1] +
    pow(((plane.getO()-planePrevPrev.getO()).Z()*d1*d2 - 
		      num * (plane.getO()-planePrev.getO()).Z() *d2/d1 +
		      num * (planePrev.getO()-planePrevPrev.getO()).Z() *d1/d2
	 ) /
		     pow(d1,2.0)/pow(d2,2.0),2.0) * oldRawCov[2]
    // And now the off-diagonal terms. This is a mess unto itself.
    // First, d^2(costh)/d(z1)d(y1)
    +
    (
     (plane.getO()-planePrev.getO()).Z() * (planePrev.getO()-planePrevPrev.getO()).Y() / pow(d1,3.)/d2 
     -
     (plane.getO()-planePrev.getO()).Y() * (planePrev.getO()-planePrevPrev.getO()).Z() / pow(d1,3.)/d2 
     +
     3.*(plane.getO()-planePrev.getO()).Y() * (plane.getO()-planePrev.getO()).Z() * num / pow(d1,5.)/d2 
     ) * tmpRawCov[3]
    // Next, d^2(costh)/d(z3)d(y3). Above w z1<->z3, y1<->y3
    +
    (
    (planePrevPrev.getO()-planePrev.getO()).Z() * (planePrev.getO()-plane.getO()).Y() / pow(d1,3.)/d2 
    -
    (planePrevPrev.getO()-planePrev.getO()).Y() * (planePrev.getO()-plane.getO()).Z() / pow(d1,3.)/d2 
    +
    3.*(planePrevPrev.getO()-planePrev.getO()).Y() * (planePrevPrev.getO()-planePrev.getO()).Z() * num / pow(d1,5.)/d2 
    ) * oldOldRawCov[3]
    // Last, d^2(costh)/d(z2)d(y2). This is an even bigger mess.
    +
    (
     ( (plane.getO()-planePrev.getO()).Y() - (planePrev.getO()-planePrevPrev.getO()).Y() ) *
     (-(plane.getO()-planePrev.getO()).Z()/pow(d1,3.)/d2 + 
       (planePrev.getO()-planePrevPrev.getO()).Z()/pow(d2,3.)/d1 
      )
     + (plane.getO()-planePrev.getO()).Y() * 
     (
      (-(planePrev.getO()-planePrevPrev.getO()).Z() + 
       (plane.getO()-planePrev.getO()).Z() 
       ) / pow(d1,3.)/d2
      - num*
      ( -3.* (plane.getO()-planePrev.getO()).Z()*d1*d2 +
	(planePrev.getO()-planePrevPrev.getO()).Z()*pow(d1,3.)/d2 
       )
      ) / pow(d1,6.) / pow(d2,2.)
     -
     (planePrev.getO()-planePrevPrev.getO()).Y()*
     (
      (-(planePrev.getO()-planePrevPrev.getO()).Z() + 
       (plane.getO()-planePrev.getO()).Z() 
       ) / pow(d2,3.)/d1
      - num*
      ( 3.* (planePrev.getO()-planePrevPrev.getO()).Z()*d1*d2 -
	(plane.getO()-planePrev.getO()).Z()*pow(d2,3.)/d1 
       )
      ) / pow(d1,2.) / pow(d2,6.)
     ) * oldRawCov[3]
    ;

  // That was delta(cos(theta))^2. I want delta(theta)^2. dTh^2 = d(cosTh)^2/sinTh^2
  Double_t theta(TMath::ACos((plane.getO()-planePrev.getO()).Unit() * (planePrev.getO()-planePrevPrev.getO()).Unit()));
  rawCov[6][6] = TMath::Min(pow(dcosTh,2.)/pow(TMath::Sin(theta),2.),pow(TMath::Pi()/2.0,2.));

  // This means I'm too close to the endpoints to have histories that allow
  // proper calculation of above rawCovs. Use below defaults instead.
  if (d1 == 0 || d2 == 0)
    {
      rawCov[3][3] = pow(0.2,2.0); // Unit Px
      rawCov[4][4] = pow(0.2,2.0); // Unit Py
      rawCov[5][5] = pow(0.2,2.0); // Unit Pz
      rawCov[6][6] = pow(0.1,2.0); // theta. 0.3/3mm, say.
      dist = 0.3;
      C = 0.0136/beta*sqrt(dist/14.0)*(1+0.038*log(dist/14.0)); 
    }

  // This forces a huge/tiny theta error, which effectively freezes/makes-us-sensitive theta, as is.
  //rawCov[6][6] = /*9.99e9*/ 0.1;
  
  TVector3 u=plane.getU();
  TVector3 v=plane.getV();
  TVector3 w=u.Cross(v);

  // Below line has been done, with code change in GFKalman that has updated
  // the plane orientation by now.
  //  TVector3 pTilde = 1.0 * (w + state[1][0] * u + state[2][0] * v);
  TVector3 pTilde = w;
  double pTildeMag = pTilde.Mag();

  jac.Zero();
  jac[6][0] = 1./C; // Should be 1/C?; Had 1 until ... 12-Feb-2013

  jac[0][3] = _U[0];
  jac[1][3] = _U[1];
  jac[2][3] = _U[2];

  jac[0][4] = _V[0];
  jac[1][4] = _V[1];
  jac[2][4] = _V[2];

  // cnp'd from RKTrackRep.cxx, line ~496
  // da{x,y,z}/du'
  jac[3][1] = 1.0/pTildeMag*(u.X()-pTilde.X()/(pTildeMag*pTildeMag)*u*pTilde);
  jac[4][1] = 1.0/pTildeMag*(u.Y()-pTilde.Y()/(pTildeMag*pTildeMag)*u*pTilde);
  jac[5][1] = 1.0/pTildeMag*(u.Z()-pTilde.Z()/(pTildeMag*pTildeMag)*u*pTilde);
  // da{x,y,z}/dv'
  jac[3][2] = 1.0/pTildeMag*(v.X()-pTilde.X()/(pTildeMag*pTildeMag)*v*pTilde);
  jac[4][2] = 1.0/pTildeMag*(v.Y()-pTilde.Y()/(pTildeMag*pTildeMag)*v*pTilde);
  jac[5][2] = 1.0/pTildeMag*(v.Z()-pTilde.Z()/(pTildeMag*pTildeMag)*v*pTilde);

  TMatrixT<Double_t> jac_orig = jac;
  TMatrixT<Double_t> jac_t = jac.T();

  TMatrixT<Double_t> result=jac_t * (rawCov * jac_orig);
  //std::cout << "hitCov="<<std::endl;
  //result.Print();

  
  oldOldRawCov = oldRawCov;
  oldRawCov = tmpRawCov;
  planePrevPrev = planePrev;
  
  return  result;
}

const genf::GFDetPlane&
genf::GFSpacepointHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  TMatrixT<Double_t> rawcoord = hit->getRawHitCoord();
  TVector3 point(rawcoord[0][0],rawcoord[1][0],rawcoord[2][0]);

  TVector3 poca,dirInPoca;
  rep->extrapolateToPoint(point,poca,dirInPoca);

  fPlane.setO(point);
  fPlane.setNormal(dirInPoca);

  return fPlane;
}

//ClassImp(GFSpacepointHitPolicy)
