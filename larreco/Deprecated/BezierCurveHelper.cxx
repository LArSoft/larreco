//////////////////////////////////////////////////////////////////////
///
/// \file   BezierCurveHelper.cxx
///
/// \brief  Helper object for interpolating between track segments
///
/// \author B J P Jones
///
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>

#include "lardataobj/RecoBase/Hit.h"
#include "larreco/Deprecated/BezierCurveHelper.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TVector.h"


namespace trkf{
//----------------------------------------------------------------------
// Constructor.
//
BezierCurveHelper::BezierCurveHelper()
{
  fCurveResolution=100;
}

BezierCurveHelper::BezierCurveHelper(int CurveRes)
{
  fCurveResolution=CurveRes;
}

//----------------------------------------------------------------------
// Destructor.
//
BezierCurveHelper::~BezierCurveHelper()
{
}


//----------------------------------------------------------------------
// Find the total length of one bezier segment connecting two seeds
// using a brute force integration
//

double BezierCurveHelper::GetSegmentLength(recob::Seed const& s1, recob::Seed const& s2)
{
  double Length=0;
  std::vector<TVector3> BezierPoints = GetBezierPoints(s1, s2, fCurveResolution);
  for(unsigned int p=0; p!=BezierPoints.size()-1; ++p)
    {
      Length+=(BezierPoints.at(p+1)-BezierPoints.at(p)).Mag();
    }
  return Length;
}


//----------------------------------------------------------------------
// Helper method for interpolation - find optimal scales for
// direction vectors
//

void BezierCurveHelper::GetDirectionScales(double * Pt1, double * Pt2, double * Dir1, double * Dir2, double * Scales)
{

  double PtSepVec[3];
  for(int i=0; i!=3; ++i) PtSepVec[i]=Pt2[i]-Pt1[i];

  double PtSep = pow(
    pow(PtSepVec[0],2) +
    pow(PtSepVec[1],2) +
    pow(PtSepVec[2],2),0.5);

  double Dir1Length = pow(Dir1[0]*Dir1[0]+
        Dir1[1]*Dir1[1]+
        Dir1[2]*Dir1[2],0.5);

  double Dir2Length = pow(Dir2[0]*Dir2[0]+
        Dir2[1]*Dir2[1]+
        Dir2[2]*Dir2[2],0.5);

  if((PtSepVec[0]*Dir1[0]+PtSepVec[1]*Dir1[1]+PtSepVec[2]*Dir1[2])>0)
    Scales[0] =     PtSep / (3.*Dir1Length);
  else
    Scales[0] = 0 - PtSep / (3.*Dir1Length);

  if((PtSepVec[0]*Dir2[0]+PtSepVec[1]*Dir2[1]+PtSepVec[2]*Dir2[2])>0)
    Scales[1] = 0 - PtSep / (3.*Dir2Length);
  else
    Scales[1] =     PtSep / (3.*Dir2Length);

}


//----------------------------------------------------------------------
// Interpolate point between two seeds and return as double[3]
//

void BezierCurveHelper::GetBezierPointXYZ(recob::Seed const& s1, recob::Seed const& s2, float s, double * xyz)
{
  TVector3 BezierPoint = GetBezierPoint(s1, s2, s);
  xyz[0]=BezierPoint[0];
  xyz[1]=BezierPoint[1];
  xyz[2]=BezierPoint[2];

}


//----------------------------------------------------------------------
// Interpolate point between two seeds and return as a TVector3
//

TVector3 BezierCurveHelper::GetBezierPoint(recob::Seed const& s1, recob::Seed const& s2, float s)
{
  return GetBezierPointCubic(s1, s2, s);
}


//-----------------------------------------------------------------------

TVector3 BezierCurveHelper::GetBezierPointQuartic(recob::Seed const& s1, recob::Seed const& s2, float s)
{
  TVector3 ReturnVec3;

  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Mid3[3], Dummy[3];
  s1.GetDirection(Dir1, Dummy);
  s2.GetDirection(Dir2, Dummy);
  s1.GetPoint(Pt1,      Dummy);
  s2.GetPoint(Pt2,      Dummy);

  if(s<=0.) return TVector3(Pt1[0],Pt1[1],Pt1[2]);
  if(s>=1.) return TVector3(Pt2[0],Pt2[1],Pt2[2]);

  TVector3 pt1(Pt1[0],Pt1[1],Pt1[2]);
  TVector3 pt2(Pt2[0],Pt2[1],Pt2[2]);
  TVector3 dir1(Dir1[0],Dir1[1],Dir1[2]);
  TVector3 dir2(Dir2[0],Dir2[1],Dir2[2]);

  double lambda1 = (pt1-pt2).Dot(dir1-dir2)/(pow(dir2.Dot(dir1),2)-dir2.Mag2()*dir1.Mag2());
  double lambda2 = (pt2-pt1).Dot(dir2-dir1)/(pow(dir1.Dot(dir2),2)-dir1.Mag2()*dir2.Mag2());


  int Sign1, Sign2;
  
  if(lambda1>0) Sign1=1;
  else Sign1=-1;

  if(lambda2>0) Sign2=1;
  else Sign2=-1;

  double ns=1.-s;
  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i] + Sign1*Dir1[i];
      Mid2[i]=0.5*(Pt1[i] + lambda1 * Dir1[i] + Pt2[i] +lambda2 * Dir2[i]);
      Mid3[i]=Pt2[i] - Sign2*Dir2[i];
      ReturnVec3[i]=
	ns*ns*ns*ns        * Pt1[i]
	+ 4.*ns*ns*ns*s    * Mid1[i]
	+ 6.*ns*ns*s*s     * Mid2[i]
	+ 4.*ns*s*s*s      * Mid3[i]
	+ s*s*s*s          * Pt2[i];
    }
  
  return ReturnVec3;
  
}


//-----------------------------------------------------------------------

TVector3 BezierCurveHelper::GetBezierPointCubic(recob::Seed const& s1, recob::Seed const& s2, float s)
{
  TVector3 ReturnVec3;

  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Dummy[3];
  s1.GetDirection(Dir1, Dummy);
  s2.GetDirection(Dir2, Dummy);
  s1.GetPoint(Pt1,      Dummy);
  s2.GetPoint(Pt2,      Dummy);

  if(s<=0.) return TVector3(Pt1[0],Pt1[1],Pt1[2]);
  if(s>=1.) return TVector3(Pt2[0],Pt2[1],Pt2[2]);
  

  double DirScales[2];
  GetDirectionScales(Pt1,Pt2,Dir1,Dir2,DirScales);

  double ns=1.-s;
  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      ReturnVec3[i]=
	ns * ns * ns *        Pt1[i]
	+ 3.* ns * ns * s *   Mid1[i]
	+ 3.* ns * s * s *    Mid2[i]
	+ s * s * s *         Pt2[i];
    }


  return ReturnVec3;



}



//----------------------------------------------------------------------
// Interpolate a whole vector of points between two seeds - more
//  efficient than doing each one separately
//


std::vector<TVector3> BezierCurveHelper::GetBezierPoints(recob::Seed const& s1, recob::Seed const& s2, int N)
{
  return GetBezierPointsCubic(s1, s2, N);
}


//--------------------------------------------

std::vector<TVector3> BezierCurveHelper::GetBezierPointsCubic(recob::Seed const& s1, recob::Seed const& s2, int N)
{

  std::vector<TVector3> ReturnVec;
  ReturnVec.resize(N);

  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Dummy[3];
  s1.GetDirection(Dir1, Dummy);
  s2.GetDirection(Dir2, Dummy);
  s1.GetPoint(Pt1,      Dummy);
  s2.GetPoint(Pt2,      Dummy);

  double DirScales[2];
  GetDirectionScales(Pt1,Pt2,Dir1,Dir2,DirScales);


  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      for(int p=0; p!=N; p++)
	{
	  double t  = float(p) / float(N-1);
	  double nt = 1.-t;
	  ReturnVec.at(p)[i] =
	    nt*nt*nt        * Pt1[i]
	    + 3.*nt*nt*t    * Mid1[i]
	    + 3.*nt*t*t     * Mid2[i]
	    + t*t*t         * Pt2[i];
	}
    }
  
  return ReturnVec;
}


//-------------------------------------------------

std::vector<TVector3> BezierCurveHelper::GetBezierPointsQuartic(recob::Seed const& s1, recob::Seed const& s2, int N)
{

  std::vector<TVector3> ReturnVec;
  ReturnVec.resize(N);

  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Mid3[3], Dummy[3];
  s1.GetDirection(Dir1, Dummy);
  s2.GetDirection(Dir2, Dummy);
  s1.GetPoint(Pt1,      Dummy);
  s2.GetPoint(Pt2,      Dummy);

  TVector3 pt1(Pt1[0],Pt1[1],Pt1[2]);
  TVector3 pt2(Pt2[0],Pt2[1],Pt2[2]);
  TVector3 dir1(Dir1[0],Dir1[1],Dir1[2]);
  TVector3 dir2(Dir2[0],Dir2[1],Dir2[2]);
  

  double lambda1 = (pt1-pt2).Dot(dir1-dir2)/(pow(dir2.Dot(dir1),2)-dir2.Mag2()*dir1.Mag2());
  double lambda2 = (pt2-pt1).Dot(dir2-dir1)/(pow(dir1.Dot(dir2),2)-dir1.Mag2()*dir2.Mag2());
   
  int Sign1, Sign2;
  
  if(lambda1>0) Sign1=1;
  else Sign1=-1;

  if(lambda2>0) Sign2=1;
  else Sign2=-1;

  for(int i=0; i!=3; i++)
    {
     
         
      Mid1[i]=Pt1[i] + Sign1*Dir1[i];
      Mid2[i]=0.5*(Pt1[i] + lambda1 * Dir1[i] + Pt2[i] +lambda2 * Dir2[i]);
      Mid3[i]=Pt2[i] - Sign2*Dir2[i];
      for(int p=0; p!=N; p++)
	{
	  double t  = float(p) / float(N-1);
	  double nt = 1.-t;
	  ReturnVec.at(p)[i] =
	    nt*nt*nt*nt        * Pt1[i]
	    + 4.*nt*nt*nt*t    * Mid1[i]
	    + 6.*nt*nt*t*t     * Mid2[i]
	    + 4.*nt*t*t*t      * Mid3[i]
	    + t*t*t*t          * Pt2[i];
	  
	}
    }
  
  return ReturnVec;
}


}
