////////////////////////////////////////////////////////////////////////
///
/// \file   BezierCurveHelper.h
///
/// \brief  Service for finding 3D track seeds
///
/// \author B J P Jones
///
/// Helper object for interpolating betweeen bezier
/// curve points, specified by seed objects
///


#ifndef BEZIERCURVEHELPER_H
#define BEZIERCURVEHELPER_H

#include <vector>

#include "lardataobj/RecoBase/Seed.h"

#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"

namespace trkf {



  class BezierCurveHelper {
  public:

    // Constructor.
    BezierCurveHelper();
    explicit BezierCurveHelper(int fCurveRes);

    // Destructor.
    ~BezierCurveHelper();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    std::vector<TVector3> GetBezierPoints(recob::Seed const& s1, recob::Seed const& s2, int N=100);
    std::vector<TVector3> GetBezierPointsQuartic(recob::Seed const& s1, recob::Seed const& s2, int N=100);
    std::vector<TVector3> GetBezierPointsCubic(recob::Seed const& s1, recob::Seed const& s2, int N=100);
    double   GetSegmentLength(recob::Seed  const& s1, recob::Seed const& s2);
    void     GetBezierPointXYZ(recob::Seed const& s1, recob::Seed const& s2, float t, double * xyz);
    TVector3 GetBezierPoint(recob::Seed const& s1, recob::Seed const& s2, float t);
    TVector3 GetBezierPointCubic(recob::Seed const& s1, recob::Seed const& s2, float t);
    TVector3 GetBezierPointQuartic(recob::Seed const& s1, recob::Seed const& s2, float t);
    void     GetDirectionScales(double * Pt1, double * Pt2, double * Dir1, double * Dir2, double *Scales);


    void SetCurveResolution(int CurveRes) {fCurveResolution=CurveRes;}
    int  GetCurveResolution()             {return fCurveResolution;}

  private:
    int fCurveResolution;


  };



}

#endif
