/**
 *  @file   PmaHit3D.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Hit 3D wrapped around recob::Hit. Adds support for PMA optimizations.
 *          See PmaTrack3D.h file for details.
 */

#ifndef PmaHit3D_h
#define PmaHit3D_h

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
namespace detinfo {
  class DetectorPropertiesData;
}

#include <cmath>

#include "TVector2.h"
#include "TVector3.h"

namespace pma {
  class Hit3D;
  class Track3D;
}

class pma::Hit3D {
  friend class Track3D;
  friend struct bTrajectory3DOrderLess;

public:
  Hit3D();
  Hit3D(detinfo::DetectorPropertiesData const& detProp, art::Ptr<recob::Hit> src);
  Hit3D(detinfo::DetectorPropertiesData const& detProp,
        unsigned int wire,
        unsigned int view,
        unsigned int tpc,
        unsigned int cryo,
        float peaktime,
        float ampl,
        float area);
  Hit3D(const pma::Hit3D& src);

  art::Ptr<recob::Hit> const&
  Hit2DPtr() const
  {
    return fHit;
  }

  TVector3 const&
  Point3D() const
  {
    return fPoint3D;
  }

  void
  SetPoint3D(const TVector3& p3d)
  {
    fPoint3D = p3d;
  }
  void
  SetPoint3D(double x, double y, double z)
  {
    fPoint3D.SetXYZ(x, y, z);
  }

  TVector2 const&
  Point2D() const noexcept
  {
    return fPoint2D;
  }
  TVector2 const&
  Projection2D() const noexcept
  {
    return fProjection2D;
  }

  unsigned int
  Cryo() const noexcept
  {
    return fCryo;
  }
  unsigned int
  TPC() const noexcept
  {
    return fTPC;
  }
  unsigned int
  View2D() const noexcept
  {
    return fPlane;
  }
  unsigned int
  Wire() const noexcept
  {
    return fWire;
  }
  float
  PeakTime() const noexcept
  {
    return fPeakTime;
  }

  float
  SummedADC() const noexcept
  {
    return fArea;
  }
  float
  GetAmplitude() const noexcept
  {
    return fAmpl;
  }
  float
  GetSigmaFactor() const noexcept
  {
    return fSigmaFactor;
  }
  void
  SetSigmaFactor(float value) noexcept
  {
    fSigmaFactor = value;
  }

  double
  Dx() const noexcept
  {
    return fDx;
  }

  double
  GetDistToProj() const
  {
    return sqrt(GetDist2ToProj());
  }
  double GetDist2ToProj() const;

  float
  GetSegFraction() const noexcept
  {
    return fSegFraction;
  }
  void
  SetProjection(const TVector2& p, float b)
  {
    fProjection2D.Set(p);
    fSegFraction = b;
  }
  void
  SetProjection(double x, double y, float b)
  {
    fProjection2D.Set(x, y);
    fSegFraction = b;
  }

  bool
  IsEnabled() const noexcept
  {
    return (fEnabled && !fOutlier);
  }
  void
  SetEnabled(bool state) noexcept
  {
    fEnabled = state;
  }

  bool
  IsOutlier() const noexcept
  {
    return fOutlier;
  }
  void
  TagOutlier(bool state) noexcept
  {
    fOutlier = state;
  }

private:
  art::Ptr<recob::Hit> fHit; // source 2D hit

  unsigned int fCryo, fTPC, fPlane, fWire;
  float fPeakTime, fAmpl, fArea;

  TVector3 fPoint3D;      // hit position in 3D space
  TVector2 fPoint2D;      // hit position in 2D wire view, scaled to [cm]
  TVector2 fProjection2D; // projection to polygonal line in 2D wire view, scaled to [cm]
  float fSegFraction;     // segment fraction set by the projection
  float fSigmaFactor;     // impact factor on the objective function

  double fDx; // dx seen by corresponding 2D hit, set during dQ/dx sequece calculation

  bool fEnabled; // used or not in the optimisation - due to various reasons
  bool fOutlier; // tagged as a not really hit of this track (like delta ray)

  pma::Track3D* fParent; // track which contains this hit
};

#endif
