/**
 *  @file   PmaSegment3D.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track segment. See PmaTrack3D.h file for details.
 */

#ifndef PmaSegment3D_h
#define PmaSegment3D_h

#include "larreco/RecoAlg/PMAlg/PmaElement3D.h"
#include "larreco/RecoAlg/PMAlg/PmaNode3D.h"
#include "larreco/RecoAlg/PMAlg/SortedObjects.h"

class TVector2;
#include "TVector3.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

namespace pma {
  class Segment3D;
  class Track3D; // only declare here to keep "parent" of Segment3D
}

class pma::Segment3D : public pma::Element3D, public pma::SortedObjectBase {
public:
  Segment3D(void) : fParent(0) {}
  Segment3D(pma::Track3D* trk, pma::Node3D* vstart, pma::Node3D* vstop);

  Vector3D Start(void) const
  {
    auto const& p = static_cast<Node3D*>(Prev())->Point3D();
    return Vector3D(p.X(), p.Y(), p.Z());
  }
  Vector3D End(void) const
  {
    auto const& p = static_cast<Node3D*>(Next())->Point3D();
    return Vector3D(p.X(), p.Y(), p.Z());
  }

  /// Distance [cm] from the 3D segment to the point 3D.
  double GetDistance2To(const TVector3& p3d) const override;

  /// Distance [cm] from the 2D point to the object's 2D projection in one of wire views.
  double GetDistance2To(const TVector2& p2d, unsigned int view) const override;

  /// Get 3D direction cosines of this segment.
  pma::Vector3D GetDirection3D(void) const override;

  /// Get 3D projection of a 2D point from the view.
  TVector3 GetProjection(const TVector2& p, unsigned int view) const;

  /// Get 3D projection of a 2D point from the view, no limitations if it falls beyond
  /// the segment endpoints.
  TVector3 GetUnconstrainedProj3D(const TVector2& p2d, unsigned int view) const override;

  /// Set hit 3D position and its 2D projection to the vertex.
  void SetProjection(pma::Hit3D& h) const override;

  /// Squared sum of half-lengths of connected 3D segments
  /// (used in the vertex position optimization).
  double Length2(void) const override;

  pma::Track3D* Parent(void) const { return fParent; }

private:
  Segment3D(const pma::Segment3D& src);

  double SumDist2Hits(void) const override;

  pma::Track3D* fParent;

  static double GetDist2(const TVector3& psrc, const TVector3& p0, const TVector3& p1);
  static double GetDist2(const TVector2& psrc, const TVector2& p0, const TVector2& p1);
};

#endif
