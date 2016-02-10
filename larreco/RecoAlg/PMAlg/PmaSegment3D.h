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

#include "RecoAlg/PMAlg/PmaElement3D.h"
#include "RecoAlg/PMAlg/PmaNode3D.h"
#include "RecoAlg/PMAlg/SortedObjects.h"

namespace pma
{
	class Segment3D;
	class Track3D; // only declare here to keep "parent" of Segment3D
}

class pma::Segment3D : public pma::Element3D, public pma::SortedObjectBase
{
public:
	Segment3D(void) : fParent(0) {}
	Segment3D(pma::Track3D* trk, pma::Node3D* vstart, pma::Node3D* vstop);
	virtual ~Segment3D(void) {}

	/// Distance [cm] from the 3D segment to the point 3D.
	virtual double GetDistance2To(const TVector3& p3d) const;

	/// Distance [cm] from the 2D point to the object's 2D projection in one of wire views.
	virtual double GetDistance2To(const TVector2& p2d, unsigned int view) const;

	/// Get 3D direction cosines.
	TVector3 GetDirection3D(void) const;

	/// Get 3D projection of a 2D point from the view.
	TVector3 GetProjection(const TVector2& p, unsigned int view) const;

	/// Get 3D projection of a 2D point from the view, no limitations if it falls beyond
	/// the segment endpoints.
	virtual TVector3 GetUnconstrainedProj3D(const TVector2& p2d, unsigned int view) const;

	/// Set hit 3D position and its 2D projection to the vertex.
	virtual void SetProjection(pma::Hit3D& h) const;

	/// Squared sum of half-lengths of connected 3D segments
	/// (used in the vertex position optimization).
	virtual double Length2(void) const;

	pma::Track3D* Parent(void) const { return fParent; }

private:
	Segment3D(const pma::Segment3D& src);

	pma::Track3D* fParent;

	static double GetDist2(const TVector3& psrc, const TVector3& p0, const TVector3& p1);
	static double GetDist2(const TVector2& psrc, const TVector2& p0, const TVector2& p1);

};

#endif

