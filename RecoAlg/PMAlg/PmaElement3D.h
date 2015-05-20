/**
 *  @file   PmaElement3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Base for 3D segments and track nodes. See PmaTrack3D.h file for details.
 *          See PmaTrack3D.h file for details.
 */

#ifndef PmaElement3D_h
#define PmaElement3D_h

#include "RecoAlg/PMAlg/PmaHit3D.h"

#include "Geometry/Geometry.h"

#include "TVector2.h"
#include "TVector3.h"

namespace pma
{
	class Element3D;

	class Track3D;
}

class pma::Element3D
{
public:
	Element3D(void);
	virtual ~Element3D(void) {}

	/// Distance [cm] from the 3D point to the object 3D.
	virtual double GetDistance2To(const TVector3& p3d) const = 0;

	/// Distance [cm] from the 2D point to the object's 2D projection in one of wire views.
	virtual double GetDistance2To(const TVector2& p2d, unsigned int view) const = 0;

	virtual void SetProjection(pma::Hit3D& h) const = 0;

	virtual double Length2(void) const = 0;
	double Length(void) const { return sqrt(Length2()); }

	pma::Hit3D& Hit(size_t index) { return *(fAssignedHits[index]); }
	void AddHit(pma::Hit3D* h)
	{
		fAssignedHits.push_back(h);
		SetProjection(*h);
	}

	size_t NHits(void) const { return fAssignedHits.size(); }
	size_t NEnabledHits(unsigned int view = geo::kUnknown) const;

	TVector3 const & ReferencePoint(size_t index) const { return *(fAssignedPoints[index]); }
	size_t NPoints(void) const { return fAssignedPoints.size(); }
	void AddPoint(TVector3* p) { fAssignedPoints.push_back(p); }

	/// Clear hits/points vectors of this element, optionally only
	/// those which are owned by given track.
	virtual void ClearAssigned(pma::Track3D* trk = 0);

	void UpdateHitParams(void);
	void UpdateProjection(void);
	void SortHits(void);

	double SumDist2(void) const;
	double SumHitsQ(unsigned int view) const { return fSumHitsQ[view]; }
	unsigned int NHits(unsigned int view) const { return fNHits[view]; }
	unsigned int NThisHits(unsigned int view) const { return fNThisHits[view]; }

	double HitsRadius3D(unsigned int view) const;

	/// Check if the vertex 3D position is fixed.
	bool IsFrozen(void) const { return fFrozen; }
	/// Fix / relese vertex 3D position.
	void SetFrozen(bool state) { fFrozen = state; }

	static float OptFactor(unsigned int view) { return fOptFactors[view]; }
	static void SetOptFactor(unsigned int view, float value) { fOptFactors[view] = value; }

protected:
	bool fFrozen;
	std::vector< pma::Hit3D* > fAssignedHits;  // 2D hits
	std::vector< TVector3* > fAssignedPoints;  // 3D peculiar points reconstructed elsewhere
	unsigned int fNThisHits[3];
	unsigned int fNHits[3];
	double fSumHitsQ[3];
	double fHitsRadius;

	static float fOptFactors[3]; // impact factors of data from various 2D views
};

#endif

