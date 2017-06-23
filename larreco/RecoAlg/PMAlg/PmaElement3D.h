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

#include "larreco/RecoAlg/PMAlg/PmaHit3D.h"

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
	/// TPC index or -1 if out of any TPC.
	int TPC(void) const { return fTPC; }
	/// Cryostat index or -1 if out of any cryostat.
	int Cryo(void) const { return fCryo; }

	/// Distance [cm] from the 3D point to the object 3D.
	virtual double GetDistance2To(const TVector3& p3d) const = 0;

	/// Distance [cm] from the 2D point to the object's 2D projection in one of wire views.
	virtual double GetDistance2To(const TVector2& p2d, unsigned int view) const = 0;

	/// Get 3D direction cosines corresponding to this element.
	virtual pma::Vector3D GetDirection3D(void) const = 0;

	virtual TVector3 GetUnconstrainedProj3D(const TVector2& p2d, unsigned int view) const = 0;

	virtual void SetProjection(pma::Hit3D& h) const = 0;

	virtual double Length2(void) const = 0;
	double Length(void) const { return sqrt(Length2()); }

	const std::vector< pma::Hit3D* > & Hits(void) const { return fAssignedHits; }

    bool HasHit(const pma::Hit3D* h) const
    {
        for (const auto a : fAssignedHits) { if (h == a) return true; }
        return false;
    } 

	pma::Hit3D& Hit(size_t index) { return *(fAssignedHits[index]); }
	void RemoveHitAt(size_t index)
	{
		if (index < fAssignedHits.size())
			fAssignedHits.erase(fAssignedHits.begin() + index);
	}
	void AddHit(pma::Hit3D* h)
	{
		fAssignedHits.push_back(h);
		SetProjection(*h);
	}

	size_t NHits(void) const { return fAssignedHits.size(); }
	size_t NEnabledHits(unsigned int view = geo::kUnknown) const;
	size_t NPrecalcEnabledHits(void) const { return fNThisHitsEnabledAll; }

	TVector3 const & ReferencePoint(size_t index) const { return *(fAssignedPoints[index]); }
	size_t NPoints(void) const { return fAssignedPoints.size(); }
	void AddPoint(TVector3* p) { fAssignedPoints.push_back(p); }

	/// Clear hits/points vectors of this element, optionally only
	/// those which are owned by given track.
	virtual void ClearAssigned(pma::Track3D* trk = 0);

	void UpdateHitParams(void);
	void UpdateProjection(void) { for (auto h : fAssignedHits) SetProjection(*h); }
	void SortHits(void);

	double SumDist2(void) const;
	double SumDist2(unsigned int view) const;
	double SumHitsQ(unsigned int view) const { return fSumHitsQ[view]; }
	unsigned int NHits(unsigned int view) const { return fNHits[view]; }
	unsigned int NThisHits(unsigned int view) const { return fNThisHits[view]; }

	double HitsRadius3D(unsigned int view) const;

	/// Check if the vertex 3D position is fixed.
	bool IsFrozen(void) const { return fFrozen; }
	/// Fix / relese vertex 3D position.
	void SetFrozen(bool state) { fFrozen = state; }

	bool SelectRndHits(size_t nmax_per_view);
	bool SelectAllHits(void);

	static float OptFactor(unsigned int view) { return fOptFactors[view]; }
	static void SetOptFactor(unsigned int view, float value) { fOptFactors[view] = value; }

protected:
	Element3D(void); // Element3D is only a common base for nodes and segments
	int fTPC, fCryo; // -1 if out of any TPC or cryostat

	bool fFrozen;
	std::vector< pma::Hit3D* > fAssignedHits;  // 2D hits
	std::vector< TVector3* > fAssignedPoints;  // 3D peculiar points reconstructed elsewhere
	size_t fNThisHits[3];
	size_t fNThisHitsEnabledAll;
	size_t fNHits[3];
	double fSumHitsQ[3];
	double fHitsRadius;

	static float fOptFactors[3]; // impact factors of data from various 2D views
};

#endif

