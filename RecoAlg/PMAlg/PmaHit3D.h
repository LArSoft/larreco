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

#include "RecoBase/Hit.h"

#include "Geometry/Geometry.h"

#include <functional>

#include "TVector2.h"
#include "TVector3.h"

namespace pma
{
	class Hit3D;
	struct bTrajectory3DOrderLess;

	class Track3D;
	
	double Dist2(const TVector2& v1, const TVector2& v2);
	double Dist2(const TVector3& v1, const TVector3& v2);
	size_t GetHitsCount(const std::vector< pma::Hit3D* >& hits, unsigned int view);
	double GetSummedADC(const std::vector< pma::Hit3D* >& hits, unsigned int view = geo::kUnknown);
	double GetSummedAmpl(const std::vector< pma::Hit3D* >& hits, unsigned int view = geo::kUnknown);

	double GetHitsRadius(const std::vector< pma::Hit3D* >& hits, bool exact = false);
}

class pma::Hit3D
{
	friend class Track3D;
	friend struct bTrajectory3DOrderLess;

public:
	Hit3D(void);
	Hit3D(const recob::Hit& src);
	Hit3D(const pma::Hit3D& src);
	virtual ~Hit3D(void) {}

	TVector3 const & Point3D(void) const { return fPoint3D; }

	void SetPoint3D(const TVector3& p3d) { fPoint3D = p3d; }

	TVector2 const & Point2D(void) const { return fPoint2D; }
	TVector2 const & Projection2D(void) const { return fProjection2D; }

	unsigned int TPC(void) const { return fTPC; }
	unsigned int View2D(void) const { return fPlane; }
	unsigned int Wire(void) const { return fWire; }
	float PeakTime(void) const { return fHit.PeakTime(); }

	float SummedADC(void) const { return fHit.SummedADC(); }
	float GetAmplitude(void) const { return fHit.PeakAmplitude(); }
	float GetSigmaFactor(void) const { return fSigmaFactor; }
	void SetSigmaFactor(float value) { fSigmaFactor = value; }

	double GetDistToProj(void) const { return sqrt(GetDist2ToProj()); }
	double GetDist2ToProj(void) const { return pma::Dist2(fPoint2D, fProjection2D); }

	void SetProjection(const TVector2& p, float b)
	{
		fProjection2D.Set(p); fSegFraction = b;
	}

	bool IsEnabled(void) const { return (fEnabled && !fOutlier); }
	void SetEnabled(bool state) { fEnabled = state; }

	bool IsOutlier(void) const { return fOutlier; }
	virtual void TagOutlier(bool state) { fOutlier = state; }

private:

	recob::Hit fHit;  // source 2D hit

	unsigned int fTPC, fPlane, fWire;

	TVector3 fPoint3D;       // hit position in 3D space
	TVector2 fPoint2D;       // hit position in 2D wire view, scaled to [cm]
	TVector2 fProjection2D;  // projection to polygonal line in 2D wire view, scaled to [cm]
	float fSegFraction;      // segment fraction set by projection
	float fSigmaFactor;      // impact factor on the objective function

	bool fEnabled; // used or not in the optimisation - due to various reasons
	bool fOutlier; // tagged as not really hit of this track (like delta ray)

};

struct pma::bTrajectory3DOrderLess :
	public std::binary_function<pma::Hit3D*, pma::Hit3D*, bool>
{
	bool operator() (pma::Hit3D* h1, pma::Hit3D* h2)
	{
		if (h1 && h2) return h1->fSegFraction < h2->fSegFraction;
		else return false;
	}
};

#endif

