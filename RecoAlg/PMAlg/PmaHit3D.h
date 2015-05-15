/**
 *  @file   PmaHit3D.h
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 */

#ifndef PmaHit3D_h
#define PmaHit3D_h

#include "RecoBase/Hit.h"

#include "TVector2.h"
#include "TVector3.h"

namespace pma
{
	double Dist2(const TVector2& v1, const TVector2& v2);
	double Dist2(const TVector3& v1, const TVector3& v2);

	class Hit3D;
	class Track3D;
}

double pma::Dist2(const TVector2& v1, const TVector2& v2)
{
	double dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y();
	return dx * dx + dy * dy;
}

double pma::Dist2(const TVector3& v1, const TVector3& v2)
{
	double dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y(), dz = v1.Z() - v2.Z();
	return dx * dx + dy * dy + dz * dz;
}

class pma::Hit3D : public recob::Hit
{
	friend class Track3D;

public:
	Hit3D(void);
	Hit3D(const recob::Hit& src);
	Hit3D(const pma::Hit3D& src);

	TVector3 const & Point3D(void) const { return fPoint3D; }

	TVector2 const & Point2D(void) const { return fPoint2D; }
	TVector2 const & Projection2D(void) const { return fProjection2D; }

	unsigned int TPC(void) const { return fTPC; }
	unsigned int View2D(void) const { return fPlane; }

	float GetAmplitude(void) const { return PeakAmplitude(); }
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

	unsigned int fTPC, fPlane;

	TVector3 fPoint3D;       // hit position in 3D space
	TVector2 fPoint2D;       // hit position in 2D wire view, scaled to [cm]
	TVector2 fProjection2D;  // projection to polygonal line in 2D wire view, scaled to [cm]
	float fSegFraction;      // segment fraction set by projection
	float fSigmaFactor;      // impact factor on the objective function

	bool fEnabled; // used or not in the optimisation - due to various reasons
	bool fOutlier; // tagged as not really hit of this track (like delta ray)

};

#endif
