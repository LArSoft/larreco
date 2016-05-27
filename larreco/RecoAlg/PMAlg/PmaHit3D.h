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

#include "lardata/RecoBase/Hit.h"

#include "larcore/Geometry/Geometry.h"

#include <functional>

#include "TVector2.h"
#include "TVector3.h"

namespace pma
{
	class Hit3D;
	struct bTrajectory3DOrderLess;

	class Track3D;
}

class pma::Hit3D
{
	friend class Track3D;
	friend struct bTrajectory3DOrderLess;

public:
	Hit3D(void);
	Hit3D(art::Ptr< recob::Hit > src);
	Hit3D(unsigned int wire, unsigned int view, unsigned int tpc, unsigned int cryo,
		float peaktime, float ampl, float area);
	Hit3D(const pma::Hit3D& src);
	virtual ~Hit3D(void) {}

	art::Ptr< recob::Hit > Hit2DPtr(void) const { return fHit; }

	TVector3 const & Point3D(void) const { return fPoint3D; }

	void SetPoint3D(const TVector3& p3d) { fPoint3D = p3d; }
	void SetPoint3D(double x, double y, double z) { fPoint3D.SetXYZ(x, y, z); }

	TVector2 const & Point2D(void) const { return fPoint2D; }
	TVector2 const & Projection2D(void) const { return fProjection2D; }

	unsigned int Cryo(void) const { return fCryo; }
	unsigned int TPC(void) const { return fTPC; }
	unsigned int View2D(void) const { return fPlane; }
	unsigned int Wire(void) const { return fWire; }
	float PeakTime(void) const { return fPeakTime; }

	float SummedADC(void) const { return fArea; }
	float GetAmplitude(void) const { return fAmpl; }
	float GetSigmaFactor(void) const { return fSigmaFactor; }
	void SetSigmaFactor(float value) { fSigmaFactor = value; }

	double Dx(void) const { return fDx; }

	double GetDistToProj(void) const { return sqrt(GetDist2ToProj()); }
	double GetDist2ToProj(void) const;

	float GetSegFraction() const { return fSegFraction; }
	void SetProjection(const TVector2& p, float b)
	{
		fProjection2D.Set(p); fSegFraction = b;
	}
	void SetProjection(double x, double y, float b)
	{
		fProjection2D.Set(x, y); fSegFraction = b;
	}

	bool IsEnabled(void) const { return (fEnabled && !fOutlier); }
	void SetEnabled(bool state) { fEnabled = state; }

	bool IsOutlier(void) const { return fOutlier; }
	virtual void TagOutlier(bool state) { fOutlier = state; }

private:

	art::Ptr< recob::Hit > fHit;  // source 2D hit

	unsigned int fCryo, fTPC, fPlane, fWire;
	float fPeakTime, fAmpl, fArea;

	TVector3 fPoint3D;       // hit position in 3D space
	TVector2 fPoint2D;       // hit position in 2D wire view, scaled to [cm]
	TVector2 fProjection2D;  // projection to polygonal line in 2D wire view, scaled to [cm]
	float fSegFraction;      // segment fraction set by the projection
	float fSigmaFactor;      // impact factor on the objective function

	double fDx;    // dx seen by corresponding 2D hit, set during dQ/dx sequece calculation

	bool fEnabled; // used or not in the optimisation - due to various reasons
	bool fOutlier; // tagged as a not really hit of this track (like delta ray)

	pma::Track3D* fParent;   // track which contains this hit
};

#endif

