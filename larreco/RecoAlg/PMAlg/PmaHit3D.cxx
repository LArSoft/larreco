/**
 *  @file   PmaHit3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Hit 3D wrapped around recob::Hit. Support for PMA optimizations.
 *          See PmaTrack3D.h file for details.
 */

#include "larreco/RecoAlg/PMAlg/PmaHit3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"


pma::Hit3D::Hit3D(void) :
	fCryo(0), fTPC(0), fPlane(0), fWire(0),
	fPeakTime(0), fAmpl(0), fArea(0),
	fPoint3D(0, 0, 0),
	fPoint2D(0, 0), fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fDx(0),
	fEnabled(true), fOutlier(false),
	fParent(0)
{
}

pma::Hit3D::Hit3D(art::Ptr< recob::Hit > src) :
	fHit(src),
	fPoint3D(0, 0, 0),
	fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fDx(0),
	fEnabled(true), fOutlier(false),
	fParent(0) // set only when pushed to track
{
	fCryo = src->WireID().Cryostat;
	fTPC = src->WireID().TPC;
	fPlane = src->WireID().Plane;
	fWire = src->WireID().Wire;

	fPeakTime = src->PeakTime();

	fAmpl = src->PeakAmplitude();
	fArea = src->SummedADC();

	fPoint2D = pma::WireDriftToCm(fWire, fPeakTime, fPlane, fTPC, fCryo);
}

pma::Hit3D::Hit3D(unsigned int wire, unsigned int view, unsigned int tpc, unsigned int cryo,
	float peaktime, float ampl, float area) :
	fPoint3D(0, 0, 0),
	fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fDx(0),
	fEnabled(false), fOutlier(false),
	fParent(0) // set only when pushed to track
{
	fCryo = cryo;
	fTPC = tpc;
	fPlane = view;
	fWire = wire;

	fPeakTime = peaktime;

	fAmpl = ampl;
	fArea = area;

	fPoint2D = pma::WireDriftToCm(fWire, fPeakTime, fPlane, fTPC, fCryo);
}

pma::Hit3D::Hit3D(const pma::Hit3D& src) :
	fHit(src.fHit),
	fCryo(src.fCryo), fTPC(src.fTPC), fPlane(src.fPlane), fWire(src.fWire),
	fPeakTime(src.fPeakTime), fAmpl(src.fAmpl), fArea(src.fArea),
	fPoint3D(src.fPoint3D),
	fPoint2D(src.fPoint2D), fProjection2D(src.fProjection2D),
	fSegFraction(src.fSegFraction), fSigmaFactor(src.fSigmaFactor),
	fDx(src.fDx),
	fEnabled(src.fEnabled), fOutlier(src.fOutlier),
	fParent(0) // set only when pushed to track
{
}

double pma::Hit3D::GetDist2ToProj(void) const
{
	return pma::Dist2(fPoint2D, fProjection2D);
}

