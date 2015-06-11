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

#include "RecoAlg/PMAlg/PmaHit3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include "Utilities/DetectorProperties.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"


pma::Hit3D::Hit3D(void) :
	fTPC(0), fPlane(0), fWire(0),
	fPoint3D(0, 0, 0),
	fPoint2D(0, 0), fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fEnabled(true), fOutlier(false)
{
}

pma::Hit3D::Hit3D(art::Ptr< recob::Hit > src) :
	fHit(src),
	fPoint3D(0, 0, 0),
	fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fEnabled(true), fOutlier(false)
{
	fTPC = src->WireID().TPC;
	fPlane = src->WireID().Plane;
	fWire = src->WireID().Wire;

	fPoint2D = pma::WireDriftToCm(fWire, src->PeakTime(), fPlane, fTPC, src->WireID().Cryostat);
}

pma::Hit3D::Hit3D(const pma::Hit3D& src) :
	fHit(src.fHit),
	fTPC(src.fTPC), fPlane(src.fPlane), fWire(src.fWire),
	fPoint3D(src.fPoint3D),
	fPoint2D(src.fPoint2D), fProjection2D(src.fProjection2D),
	fSegFraction(src.fSegFraction), fSigmaFactor(src.fSigmaFactor),
	fEnabled(src.fEnabled), fOutlier(src.fOutlier)
{
}

double pma::Hit3D::GetDist2ToProj(void) const
{
	return pma::Dist2(fPoint2D, fProjection2D);
}

