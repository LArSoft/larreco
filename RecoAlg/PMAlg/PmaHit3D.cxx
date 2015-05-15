/**
 *  @file   PmaHit3D.h
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 */

#include "RecoAlg/PMAlg/PmaHit3D.h"

#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

pma::Hit3D::Hit3D(void) :
	fTPC(0), fPlane(0),
	fPoint3D(0, 0, 0),
	fPoint2D(0, 0), fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fEnabled(true), fOutlier(false)
{
}

pma::Hit3D::Hit3D(const recob::Hit& src) :
	recob::Hit(src),
	fPoint3D(0, 0, 0),
	fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fEnabled(true), fOutlier(false)
{
	fTPC = src.WireID().TPC;
	fPlane = src.WireID().Plane;

	art::ServiceHandle<geo::Geometry> geom;
	art::ServiceHandle<util::DetectorProperties> detprop;

	double wpitch = geom->TPC(fTPC).Plane(fPlane).WirePitch();
	double dpitch = fabs(detprop->GetXTicksCoefficient(fTPC, 0));

	fPoint2D.Set(wpitch * src.WireID().Wire, dpitch * src.PeakTime());
}

pma::Hit3D::Hit3D(const pma::Hit3D& src) :
	recob::Hit(src),
	fTPC(src.fTPC), fPlane(src.fPlane),
	fPoint3D(src.fPoint3D),
	fPoint2D(src.fPoint2D), fProjection2D(src.fProjection2D),
	fSegFraction(src.fSegFraction), fSigmaFactor(src.fSigmaFactor),
	fEnabled(src.fEnabled), fOutlier(src.fOutlier)
{
}

