/**
 *  @file   PmaHit3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Hit 3D wrapped around recob::Hit. Adds support for PMA optimizations.
 */

#include "RecoAlg/PMAlg/PmaHit3D.h"

#include "Utilities/DetectorProperties.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

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

size_t pma::GetHitsCount(const std::vector< pma::Hit3D >& hits, unsigned int view)
{
	size_t n = 0;
	for (auto const& hit : hits)
		if ((view == geo::kUnknown) || (view == hit.View2D())) n++;
	return n;
}

double pma::GetSummedADC(const std::vector< pma::Hit3D >& hits, unsigned int view)
{
	double sum = 0.0;
	for (auto const& hit : hits)
		if ((view == geo::kUnknown) || (view == hit.View2D()))
			sum += hit.SummedADC();
	return sum;
}

double pma::GetSummedAmpl(const std::vector< pma::Hit3D >& hits, unsigned int view)
{
	double sum = 0.0;
	for (auto const& hit : hits)
		if ((view == geo::kUnknown) || (view == hit.View2D()))
			sum += hit.PeakAmplitude();
	return sum;
}

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

