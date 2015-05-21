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

size_t pma::GetHitsCount(const std::vector< pma::Hit3D* >& hits, unsigned int view)
{
	size_t n = 0;
	for (auto const& hit : hits)
		if ((view == geo::kUnknown) || (view == hit->View2D())) n++;
	return n;
}

double pma::GetSummedADC(const std::vector< pma::Hit3D* >& hits, unsigned int view)
{
	double sum = 0.0;
	for (auto const& hit : hits)
		if ((view == geo::kUnknown) || (view == hit->View2D()))
			sum += hit->SummedADC();
	return sum;
}

double pma::GetSummedAmpl(const std::vector< pma::Hit3D* >& hits, unsigned int view)
{
	double sum = 0.0;
	for (auto const& hit : hits)
		if ((view == geo::kUnknown) || (view == hit->View2D()))
			sum += hit->GetAmplitude();
	return sum;
}

double pma::GetHitsRadius3D(const std::vector< pma::Hit3D* >& hits, bool exact)
{
	if (!exact && (hits.size() < 5)) return 0.0;

	if (hits.size() == 0) return 0.0;

	TVector3 mean(0, 0, 0);
	for (size_t i = 0; i < hits.size(); i++)
	{
		mean += hits[i]->Point3D();
	}
	mean *= (1.0 / hits.size());

	double r2, max_r2 = pma::Dist2(hits.front()->Point3D(), mean);
	for (size_t i = 1; i < hits.size(); i++)
	{
		r2 = pma::Dist2(hits[i]->Point3D(), mean);
		if (r2 > max_r2) max_r2 = r2;
	}
	return sqrt(max_r2);
}

double pma::GetHitsRadius2D(const std::vector< pma::Hit3D* >& hits, bool exact)
{
	if (!exact && (hits.size() < 5)) return 0.0;

	if (hits.size() == 0) return 0.0;

	TVector2 mean(0, 0);
	for (size_t i = 0; i < hits.size(); i++)
	{
		mean += hits[i]->Point2D();
	}
	mean *= (1.0 / hits.size());

	double r2, max_r2 = pma::Dist2(hits.front()->Point2D(), mean);
	for (size_t i = 1; i < hits.size(); i++)
	{
		r2 = pma::Dist2(hits[i]->Point2D(), mean);
		if (r2 > max_r2) max_r2 = r2;
	}
	return sqrt(max_r2);
}

pma::Hit3D::Hit3D(void) :
	fTPC(0), fPlane(0), fWire(0),
	fPoint3D(0, 0, 0),
	fPoint2D(0, 0), fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fEnabled(true), fOutlier(false)
{
}

pma::Hit3D::Hit3D(const recob::Hit& src) :
	fHit(src),
	fPoint3D(0, 0, 0),
	fProjection2D(0, 0),
	fSegFraction(0), fSigmaFactor(1),
	fEnabled(true), fOutlier(false)
{
	fTPC = src.WireID().TPC;
	fPlane = src.WireID().Plane;
	fWire = src.WireID().Wire;

	art::ServiceHandle<geo::Geometry> geom;
	art::ServiceHandle<util::DetectorProperties> detprop;

	double wpitch = geom->TPC(fTPC).Plane(fPlane).WirePitch();
	double dpitch = fabs(detprop->GetXTicksCoefficient(fTPC, 0));

	fPoint2D.Set(wpitch * fWire, dpitch * src.PeakTime());
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

