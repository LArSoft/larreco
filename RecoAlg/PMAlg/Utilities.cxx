/**
 *  @file   Utilities.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Some geometrical functions and sorting helpers.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/Utilities.h"
#include "RecoAlg/PMAlg/PmaHit3D.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

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

double pma::GetSegmentProjVector(const TVector2& p, const TVector2& p0, const TVector2& p1)
{
	TVector2 v0(p); v0 -= p0;
	TVector2 v1(p1); v1 -= p0;

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();
	double mag = v0Norm * v1Norm;
	double cosine = 0.0;
	if (mag != 0.0) cosine = v0 * v1 / mag;

	return v0Norm * cosine / v1Norm;
}

double pma::GetSegmentProjVector(const TVector3& p, const TVector3& p0, const TVector3& p1)
{
	TVector3 v0(p); v0 -= p0;
	TVector3 v1(p1); v1 -= p0;

	double v0Norm = v0.Mag();
	double v1Norm = v1.Mag();
	double mag = v0Norm * v1Norm;
	double cosine = 0.0;
	if (mag != 0.0) cosine = v0 * v1 / mag;

	return v0Norm * cosine / v1Norm;
}

TVector2 pma::GetProjectionToSegment(const TVector2& p, const TVector2& p0, const TVector2& p1)
{
	TVector2 v1(p1); v1 -= p0;

	double b = GetSegmentProjVector(p, p0, p1);
	TVector2 r(p0);
	r += (v1 * b);
	return r;
}

TVector3 pma::GetProjectionToSegment(const TVector3& p, const TVector3& p0, const TVector3& p1)
{
	TVector3 v1(p1); v1 -= p0;

	double b = GetSegmentProjVector(p, p0, p1);
	TVector3 r(p0);
	r += (v1 * b);
	return r;
}

bool pma::bTrajectory3DOrderLess::operator() (pma::Hit3D* h1, pma::Hit3D* h2)
{
	if (h1 && h2) return h1->fSegFraction < h2->fSegFraction;
	else return false;
}

bool pma::bTrajectory3DDistLess::operator() (pma::Hit3D* h1, pma::Hit3D* h2)
{
	if (h1 && h2) return h1->GetDist2ToProj() < h2->GetDist2ToProj();
	else return false;
}

pma::bSegmentProjLess::bSegmentProjLess(const TVector3& s0, const TVector3& s1) :
	segStart(s0), segStop(s1)
{
	if (s0 == s1) mf::LogError("pma::bSegmentProjLess") << "Vectors equal!";
}

