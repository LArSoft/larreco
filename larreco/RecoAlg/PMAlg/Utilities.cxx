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

#include "DetectorInfoServices/DetectorPropertiesService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TVectorT.h"
#include "TMatrixT.h"

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

double pma::SolveLeastSquares3D(const std::vector< std::pair<TVector3, TVector3> >& lines, TVector3& result)
{
	// RS: please, ask me if you need examples/explanation of formulas as they
	// are not easy to derive from the code solely; I have Mathcad sources that
	// were used to test the solving method, weighting, etc.

	result.SetXYZ(0., 0., 0.);
	if (lines.size() < 2)
	{
		mf::LogError("pma::SolveLeastSquares3D") << "Need min. two lines.";
		return -1.0;
	}

	double m;
	std::vector< TVectorT<double> > U, P;
	for (size_t v = 0; v < lines.size(); v++)
	{
		TVector3 point = lines[v].first;
		TVector3 dir = lines[v].second;
		dir -= point; m = dir.Mag();
		if (m > 0.0)
		{
			dir *= 1.0 / m;

			P.push_back(TVectorT<double>(3));
			P.back()[0] = point.X(); P.back()[1] = point.Y(); P.back()[2] = point.Z();

			U.push_back(TVectorT<double>(3));
			U.back()[0] = dir.X(); U.back()[1] = dir.Y(); U.back()[2] = dir.Z();
		}
		else mf::LogWarning("pma::SolveLeastSquares3D") << "Line undefined.";
	}
	if (P.size() < 2)
	{
		mf::LogError("pma::SolveLeastSquares3D") << "Need min. two lines.";
		return -1.0;
	}

	TVectorT<double> x(3), y(3), w(3);
	TMatrixT<double> A(3, 3);
	double ur, uc, pc;
	double s_uc2[3], s_ur_uc[3];
	double s_p_uc2[3], s_p_ur_uc[3];

	w[0] = 1.0; w[1] = 1.0; w[2] = 1.0;
	for (size_t r = 0; r < 3; r++)
	{
		y[r] = 0.0;
		for (size_t c = 0; c < 3; c++)
		{
			s_uc2[c] = 0.0; s_ur_uc[c] = 0.0;
			s_p_uc2[c] = 0.0; s_p_ur_uc[c] = 0.0;

			for (size_t v = 0; v < P.size(); v++)
			{
				//w[1] = fWeights[v]; // to remember that individual coordinates can be supressed...
				//w[2] = fWeights[v];

				ur = U[v][r];
				uc = U[v][c];
				pc = P[v][c];

				s_uc2[c] += w[r] * w[c] * (1 - uc * uc);
				s_p_uc2[c] += w[r] * w[r] * pc * (1 - uc * uc);

				s_ur_uc[c] += w[r] * w[c] * ur * uc;
				s_p_ur_uc[c] += w[r] * w[r] * pc * ur * uc;
			}

			if (r == c)
			{
				y[r] += s_p_uc2[c];
				A(r, c) = s_uc2[c];
			}
			else
			{
				y[r] -= s_p_ur_uc[c];
				A(r, c) = -s_ur_uc[c];
			}
		}
	}
	x = A.InvertFast() * y;

	result.SetXYZ(x[0], x[1], x[2]);

	TVector3 pproj;
	double dx, dy, dz, mse = 0.0;
	for (size_t v = 0; v < lines.size(); v++)
	{
		pproj = pma::GetProjectionToSegment(result, lines[v].first, lines[v].second);

		dx = result.X() - pproj.X(); // dx, dy, dz and the result point can be weighted
		dy = result.Y() - pproj.Y(); // here (linearly) by each line uncertainty
		dz = result.Z() - pproj.Z();

		mse += dx * dx + dy * dy + dz * dz;
	}
	return mse / lines.size();
}

TVector2 pma::GetProjectionToPlane(const TVector3& p, unsigned int view, unsigned int tpc, unsigned int cryo)
{
	art::ServiceHandle<geo::Geometry> geom;

	return TVector2(
		geom->TPC(tpc, cryo).Plane(view).WirePitch() * geom->WireCoordinate(p.Y(), p.Z(), view, tpc, cryo),
		p.X()
	);
}

TVector2 pma::GetVectorProjectionToPlane(const TVector3& v, unsigned int view, unsigned int tpc, unsigned int cryo)
{
	TVector3 v0_3d(0., 0., 0.);
	TVector2 v0_2d = GetProjectionToPlane(v0_3d, view, tpc, cryo);
	TVector2 v1_2d = GetProjectionToPlane(v, view, tpc, cryo);

	return v1_2d - v0_2d;
}

TVector2 pma::WireDriftToCm(unsigned int wire, float drift, unsigned int view, unsigned int tpc, unsigned int cryo)
{
	art::ServiceHandle<geo::Geometry> geom;
	const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

	return TVector2(
		geom->TPC(tpc, cryo).Plane(view).WirePitch() * wire,
		detprop->ConvertTicksToX(drift, view, tpc, cryo)
	);
}

TVector2 pma::CmToWireDrift(float xw, float yd, unsigned int view, unsigned int tpc, unsigned int cryo)
{
	art::ServiceHandle<geo::Geometry> geom;
	const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

	return TVector2(
		xw / geom->TPC(tpc, cryo).Plane(view).WirePitch(),
		detprop->ConvertXToTicks(yd, view, tpc, cryo)
	);
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

