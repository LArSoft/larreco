/**
 *  @file   PmaSegment3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track segment. See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaSegment3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::Segment3D::Segment3D(pma::Track3D* trk, pma::Node3D* vstart, pma::Node3D* vstop) :
	SortedObjectBase(vstart, vstop),
	fParent(trk)
{
	if (vstart->TPC() == vstop->TPC()) fTPC = vstart->TPC();
	if (vstart->Cryo() == vstop->Cryo()) fCryo = vstart->Cryo();
}

double pma::Segment3D::GetDistance2To(const TVector3& p3d) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	return GetDist2(p3d, v0->Point3D(), v1->Point3D());
}

double pma::Segment3D::GetDistance2To(const TVector2& p2d, unsigned int view) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	return GetDist2(p2d, v0->Projection2D(view), v1->Projection2D(view));
}

TVector3 pma::Segment3D::GetProjection(const TVector2& p, unsigned int view) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	TVector2 v0(p);
	v0 -= vStart->Projection2D(view);

	TVector2 v1(vStop->Projection2D(view));
	v1 -= vStart->Projection2D(view);
	
	TVector3 v3d(vStop->Point3D());
	v3d -= vStart->Point3D();
	
	TVector3 v3dStart(vStart->Point3D());
	TVector3 v3dStop(vStop->Point3D());

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();

	TVector3 result(0, 0, 0);
	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0 * v1 / mag;
		double b = v0Norm * cosine / v1Norm;

		if (b < 1.0)
		{
			result = v3dStart;
			if (b > 0.0) result += (v3d * b);
		}
		else result = v3dStop;
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection.";

		result = v3dStart;
		result += v3dStop;
		result *= 0.5;
	}
	return result;
}

TVector3 pma::Segment3D::GetUnconstrainedProj3D(const TVector2& p2d, unsigned int view) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	TVector2 v0(p2d);
	v0 -= vStart->Projection2D(view);

	TVector2 v1(vStop->Projection2D(view));
	v1 -= vStart->Projection2D(view);

	TVector3 v3d(vStop->Point3D());
	v3d -= vStart->Point3D();

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0 * v1 / mag;
		double b = v0Norm * cosine / v1Norm;

		return vStart->Point3D() + (v3d * b);
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection." << std::endl;

		v3d = vStart->Point3D();
		v3d += vStop->Point3D();
		v3d *= 0.5;

		return v3d;
	}
}

void pma::Segment3D::SetProjection(pma::Hit3D& h) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	TVector2 v0(h.Point2D());
	v0 -= vStart->Projection2D(h.View2D());

	TVector2 v1(vStop->Projection2D(h.View2D()));
	v1 -= vStart->Projection2D(h.View2D());

	TVector3 v3d(vStop->Point3D());
	v3d -= vStart->Point3D();

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();

	TVector2 p(vStart->Projection2D(h.View2D()));

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0 * v1 / mag;
		double b = v0Norm * cosine / v1Norm;

		p += (v1 * b);

		h.SetProjection(p, (float)b);
		h.SetPoint3D(vStart->Point3D() + (v3d * b));
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection." << std::endl;

		p += vStop->Projection2D(h.View2D());
		p *= 0.5F; h.SetProjection(p, 0.0F);

		v3d = vStart->Point3D();
		v3d += vStop->Point3D();
		v3d *= 0.5;
		h.SetPoint3D(v3d);
	}
}

double pma::Segment3D::Length2(void) const
{
	if (prev && next)
		return pma::Dist2( ((pma::Node3D*)prev)->Point3D(), ((pma::Node3D*)next)->Point3D() );
	else
	{
		mf::LogError("pma::Segment3D") << "Segment endpoints not set.";
		return 0.0;
	}
}


double pma::Segment3D::GetDist2(const TVector3& psrc, const TVector3& p0, const TVector3& p1)
{
	TVector3 v0(psrc); v0 -= p0;
	TVector3 v1(p1);   v1 -= p0;

	TVector3 v2(psrc); v2 -= p1;
	TVector3 v3(v1);   v3 *= -1.0;

	double v0Norm2 = v0.Mag2();
	double v1Norm2 = v1.Mag2();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm2 < eps)
	{
		mf::LogWarning("pma::Segment3D") << "Short segment or its projection.";
		v1 = p0; v1 += p1; v1 *= 0.5;
		return pma::Dist2(v1, psrc);
	}

	double mag01 = sqrt(v0Norm2 * v1Norm2);
	double cosine01 = 0.0;
	if (mag01 != 0.0) cosine01 = v0 * v1 / mag01;

	double v2Norm2 = v2.Mag2();
	double mag23 = sqrt(v2Norm2 * v3.Mag2());
	double cosine23 = 0.0;
	if (mag23 != 0.0) cosine23 = v2 * v3 / mag23;

	double result = 0.0;
	if ((cosine01 > 0.0) && (cosine23 > 0.0))
	{
		result = (1.0 - cosine01 * cosine01) * v0Norm2;
	}
	else // increase distance to prefer hit assigned to the vertex, not segment
	{
		if (cosine01 <= 0.0) result = 1.0001 * v0Norm2;
		else result = 1.0001 * v2Norm2;
	}

	if (result >= 0.0) return result;
	else return 0.0;
}

double pma::Segment3D::GetDist2(const TVector2& psrc, const TVector2& p0, const TVector2& p1)
{
	TVector2 v0(psrc); v0 -= p0;
	TVector2 v1(p1);   v1 -= p0;

	TVector2 v2(psrc); v2 -= p1;
	TVector2 v3(v1);   v3 *= -1.0;

	double v0Norm2 = v0.Mod2();
	double v1Norm2 = v1.Mod2();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm2 < eps)
	{
		mf::LogVerbatim("pma::Segment3D") << "Short segment or its projection.";
		v1 = p0; v1 += p1; v1 *= 0.5;
		return pma::Dist2(v1, psrc);
	}

	double mag01 = sqrt(v0Norm2 * v1Norm2);
	double cosine01 = 0.0;
	if (mag01 != 0.0) cosine01 = v0 * v1 / mag01;

	double v2Norm2 = v2.Mod2();
	double mag23 = sqrt(v2Norm2 * v3.Mod2());
	double cosine23 = 0.0;
	if (mag23 != 0.0) cosine23 = v2 * v3 / mag23;

	double result = 0.0;
	if ((cosine01 > 0.0) && (cosine23 > 0.0))
	{
		result = (1.0 - cosine01 * cosine01) * v0Norm2;
	}
	else // increase distance to prefer hit assigned to the vertex, not segment
	{
		if (cosine01 <= 0.0) result = 1.0001 * v0Norm2;
		else result = 1.0001 * v2Norm2;
	}

	if (result >= 0.0) return result;
	else return 0.0;
}

