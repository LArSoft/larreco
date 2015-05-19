/**
 *  @file   PmaNode3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track node.
 */

#include "RecoAlg/PMAlg/PmaNode3D.h"
#include "RecoAlg/PMAlg/PmaSegment3D.h"

#include "Geometry/TPCGeo.h"

pma::Node3D::Node3D(void) :
	fTPC(0), fCryo(0),
	fMinX(0), fMaxX(0),
	fMinY(0), fMaxY(0),
	fMinZ(0), fMaxZ(0),
	fPoint3D(0, 0, 0)
{
	fProj2D[0].Set(0);
	fProj2D[1].Set(0);
	fProj2D[2].Set(0);
}

pma::Node3D::Node3D(const TVector3& p3d, unsigned int tpc, unsigned int cryo) :
	fTPC(tpc), fCryo(cryo)
{
	const auto& tpcGeo = fGeom->TPC(tpc, cryo);
	fMinX = tpcGeo.MinX(); fMaxX = tpcGeo.MaxX();
	fMinY = tpcGeo.MinY(); fMaxY = tpcGeo.MaxY();
	fMinZ = tpcGeo.MinZ(); fMaxZ = tpcGeo.MaxZ();
	SetPoint3D(p3d);
}

void pma::Node3D::LimitPoint3D(float margin)
{
	if (fPoint3D.X() < fMinX - margin) fPoint3D.SetX(fMinX - margin);
	if (fPoint3D.X() > fMaxX + margin) fPoint3D.SetX(fMaxX + margin);

	if (fPoint3D.Y() < fMinY - margin) fPoint3D.SetY(fMinY - margin);
	if (fPoint3D.Y() > fMaxY + margin) fPoint3D.SetY(fMaxY + margin);

	if (fPoint3D.Z() < fMinZ - margin) fPoint3D.SetZ(fMinZ - margin);
	if (fPoint3D.Z() > fMaxZ + margin) fPoint3D.SetZ(fMaxZ + margin);
}

void pma::Node3D::UpdateProj2D(void)
{
	fProj2D[0].Set(
		fGeom->WireCoordinate(fPoint3D.Y(), fPoint3D.Z(), geo::kU, fTPC, fCryo),
		fDetProp->ConvertXToTicks(fPoint3D.X(), geo::kU, fTPC, fCryo)
	);

	fProj2D[1].Set(
		fGeom->WireCoordinate(fPoint3D.Y(), fPoint3D.Z(), geo::kV, fTPC, fCryo),
		fDetProp->ConvertXToTicks(fPoint3D.X(), geo::kV, fTPC, fCryo)
	);

	fProj2D[2].Set(
		fGeom->WireCoordinate(fPoint3D.Y(), fPoint3D.Z(), geo::kZ, fTPC, fCryo),
		fDetProp->ConvertXToTicks(fPoint3D.X(), geo::kZ, fTPC, fCryo)
	);
}

void pma::Node3D::SetPoint3D(const TVector3& p3d)
{
	fPoint3D = p3d;
	LimitPoint3D();
	UpdateProj2D();
}

double pma::Node3D::GetDistance2To(const TVector3& p3d) const
{
	return pma::Dist2(fPoint3D, p3d);
}

double pma::Node3D::GetDistance2To(const TVector2& p2d, unsigned int view) const
{
	return pma::Dist2(fProj2D[view], p2d);
}

void pma::Node3D::SetProjection(pma::Hit3D& h) const
{
	TVector2 gstart;
	TVector3 g3d;
	if (prev)
	{
		pma::Node3D* vtx = static_cast< pma::Node3D* >(prev->Prev());
		gstart = vtx->Projection2D(h.View2D());
		if (!next) g3d = vtx->Point3D();
	}
	else if (next)
	{
		pma::Node3D* vtx = static_cast< pma::Node3D* >(next->Next());
		gstart = Projection2D(h.View2D());
		gstart -= vtx->Projection2D(h.View2D()) - Projection2D(h.View2D());
		if (!prev)
		{
			g3d = fPoint3D;
			g3d -= vtx->Point3D() - fPoint3D;
		}
	}
	else
	{
		std::cout << "Isolated vertex." << std::endl;
		TVector2 p(Projection2D(h.View2D()));
		h.SetProjection(p, 0.0F);
		h.SetPoint3D(fPoint3D);
		return;
	}

	TVector2 v0(h.Point2D());
	v0 -= Projection2D(h.View2D());

	TVector2 v1(gstart);
	v1 -= Projection2D(h.View2D());

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();
	double mag = v0Norm * v1Norm;
	double cosine = 0.0;
	if (mag != 0.0) cosine = v0 * v1 / mag;

	TVector2 p(Projection2D(h.View2D()));

	if (prev && next)
	{
		pma::Node3D* vNext = static_cast< pma::Node3D* >(next->Next());
		TVector2 vN(vNext->Projection2D(h.View2D()));
		vN -= Projection2D(h.View2D());

		mag = v0Norm * vN.Mod();
		double cosineN = 0.0;
		if (mag != 0.0) cosineN = v0 * vN / mag;

		// hit on the previous segment side, sorting on the -cosine(prev_seg, point)  /max.val. = 1/
		if (cosineN <= cosine) h.SetProjection(p, -(float)cosine);
		// hit on the next segment side, sorting on the 1+cosine(next_seg, point)  /min.val. = 1/
		else h.SetProjection(p, 2.0F + (float)cosineN);

		h.SetPoint3D(fPoint3D);
	}
	else
	{
		float b = (float)(v0Norm * cosine / v1Norm);
		p += (v1 * b);
		h.SetProjection(p, -b);

		g3d -= fPoint3D;
		h.SetPoint3D(fPoint3D + (g3d * b));
	}
}

double pma::Node3D::Length2(void) const
{
	double l = 0.0;
	if (next) l += (static_cast< pma::Segment3D* >(next))->Length();
	if (prev) l += (static_cast< pma::Segment3D* >(prev))->Length();

	if (next && prev) return 0.25 * l * l;
	else return l * l;
}

double pma::Node3D::SegmentCos(void) const
{
	if (prev && next)
	{
		pma::Node3D* vStop1 = static_cast< pma::Node3D* >(prev->Prev());
		pma::Node3D* vStop2 = static_cast< pma::Node3D* >(next->Next());
		TVector3 v1(vStop1->fPoint3D); v1 -= fPoint3D;
		TVector3 v2(vStop2->fPoint3D); v2 -= fPoint3D;
		double mag = sqrt(v1.Mag2() * v2.Mag2());
		double cosine = 0.0;
		if (mag != 0.0) cosine = v1 * v2 / mag;
		return cosine;
	}
	else
	{
		std::cout << "PlVertex3D::SegmentCos(): neighbours not initialized." << std::endl;
		return -1.0;
	}
}

double pma::Node3D::SegmentCosWirePlane(void) const
{
	if (prev && next)
	{
		pma::Node3D* vStop1 = static_cast< pma::Node3D* >(prev->Prev());
		pma::Node3D* vStop2 = static_cast< pma::Node3D* >(next->Next());
		TVector2 v1, v2;
		v1.Set( vStop1->fPoint3D.Y() - fPoint3D.Y(),
		        vStop1->fPoint3D.Z() - fPoint3D.Z());
		v2.Set( vStop2->fPoint3D.Y() - fPoint3D.Y(),
		        vStop2->fPoint3D.Z() - fPoint3D.Z());
		double mag = sqrt(v1.Mod2() * v2.Mod2());
		double cosine = 0.0;
		if (mag != 0.0) cosine = v1 * v2 / mag;
		return cosine;
	}
	else
	{
		std::cout << "PlVertex3D::SegmentCosZX(): neighbours not initialized." << std::endl;
		return -1.0;
	}
}

double pma::Node3D::SegmentCosHorizontal(void) const
{
	if (prev && next)
	{
		pma::Node3D* vStop1 = static_cast< pma::Node3D* >(prev->Prev());
		pma::Node3D* vStop2 = static_cast< pma::Node3D* >(next->Next());
		TVector2 v1, v2;
		v1.Set( vStop1->fPoint3D.X() - fPoint3D.X(),
		        vStop1->fPoint3D.Z() - fPoint3D.Z());
		v2.Set( vStop2->fPoint3D.X() - fPoint3D.X(),
		        vStop2->fPoint3D.Z() - fPoint3D.Z());
		double mag = sqrt(v1.Mod2() * v2.Mod2());
		double cosine = 0.0;
		if (mag != 0.0) cosine = v1 * v2 / mag;
		return cosine;
	}
	else
	{
		std::cout << "PlVertex3D::SegmentCosZY(): neighbours not initialized." << std::endl;
		return -1.0;
	}
}

double pma::Node3D::GetObjFunction(float penaltyValue, float endSegWeight) const
{
	return 0.;
}

void pma::Node3D::Optimize(float penaltyValue, float endSegWeight)
{
}

void pma::Node3D::ClearAssigned(pma::Track3D* trk)
{
	if (!trk)
	{
		// like in the base class:
		fAssignedPoints.clear();
		fAssignedHits.clear();
	}
/*	else
	{
		std::vector< AF::PolygonalLine3D* > to_check;
		AF::PlSegment3D* seg;
		if (Prev())
		{
			seg = static_cast< AF::PlSegment3D* >(Prev());
			if (seg->Parent() != pl) to_check.push_back(seg->Parent());
		}
		for (unsigned int i = 0; i < NextCount(); i++)
		{
			seg = static_cast< AF::PlSegment3D* >(Next(i));
			if (seg->Parent() != pl) to_check.push_back(seg->Parent());
		}
			
		unsigned int p = 0;
		while (p < fAssignedPoints.size())
		{
			bool found = false;
			for (unsigned int t = 0; t < to_check.size(); t++)
				if (to_check[t]->HasRefPoint(fAssignedPoints[p]))
				{
					found = true; break;
				}

			if (!found) fAssignedPoints.erase(fAssignedPoints.begin() + p);
			else p++;
		}

		unsigned int h = 0;
		while (h < fAssignedHits.size())
		{
			bool found = false;
			AF::PlHit3D* hit = fAssignedHits[h];

			for (unsigned int t = 0; (t < to_check.size()) && !found; t++)
				for (unsigned int i = 0; i < to_check[t]->size(); i++)
				{
					AF::PlHit3D* plHit = static_cast< AF::PlHit3D* >((*(to_check[t]))[i]);
					if (hit == plHit)
					{
						found = true; break;
					}
				}

			if (!found) fAssignedHits.erase(fAssignedHits.begin() + h);
			else h++;
		}
	}
*/
	fHitsRadius = 0.0F;
}

