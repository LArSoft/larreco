/**
 *  @file   PmaNode3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track node. See PmaTrack3D.h file for details.
 */

#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// Fixed optimization directions:     X      Y      Z
bool pma::Node3D::fGradFixed[3] = { false, false, false };

double pma::Node3D::fMargin = 3.0;

pma::Node3D::Node3D(void) :
	fMinX(0), fMaxX(0),
	fMinY(0), fMaxY(0),
	fMinZ(0), fMaxZ(0),
	fPoint3D(0, 0, 0),
	fDriftOffset(0),
	fIsVertex(false)
{
	fTPC = 0; fCryo = 0;

	fProj2D[0].Set(0);
	fProj2D[1].Set(0);
	fProj2D[2].Set(0);
}

pma::Node3D::Node3D(const TVector3& p3d, unsigned int tpc, unsigned int cryo, bool vtx, double xshift) :
    fDriftOffset(xshift),
	fIsVertex(vtx)
{
	fTPC = tpc; fCryo = cryo;

	const auto& tpcGeo = fGeom->TPC(tpc, cryo);

	unsigned int lastPlane = geo::kZ;
	while ((lastPlane > 0) && !tpcGeo.HasPlane(lastPlane)) lastPlane--;

	auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	fMinX = detprop->ConvertTicksToX(0, lastPlane, tpc, cryo);
	fMaxX = detprop->ConvertTicksToX(detprop->NumberTimeSamples() - 1, lastPlane, tpc, cryo);
	if (fMaxX < fMinX) { double tmp = fMaxX; fMaxX = fMinX; fMinX = tmp; }

	fMinY = tpcGeo.MinY(); fMaxY = tpcGeo.MaxY();
	fMinZ = tpcGeo.MinZ(); fMaxZ = tpcGeo.MaxZ();

	for (size_t i = 0; i < 3; i++)
	{
		if (tpcGeo.HasPlane(i)) fWirePitch[i] = tpcGeo.Plane(i).WirePitch();
		else fWirePitch[i] = 0.0;
	}

	SetPoint3D(p3d);
}

double pma::Node3D::GetDistToWall(void) const
{
	double d, dmin = fPoint3D.X() - fMinX;
	d = fMaxX - fPoint3D.X();
	if (d < dmin) dmin = d;

	d = fPoint3D.Y() - fMinY;
	if (d < dmin) dmin = d;
	d = fMaxY - fPoint3D.Y();
	if (d < dmin) dmin = d;

	d = fPoint3D.Z() - fMinZ;
	if (d < dmin) dmin = d;
	d = fMaxZ - fPoint3D.Z();
	if (d < dmin) dmin = d;

	return dmin;
}

bool pma::Node3D::SameTPC(const TVector3& p3d, float margin) const
{
	if (((fMinX - margin) <= p3d.X()) && (p3d.X() <= (fMaxX + margin)) &&
	    ((fMinY - margin) <= p3d.Y()) && (p3d.Y() <= (fMaxY + margin)) &&
	    ((fMinZ - margin) <= p3d.Z()) && (p3d.Z() <= (fMaxZ + margin))) return true;
	else return false;
}

bool pma::Node3D::LimitPoint3D(void)
{
	bool trimmed = false;

	if (fPoint3D.X() < fMinX - fMargin) { fPoint3D.SetX(fMinX - fMargin); trimmed = true; }
	if (fPoint3D.X() > fMaxX + fMargin) { fPoint3D.SetX(fMaxX + fMargin); trimmed = true; }

	if (fPoint3D.Y() < fMinY - fMargin) { fPoint3D.SetY(fMinY - fMargin); trimmed = true; }
	if (fPoint3D.Y() > fMaxY + fMargin) { fPoint3D.SetY(fMaxY + fMargin); trimmed = true; }

	if (fPoint3D.Z() < fMinZ - fMargin) { fPoint3D.SetZ(fMinZ - fMargin); trimmed = true; }
	if (fPoint3D.Z() > fMaxZ + fMargin) { fPoint3D.SetZ(fMaxZ + fMargin); trimmed = true; }

	return trimmed;
}

void pma::Node3D::UpdateProj2D(void)
{
	fProj2D[0].Set(
		fWirePitch[0] * fGeom->WireCoordinate(fPoint3D.Y(), fPoint3D.Z(), geo::kU, fTPC, fCryo),
		fPoint3D.X() - fDriftOffset
	);

	fProj2D[1].Set(
		fWirePitch[1] * fGeom->WireCoordinate(fPoint3D.Y(), fPoint3D.Z(), geo::kV, fTPC, fCryo),
		fPoint3D.X() - fDriftOffset
	);
	
	fProj2D[2].Set(
		fWirePitch[2] * fGeom->WireCoordinate(fPoint3D.Y(), fPoint3D.Z(), geo::kZ, fTPC, fCryo),
		fPoint3D.X() - fDriftOffset
	);
}

bool pma::Node3D::SetPoint3D(const TVector3& p3d)
{
	fPoint3D = p3d;

	bool accepted = !LimitPoint3D();
	UpdateProj2D();

	return accepted;
}

double pma::Node3D::GetDistance2To(const TVector3& p3d) const
{
	return pma::Dist2(fPoint3D, p3d);
}

double pma::Node3D::GetDistance2To(const TVector2& p2d, unsigned int view) const
{
	return pma::Dist2(fProj2D[view], p2d);
}

pma::Vector3D pma::Node3D::GetDirection3D(void) const
{
    pma::Element3D* seg = 0;
    if (next) { seg = dynamic_cast< pma::Element3D* >(next); }
    else if (prev) { seg = dynamic_cast< pma::Element3D* >(prev); }
    else { throw cet::exception("Node3D") << "Isolated vertex." << std::endl; }

    if (seg) { return seg->GetDirection3D(); }
    else { throw cet::exception("Node3D") << "Corrupted vertex." << std::endl; }
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
		mf::LogError("pma::Node3D") << "Isolated vertex.";
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
		if (fFrozen) // limit 3D positions to outermose node if frozen
		{
			h.SetPoint3D(fPoint3D);
		}
		else // or set 3D positions along the line of outermost segment
		{
			g3d -= fPoint3D;
			h.SetPoint3D(fPoint3D + (g3d * b));

			p += (v1 * b);
		}
		h.SetProjection(p, -b);
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
		mf::LogError("pma::Node3D") << "pma::Node3D::SegmentCos(): neighbours not initialized.";
		return -1.0;
	}
}

double pma::Node3D::SegmentCosWirePlane(void) const
{
	if (prev && next)
	{
		pma::Node3D* vStop1 = static_cast< pma::Node3D* >(prev->Prev());
		pma::Node3D* vStop2 = static_cast< pma::Node3D* >(next->Next());
		TVector2 v1(vStop1->fPoint3D.Y() - fPoint3D.Y(), vStop1->fPoint3D.Z() - fPoint3D.Z());
		TVector2 v2(vStop2->fPoint3D.Y() - fPoint3D.Y(), vStop2->fPoint3D.Z() - fPoint3D.Z());
		double mag = sqrt(v1.Mod2() * v2.Mod2());
		double cosine = 0.0;
		if (mag != 0.0) cosine = v1 * v2 / mag;
		return cosine;
	}
	else
	{
		mf::LogError("pma::Node3D") << "pma::Node3D::SegmentCosZX(): neighbours not initialized.";
		return -1.0;
	}
}

double pma::Node3D::SegmentCosTransverse(void) const
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
		mf::LogError("pma::Node3D") << "pma::Node3D::SegmentCosZY(): neighbours not initialized.";
		return -1.0;
	}
}

// *** Note: should be changed / generalized for horizontal wire planes (e.g. 2-phase LAr). ***
double pma::Node3D::EndPtCos2Transverse(void) const
{
	if (prev && next)
	{
		pma::Node3D* vStart = static_cast< pma::Node3D* >(prev->Prev());
		pma::Node3D* vStop = static_cast< pma::Node3D* >(next->Next());

		double dy = vStop->Point3D().X() - vStart->Point3D().X();
		double dz = vStop->Point3D().Z() - vStart->Point3D().Z();
		double len2 = dy * dy + dz * dz;
		double cosine2 = 0.0;
		if (len2 > 0.0) cosine2 = dz * dz / len2;
		return cosine2;
	}
	else return 0.0;
}

double pma::Node3D::PiInWirePlane(void) const
{
	if (prev && NextCount())
	{
		pma::Segment3D* seg0 = dynamic_cast< pma::Segment3D* >(prev);
		pma::Segment3D* seg1 = dynamic_cast< pma::Segment3D* >(Next(0));
		unsigned int nInd1 = NHits(geo::kU) + seg0->NHits(geo::kU) + seg1->NHits(geo::kU);

		if (fHitsRadius > 0.0F)
			return (1.0 + SegmentCosWirePlane()) * fHitsRadius * fHitsRadius / (4 * nInd1 + 1.0);
		else return (1.0 + SegmentCosWirePlane()) * Length2() / (4 * nInd1 + 1.0);
	}
	else return 0.0;
}

// Constraint on two segments angle in projection to plane parallel to wire plane, suppressed by
// the orientation in the plane transverse to wire plane (only sections with low variation of
// drift time are penalized with this constraint); PiInWirePlane() components are reduced if
// there are Ind1 / geo::kU hits which add information to the object shape.
// *** Note: should be changed / generalized for horizontal wire planes (e.g. 2-phase LAr). ***
double pma::Node3D::PenaltyInWirePlane(void) const
{
	if (fIsVertex) return 0.0;

	unsigned int nseg = 1;
	double penalty = PiInWirePlane();
	pma::Node3D* v;
	if (next)
	{
		v = static_cast< pma::Node3D* >(next->Next());
		penalty += v->PiInWirePlane(); nseg++;
	}
	if (prev)
	{
		v = static_cast< pma::Node3D* >(prev->Prev());
		penalty += v->PiInWirePlane(); nseg++;
	}
	if (penalty > 0.0) return pow(EndPtCos2Transverse(), 10) * penalty / nseg;
	else return 0.0;
}

bool pma::Node3D::IsBranching(void) const
{
	size_t nnext = NextCount();
	if (nnext > 1) return true; // 1 trk -> vtx -> n*trk

	if (prev && nnext)
	{
		pma::Segment3D* segPrev = static_cast< pma::Segment3D* >(prev);
		pma::Segment3D* segNext = static_cast< pma::Segment3D* >(Next(0));
		if (segNext->Parent() != segPrev->Parent()) // 1 trk -> vtx -> 1 trk
			return true;
	}
	return false;
}

bool pma::Node3D::IsTPCEdge(void) const
{
	if (prev && (NextCount() == 1))
	{
		pma::Segment3D* segPrev = static_cast< pma::Segment3D* >(prev);
		pma::Segment3D* segNext = static_cast< pma::Segment3D* >(Next(0));

		if ((segPrev->TPC() < 0) || (segNext->TPC() < 0)) return true;
	}
	return false;
}

std::vector< pma::Track3D* > pma::Node3D::GetBranches(void) const
{
	std::vector< pma::Track3D* > branches;
	if (NextCount())
	{
		branches.reserve(NextCount());
		for (size_t i = 0; i < NextCount(); ++i)
		{
			pma::Segment3D* seg = static_cast< pma::Segment3D* >(Next(i));
			branches.push_back(seg->Parent());
		}
	}
	return branches;
}

double pma::Node3D::Pi(float endSegWeight, bool doAsymm) const
{
	if (fIsVertex) return 0.0;

	if (prev && NextCount())
	{
		pma::Segment3D* segPrev = static_cast< pma::Segment3D* >(prev);
		pma::Segment3D* segNext = static_cast< pma::Segment3D* >(Next(0));

		double scale = 1.0;
		if ((segPrev->TPC() < 0) || (segNext->TPC() < 0)) scale = 0.5; // lower penalty on segments between tpc's

		double segCos = SegmentCos();

		double lAsymmFactor = 0.0;
		if (doAsymm)
		{
			double lPrev = segPrev->Length();
			double lNext = segNext->Length();
			double lSum = lPrev + lNext;
			if (lSum > 0.1)
			{
				double lAsymm = (1.0 - segCos) * (lPrev - lNext) / lSum;
				lAsymmFactor = 0.05 * lAsymm * lAsymm;
			}
		}

		if (fHitsRadius > 0.0F) return scale * (1.0 + segCos + lAsymmFactor) * fHitsRadius * fHitsRadius;
		else return scale * (1.0 + segCos + lAsymmFactor) * Length2();
	}
	else
	{
		double pi_result = 0.0;
		unsigned int nSeg = 0;
		pma::Segment3D* seg = 0;
		if (prev)
		{
			seg = static_cast< pma::Segment3D* >(prev);

			SortedObjectBase* prevVtx = seg->Prev();
			if (prevVtx->Prev()) nSeg++;
			nSeg += prevVtx->NextCount();
		}
		else if (next)
		{
			seg = static_cast< pma::Segment3D* >(next);
			
			SortedObjectBase* nextVtx = seg->Next(0);
			nSeg += nextVtx->NextCount() + 1;
		}
		else
		{
			mf::LogWarning("pma::Node3D") << "pma::Node3D::Pi(): an isolated vertex?";
			return 0.0;
		}
		if (nSeg == 1) pi_result = endSegWeight * seg->Length2();
		return pi_result;
	}
}

double pma::Node3D::Penalty(float endSegWeight) const
{
	unsigned int nseg = 1;
	double penalty = Pi(endSegWeight, true);

	pma::Node3D* v;
	for (unsigned int i = 0; i < NextCount(); i++)
	{
		v = static_cast< pma::Node3D* >(Next(i)->Next());
		penalty += v->Pi(endSegWeight, false); nseg++;
	}
	if (prev)
	{
		v = static_cast< pma::Node3D* >(prev->Prev());
		penalty += v->Pi(endSegWeight, false); nseg++;
	}
	return penalty / nseg;
}

double pma::Node3D::Mse(void) const
{
	unsigned int nhits = NPrecalcEnabledHits(); //NEnabledHits();
	double mse = SumDist2();

	pma::Segment3D* seg;
	for (unsigned int i = 0; i < NextCount(); i++)
	{
		seg = static_cast< pma::Segment3D* >(Next(i));
		nhits += seg->NPrecalcEnabledHits(); //NEnabledHits();
		mse += seg->SumDist2();
	}
	if (prev)
	{
		seg = static_cast< pma::Segment3D* >(prev);
		nhits += seg->NPrecalcEnabledHits(); //NEnabledHits();
		mse += seg->SumDist2();
	}
	if (!nhits) return 0.0;
	else return mse / nhits;
}

double pma::Node3D::GetObjFunction(float penaltyValue, float endSegWeight) const
{
	return Mse() + penaltyValue * (Penalty(endSegWeight) + PenaltyInWirePlane());
}

double pma::Node3D::MakeGradient(float penaltyValue, float endSegWeight)
{	
	double l1 = 0.0, l2 = 0.0, minLength2 = 0.0;
	TVector3 tmp(fPoint3D), gpoint(fPoint3D);

	pma::Segment3D* seg;
	if (prev)
	{
		seg = static_cast< pma::Segment3D* >(prev);
		l1 = seg->Length2();
	}
	if (next)
	{
		seg = static_cast< pma::Segment3D* >(next);
		l2 = seg->Length2();
	}
	if ((l1 > 0.01) && (l1 < l2)) minLength2 = l1;
	else if ((l2 > 0.01) && (l2 < l1)) minLength2 = l2;
	else minLength2 = 0.01;

	double dxi = 0.001 * sqrt(minLength2);

	if (dxi < 6.0E-37) return 0.0;

	double gi, g0, gz;
	gz = g0 = GetObjFunction(penaltyValue, endSegWeight);

	//if (fQPenaltyFactor > 0.0F) gz += fQPenaltyFactor * QPenalty(); <----------------------- maybe later..

	if (!fGradFixed[0]) // gradX
	{
		gpoint[0] = tmp[0] + dxi;
		SetPoint3D(gpoint);
		gi = GetObjFunction(penaltyValue, endSegWeight);
		fGradient[0] = (g0 - gi) / dxi;

		gpoint[0] = tmp[0] - dxi;
		SetPoint3D(gpoint);
		gi = GetObjFunction(penaltyValue, endSegWeight);
		fGradient[0] = 0.5 * (fGradient[0] + (gi - g0) / dxi);

		gpoint[0] = tmp[0];
	}

	if (!fGradFixed[1]) // gradY
	{
		gpoint[1] = tmp[1] + dxi;
		SetPoint3D(gpoint);
		gi = GetObjFunction(penaltyValue, endSegWeight);
		fGradient[1] = (g0 - gi) / dxi;

		gpoint[1] = tmp[1] - dxi;
		SetPoint3D(gpoint);
		gi = GetObjFunction(penaltyValue, endSegWeight);
		fGradient[1] = 0.5 * (fGradient[1] + (gi - g0) / dxi);

		gpoint[1] = tmp[1];
	}

	if (!fGradFixed[2]) // gradZ
	{
		gpoint[2] = tmp[2] + dxi;
		SetPoint3D(gpoint);
		gi = GetObjFunction(penaltyValue, endSegWeight);
		//if (fQPenaltyFactor > 0.0F) gi += fQPenaltyFactor * QPenalty();
		fGradient[2] = (gz - gi) / dxi;

		gpoint[2] = tmp[2] - dxi;
		SetPoint3D(gpoint);
		gi = GetObjFunction(penaltyValue, endSegWeight);
		//if (fQPenaltyFactor > 0.0F) gi += fQPenaltyFactor * QPenalty();
		fGradient[2] = 0.5 * (fGradient[2] + (gi - gz) / dxi);

		gpoint[2] = tmp[2];
	}

	SetPoint3D(tmp);
	if (fGradient.Mag2() < 6.0E-37) return 0.0;

	return g0;
}

double pma::Node3D::StepWithGradient(float alfa, float tol, float penalty, float weight)
{
	unsigned int steps = 0;
	double t, t1, t2, t3, g, g0, g1, g2, g3, p1, p2;
	double eps = 6.0E-37, zero_tol = 1.0E-15;
	TVector3 tmp(fPoint3D), gpoint(fPoint3D);

	g = MakeGradient(penalty, weight);
	if (g < zero_tol) return 0.0;
	g0 = g;

	//**** first three points ****//
	alfa *= 0.8F;
	t2 = 0.0; g2 = g;
	t3 = 0.0; g3 = g;
	do
	{
		t1 = t2; g1 = g2;
		t2 = t3; g2 = g3;

		alfa *= 1.25F;
		t3 += alfa;
		gpoint = tmp;
		gpoint += (fGradient * t3);
		if (!SetPoint3D(gpoint)) // stepped out of allowed volume
		{
			//std::cout << "****  SetPoint trimmed 1 ****" << std::endl;
			g3 = GetObjFunction(penalty, weight);
			if (g3 < g2) return (g0 - g3) / g3;   // exit with the node at the border
			else { SetPoint3D(tmp); return 0.0; } // exit with restored original position
		}
		
		g3 = GetObjFunction(penalty, weight);
		
		if (g3 < zero_tol) return 0.0;

		if (++steps > 1000) { SetPoint3D(tmp); return 0.0; }

	} while (g3 < g2);
	//****************************//

	//*** first step overshoot ***//
	if (steps == 1)
	{
		t2 = t3; g2 = g3;
		do
		{
			t3 = t2; g3 = g2;
			t2 = (t1 * g3 + t3 * g1) / (g1 + g3);

			// small shift...
			t2 = 0.05 * t3 + 0.95 * t2;

			// break: starting point is at the minimum
			//if (t2 == t1) { SetPoint3D(tmp); return 0.0F; }

			// break: starting point is very close to the minimum
			if (fabs(t2 - t1) < tol) { SetPoint3D(tmp); return 0.0; }

			gpoint = tmp;
			gpoint += (fGradient * t2);
			if (!SetPoint3D(gpoint)) // select the best point to exit
			{
				//std::cout << "****  SetPoint trimmed 2 ****" << std::endl;
				g2 = GetObjFunction(penalty, weight);
				if (g2 < g0) return (g0 - g2) / g2;   // exit with the node at the border
				else if (g1 < g0)
				{
					gpoint = tmp; gpoint += (fGradient * t1);
					return (g0 - g1) / g1;
				}
				else if (g3 < g0)
				{
					gpoint = tmp; gpoint += (fGradient * t3);
					return (g0 - g3) / g3;
				}
				else { SetPoint3D(tmp); return 0.0; }
			}
			g2 = GetObjFunction(penalty, weight);

			if (g2 < zero_tol) return 0.0;
			steps++;

		} while (g2 >= g1);
	}
	//****************************//

	while (fabs(t1 - t3) > tol)
	{
		//*** 3-point minimization ***//
		if ((fabs(t2 - t1) < eps) || (fabs(t2 - t3) < eps))
			break; // minimum on the edge
		if ((fabs(g2 - g1) < eps) && (fabs(g2 - g3) < eps))
			break; // ~singularity

		p1 = (t2 - t1) * (g2 - g3);
		p2 = (t2 - t3) * (g2 - g1);
		if (fabs(p1 - p2) < eps) break; // ~linearity

		t = t2 + ((t2 - t1) * p1 - (t2 - t3) * p2) / (2 * (p2 - p1));
		if ((t <= t1) || (t >= t3))
			t = (t1 * g3 + t3 * g1) / (g1 + g3);

		gpoint = tmp;
		gpoint += (fGradient * t);
		if (!SetPoint3D(gpoint)) // select the best point to exit
		{
			//std::cout << "****  SetPoint trimmed 3 ****" << std::endl;
			g = GetObjFunction(penalty, weight);
			if ((g < g0) && (g < g1) && (g < g3)) return (g0 - g) / g;   // exit with the node at the border
			else if ((g1 < g0) && (g1 < g3))
			{
				gpoint = tmp; gpoint += (fGradient * t1);
				return (g0 - g1) / g1;
			}
			else if (g3 < g0)
			{
				gpoint = tmp; gpoint += (fGradient * t3);
				return (g0 - g3) / g3;
			}
			else { SetPoint3D(tmp); return 0.0; }
		}
		
		g = GetObjFunction(penalty, weight);
		if (g < zero_tol) return 0.0;
		steps++;
		//****************************//

		//*** select next 3 points ***//
		if (fabs(t - t2) < 0.2 * tol) break; // start in a new direction
		if (g < g2)
		{
			if (t < t2) { t3 = t2; g3 = g2; }
			else { t1 = t2; g1 = g2; }
			t2 = t; g2 = g;
		}
		else
		{
			if (t < t2) { t1 = t; g1 = g; }
			else { t3 = t; g3 = g; }
		}
		//****************************//
	}

	return (g0 - g) / g;
}

void pma::Node3D::Optimize(float penaltyValue, float endSegWeight)
{
	if (!fFrozen)
	{
		double dg = StepWithGradient(0.1F, 0.002F, penaltyValue, endSegWeight);
		if (dg > 0.01) dg = StepWithGradient(0.03F, 0.0001F, penaltyValue, endSegWeight);
		if (dg > 0.0) dg = StepWithGradient(0.03F, 0.0001F, penaltyValue, endSegWeight);
	}
}

void pma::Node3D::ClearAssigned(pma::Track3D* trk)
{
	if (!trk)
	{
		// like in the base class:
		fAssignedPoints.clear();
		fAssignedHits.clear();
	}
	else
	{
		std::vector< pma::Track3D* > to_check;
		pma::Segment3D* seg;
		if (Prev())
		{
			seg = static_cast< pma::Segment3D* >(Prev());
			if (seg->Parent() != trk) to_check.push_back(seg->Parent());
		}
		for (unsigned int i = 0; i < NextCount(); i++)
		{
			seg = static_cast< pma::Segment3D* >(Next(i));
			if (seg->Parent() != trk) to_check.push_back(seg->Parent());
		}
			
		unsigned int p = 0;
		while (p < fAssignedPoints.size())
		{
			bool found = false;
			for (size_t t = 0; t < to_check.size(); t++)
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
			pma::Hit3D* hit = fAssignedHits[h];

			for (size_t t = 0; (t < to_check.size()) && !found; t++)
				for (size_t i = 0; i < to_check[t]->size(); i++)
				{
					pma::Hit3D* pmaHit = static_cast< pma::Hit3D* >((*(to_check[t]))[i]);
					if (hit == pmaHit)
					{
						found = true; break;
					}
				}

			if (!found) fAssignedHits.erase(fAssignedHits.begin() + h);
			else h++;
		}
	}

	fHitsRadius = 0.0F;
}

