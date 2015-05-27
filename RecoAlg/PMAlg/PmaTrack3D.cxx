/**
 *  @file   PmaTrack3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaTrack3D.h"

pma::Track3D::Track3D(void) :
	fMaxHitsPerSeg(70),
	fPenaltyFactor(2.0F),
	fMaxSegStopFactor(3.0F),

	fSegStopValue(2), fMinSegStop(2), fMaxSegStop(2),
	fSegStopFactor(0.2F),
	fPenaltyValue(0.1F),
	fEndSegWeight(0.05F),
	fHitsRadius(1.0F)
{
}

pma::Track3D::~Track3D(void)
{
	for (size_t i = 0; i < fHits.size(); i++) delete fHits[i];
	for (size_t i = 0; i < fAssignedPoints.size(); i++) delete fAssignedPoints[i];

	for (size_t i = 0; i < fSegments.size(); i++) delete fSegments[i];
	for (size_t i = 0; i < fNodes.size(); i++)
		if (!fNodes[i]->NextCount() && !fNodes[i]->Prev()) delete fNodes[i];
}

void pma::Track3D::Initialize(float initEndSegW)
{
	if (InitFromRefPoints()) std::cout << "pma::Track3D initialized with 3D reference points." << std::endl;
	else
	{
		if (InitFromHits(initEndSegW)) std::cout << "pma::Track3D initialized with hit positions." << std::endl;
		else { InitFromMiddle(); std::cout << "pma::Track3D initialized in the module center." <<std::endl; }
	}
	UpdateHitsRadius();
}

bool pma::Track3D::PCEndpoints(TVector2 & start, TVector2 & stop, unsigned int view) const
{
	unsigned int nHits = 0;
	TVector2 mean(0., 0.), stdev(0., 0.), p(0., 0.);
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		if (hit->View2D() != view) continue;

		p.Set((double)hit->Wire(), (double)hit->PeakTime());
		mean += p;
		p.Set( p.X()*p.X(), p.Y()*p.Y() );
		stdev += p;
		nHits++;
	}
	if (nHits < 3) return false;
	stdev *= 1.0 / nHits;
	mean *= 1.0 / nHits;
	p = mean;
	p.Set( p.X()*p.X(), p.Y()*p.Y() );
	stdev -= p;

	double sx = stdev.X(), sy = stdev.Y();
	if (sx >= 0.0) sx = sqrt(sx);
	else sx = 0.0;
	if (sy >= 0.0) sy = sqrt(sy);
	else sy = 0.0;
	stdev.Set(sx, sy);

	double scale = 2.0 * stdev.Mod();
	double iscale = 1.0 / scale;

	unsigned int max_index = 0;
	double norm2, max_norm2 = 0.0;
	std::vector< TVector2 > data;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		if (hit->View2D() != view) continue;

		p.Set((double)hit->Wire(), (double)hit->PeakTime());
		p -= mean; p *= iscale;
		data.push_back(p);

		norm2 = p.Mod2();
		if (norm2 > max_norm2)
		{
			max_index = data.size() - 1;
			max_norm2 = norm2;
		}
	}
	if (max_norm2 == 0.0) return false;

	double y = 0.0, kappa = 1.0, prev_kappa, kchg = 1.0;
	TVector2 w(data[max_index]);

	while (kchg > 0.0001)
		for (size_t i = 0; i < data.size(); i++)
		{
			y = data[i] * w;
			w += (y/kappa) * (data[i] - y*w);

			prev_kappa = kappa;
			kappa += y*y;
			kchg = fabs((kappa - prev_kappa) / prev_kappa);
		}
	w *= 1.0 / w.Mod();

	TVector2 pt(0., 0.), v1(w), v2(w);
	v1 *= scale; v1 += mean;
	v2 *= -scale; v2 += mean;
	double bmin = 0.5, bmax = 0.5, b;
	for (size_t i = 0; i < data.size(); i++)
	{
		data[i] *= scale; data[i] += mean;
		pt = data[i];

		b = pma::GetSegmentProjVector(pt, v1, v2);
		if (b > bmax) { bmax = b; stop = pma::GetProjectionToSegment(pt, v1, v2); }
		if (b < bmin) { bmin = b; start = pma::GetProjectionToSegment(pt, v1, v2); }
	}
	return true;
}

void pma::Track3D::ClearNodes(void)
{
	for (size_t i = 0; i < fNodes.size(); i++) delete fNodes[i];
	fNodes.clear();
}

bool pma::Track3D::InitFromHits(float initEndSegW)
{
	art::ServiceHandle<util::DetectorProperties> detprop;
	art::ServiceHandle<geo::Geometry> geom;

	bool useColl = true, useInd2 = true, tryInd1 = false;

	const unsigned int nColl = NHits(geo::kZ);
	if (nColl < 3) { useColl = false; tryInd1 = true; }

	const unsigned int nInd2 = NHits(geo::kV);
	if (nInd2 < 3) { useInd2 = false; tryInd1 = true; }

	const unsigned int nInd1 = NHits(geo::kU);
	if (nInd1 < 3) tryInd1 = false;

	if (tryInd1 && !(useColl || useInd2)) return false;
	else if (!tryInd1 && !(useColl && useInd2)) return false;

	int tpc = TPCs().front(); // just the first tpc, many tpc's are ok, but need to generalize code
	
	if (Cryos().size() > 1) return false; // would need to generalize code even more if many cryostats
	int cryo = Cryos().front(); // just the first cryo

	ClearNodes();

	float wtmp = fEndSegWeight;
	fEndSegWeight = initEndSegW;

	TVector2 p00(0.0F, 0.0F);

	TVector2 ptCollStart(p00), ptCollStop(p00);
	if (useColl) useColl = PCEndpoints(ptCollStart, ptCollStop, geo::kZ);

	TVector2 ptInd2Start(p00), ptInd2Stop(p00);
	if (useInd2) useInd2 = PCEndpoints(ptInd2Start, ptInd2Stop, geo::kV);

	TVector2 ptInd1Start(p00), ptInd1Stop(p00);
	if (tryInd1 && !PCEndpoints(ptInd1Start, ptInd1Stop, geo::kU))
	{
		return false;
	}

	// endpoints for the first combination:
	TVector3 v3d_1(0., 0., 0.);
	TVector3 v3d_2(0., 0., 0.);
	// endpoints for the inverted combination
	TVector3 v3d_3(0., 0., 0.);
	TVector3 v3d_4(0., 0., 0.);

	unsigned int wireV_idx, wireZ_idx;
	double x, y, z;
	if (useColl && useInd2)
	{
		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStart.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd2Start.Y(), geo::kV, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd2Start.X());
		wireZ_idx = (unsigned int)round(ptCollStart.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kV, geo::kZ, cryo, tpc, y, z);
		v3d_1.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStop.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd2Stop.Y(), geo::kV, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd2Stop.X());
		wireZ_idx = (unsigned int)round(ptCollStop.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kV, geo::kZ, cryo, tpc, y, z);
		v3d_2.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStart.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd2Stop.Y(), geo::kV, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd2Stop.X());
		wireZ_idx = (unsigned int)round(ptCollStart.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kV, geo::kZ, cryo, tpc, y, z);
		v3d_3.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStop.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd2Start.Y(), geo::kV, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd2Start.X());
		wireZ_idx = (unsigned int)round(ptCollStop.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kV, geo::kZ, cryo, tpc, y, z);
		v3d_4.SetXYZ(x, y, z);

		std::cout << "Coll - Ind2." << std::endl;
	}
	else if (useColl && tryInd1)
	{
		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStart.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Start.X());
		wireZ_idx = (unsigned int)round(ptCollStart.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_1.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStop.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Stop.X());
		wireZ_idx = (unsigned int)round(ptCollStop.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_2.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStart.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Stop.X());
		wireZ_idx = (unsigned int)round(ptCollStart.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_3.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStop.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Start.X());
		wireZ_idx = (unsigned int)round(ptCollStop.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_4.SetXYZ(x, y, z);

		std::cout << "Coll - Ind1." << std::endl;
	}
	else if (useInd2 && tryInd1)
	{
		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Start.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Start.X());
		wireZ_idx = (unsigned int)round(ptInd2Start.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_1.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Stop.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Stop.X());
		wireZ_idx = (unsigned int)round(ptInd2Stop.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_2.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Start.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Stop.X());
		wireZ_idx = (unsigned int)round(ptInd2Start.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_3.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Stop.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireV_idx = (unsigned int)round(ptInd1Start.X());
		wireZ_idx = (unsigned int)round(ptInd2Stop.X());
		geom->IntersectionPoint(wireV_idx, wireZ_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_4.SetXYZ(x, y, z);

		std::cout << "Ind2 - Ind1." << std::endl;
	}
	else return false;

	// try objective function on the first combination
	if (fAssignedPoints.size() == 1)
	{
		TVector3* refpt = fAssignedPoints[0];
		double d1 = pma::Dist2(v3d_1, *refpt);
		double d2 = pma::Dist2(v3d_2, *refpt);
		if (d1 < d2) v3d_1 = *refpt;
		else v3d_2 = *refpt;
	}
	AddNode(v3d_1, tpc, cryo);
	AddNode(v3d_2, tpc, cryo);

	RebuildSegments();
	MakeProjection();
	UpdateHitsRadius();
	
	double g1 = Optimize(0, 0.0001F);
	std::cout << "  g1 = " << g1 << std::endl;
	//--------------------------------------------------

	// try objective function on the inverted combination
	if (fAssignedPoints.size() == 1)
	{
		TVector3* refpt = fAssignedPoints[0];
		double d1 = pma::Dist2(v3d_3, *refpt);
		double d2 = pma::Dist2(v3d_4, *refpt);
		if (d1 < d2) v3d_3 = *refpt;
		else v3d_4 = *refpt;
	}
	ClearNodes();
	AddNode(v3d_3, tpc, cryo);
	AddNode(v3d_4, tpc, cryo);

	RebuildSegments();
	MakeProjection();
	UpdateHitsRadius();

	double g2 = Optimize(0, 0.0001F);
	std::cout << "  g2 = " << g2 << std::endl;
	//--------------------------------------------------

	// compare objective functions of the two combinations
	if (g1 < g2)
	{
		ClearNodes();
		AddNode(v3d_1, tpc, cryo);
		AddNode(v3d_2, tpc, cryo);

		RebuildSegments();
		MakeProjection();
		UpdateHitsRadius();
		std::cout << "First combination good." << std::endl;
	}
	else
	{
		std::cout << "Inverted combination good." << std::endl;
	}
	//--------------------------------------------------
/*
	// try to correct if the optimization has converged to an extremely short segment
	if (((g1 > 1.0e+19) && (g2 > 1.0e+19)) ||
	    (fVertices[0]->GetDistanceTo(*(fVertices[1])) < 0.3) ||
	    ((nColl > 1) && (fSegments[0]->HitsRadius3D(T600::kViewColl) < 0.4F)) ||
	    ((nInd2 > 1) && (fSegments[0]->HitsRadius3D(T600::kViewInd2) < 0.4F)))
	{
		LOG_WARN("Poor initialization - will try again.\n");

		AF::PlHit3D* hit0_a = static_cast< AF::PlHit3D* >(front());
		AF::PlHit3D* hit0_b = NULL;
		AF::PlHit3D* hit = NULL;
		int diff, minDiff, minDrift = hit0_a->SourceHit2D().iDrift;
		for (unsigned int i = 1; i < size(); i++)
		{
			hit = static_cast< AF::PlHit3D* >((*this)[i]);
			if (hit->SourceHit2D().iDrift < minDrift)
			{
				minDrift = hit->SourceHit2D().iDrift;
				hit0_a = hit;
			}
		}
		minDiff = 5000;
		for (unsigned int i = 1; i < size(); i++)
		{
			hit = static_cast< AF::PlHit3D* >((*this)[i]);
			diff = (int)fabs(hit->SourceHit2D().iDrift - minDrift);
			if ((diff < minDiff) && (hit->View2D() != hit0_a->View2D()))
			{
				minDiff = diff; hit0_b = hit;
			}
		}

		AF::PlHit3D* hit1_a = static_cast< AF::PlHit3D* >(front());
		AF::PlHit3D* hit1_b = NULL;
		int maxDiff, maxDrift = hit1_a->SourceHit2D().iDrift;
		for (unsigned int i = 1; i < size(); i++)
		{
			hit = static_cast< AF::PlHit3D* >((*this)[i]);
			if (hit->SourceHit2D().iDrift > maxDrift)
			{
				maxDrift = hit->SourceHit2D().iDrift;
				hit1_a = hit;
			}
		}
		minDiff = 5000;
		for (unsigned int i = 1; i < size(); i++)
		{
			hit = static_cast< AF::PlHit3D* >((*this)[i]);
			diff = (int)fabs(hit->SourceHit2D().iDrift - maxDrift);
			if ((diff < minDiff) && (hit->View2D() != hit1_a->View2D()))
			{
				minDiff = diff; hit1_b = hit;
			}
		}

		if (hit0_a && hit0_b && hit1_a && hit1_b)
		{
			AF::RefPoint refpt0_a((float)hit0_a->SourceHit2D().iWire, (float)hit0_a->SourceHit2D().iDrift);
			refpt0_a.SetModuleView(fModule, hit0_a->View2D());

			AF::RefPoint refpt0_b((float)hit0_b->SourceHit2D().iWire, (float)hit0_b->SourceHit2D().iDrift);
			refpt0_b.SetModuleView(fModule, hit0_b->View2D());

			AF::RefPoint refpt1_a((float)hit1_a->SourceHit2D().iWire, (float)hit1_a->SourceHit2D().iDrift);
			refpt1_a.SetModuleView(fModule, hit1_a->View2D());

			AF::RefPoint refpt1_b((float)hit1_b->SourceHit2D().iWire, (float)hit1_b->SourceHit2D().iDrift);
			refpt1_b.SetModuleView(fModule, hit1_b->View2D());

			fVertices[0]->SetPoint3D(Point3DTool::Get3DCm(&refpt0_a, &refpt0_b));
			fVertices[1]->SetPoint3D(Point3DTool::Get3DCm(&refpt1_a, &refpt1_b));
			MakeProjection();
			UpdateHitsRadius();
			Optimize(0, 0.0001F);
		}
		else
		{
			LOG_WARN("Good hits not found.\n");
			fEndSegWeight = wtmp;
			return false;
		}
	}

	if (fVertices[0]->GetDistanceTo(*(fVertices[1])) < 0.3)
	{
		LOG_WARN("Short initial segment.\n");
		fEndSegWeight = wtmp;
		return false;
	}
*/
	fEndSegWeight = wtmp;
	return true;
}

bool pma::Track3D::InitFromRefPoints(void)
{
	return false;
}

void pma::Track3D::InitFromMiddle(void)
{
}

bool pma::Track3D::push_back(recob::Hit const* hit)
{
	for (size_t i = 0; i < fHits.size(); i++)
	{
		pma::Hit3D* hit_i = fHits[i];
		if (hit_i->fHit == hit) return false;
	}
	fHits.push_back(new pma::Hit3D(hit));
	return true;
}

void pma::Track3D::AddHits(const std::vector< recob::Hit const* >& hits)
{
	for (auto const& hit : hits) push_back(hit);
}

unsigned int pma::Track3D::NHits(unsigned int view) const
{
	unsigned int n = 0;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		if (hit->View2D() == view) n++;
	}
	return n;
}

unsigned int pma::Track3D::NEnabledHits(unsigned int view) const
{
	unsigned int n = 0;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		if (hit->IsEnabled() &&
		    ((view == geo::kUnknown) || (view == hit->View2D()))) n++;
	}
	return n;
}

std::vector< int > pma::Track3D::TPCs(void) const
{
	std::vector< int > tpc_idxs;
	for (size_t i = 0; i < size(); i++)
	{
		int tpc = (*this)[i]->TPC();

		bool found = false;
		for (size_t j = 0; j < tpc_idxs.size(); j++)
			if (tpc_idxs[j] == tpc) { found = true; break; }

		if (!found) tpc_idxs.push_back(tpc);
	}
	return tpc_idxs;
}

std::vector< int > pma::Track3D::Cryos(void) const
{
	std::vector< int > cryo_idxs;
	for (size_t i = 0; i < size(); i++)
	{
		int cryo = (*this)[i]->Cryo();

		bool found = false;
		for (size_t j = 0; j < cryo_idxs.size(); j++)
			if (cryo_idxs[j] == cryo) { found = true; break; }

		if (!found) cryo_idxs.push_back(cryo);
	}
	return cryo_idxs;
}

pma::Segment3D* pma::Track3D::NextSegment(pma::Node3D* vtx) const
{
	pma::Segment3D* seg = 0;
	unsigned int nCount = vtx->NextCount();
	unsigned int k = 0;
	while (k < nCount)
	{
		seg = static_cast< pma::Segment3D* >(vtx->Next(k));
		if (seg && (seg->Parent() == this)) return seg;
		k++;
	}
	return 0;
}

pma::Segment3D* pma::Track3D::PrevSegment(pma::Node3D* vtx) const
{
	if (vtx->Prev())
	{
		pma::Segment3D* seg = static_cast< pma::Segment3D* >(vtx->Prev());
		if (seg->Parent() == this) return seg;
	}
	return 0;
}

void pma::Track3D::AddNode(TVector3 const & p3d, unsigned int tpc, unsigned int cryo)
{
	pma::Node3D* vtx = new pma::Node3D(p3d, tpc, cryo);
	fNodes.push_back(vtx);

	if (fNodes.size() > 1) RebuildSegments();
}

bool pma::Track3D::AddNode(void)
{
	pma::Segment3D* seg;
	pma::Segment3D* maxSeg = NULL;

	size_t si = 0;
	while (si < fSegments.size())
	{
		if (!fSegments[si]->IsFrozen())
		{
			maxSeg = fSegments[si];
			break;
		}
		else si++;
	}
	if (!maxSeg) return false;

	unsigned int nHitsByView[3];
	unsigned int nHits, maxHits = 0;
	unsigned int vIndex = 0, segHits, maxSegHits = 0;
	float segLength, maxLength = maxSeg->Length();
	for (unsigned int i = si + 1; i < fNodes.size(); i++)
	{
		seg = static_cast< pma::Segment3D* >(fNodes[i]->Prev());
		if (seg->IsFrozen()) continue;

		nHitsByView[0] = seg->NEnabledHits(geo::kU);
		nHitsByView[1] = seg->NEnabledHits(geo::kV);
		nHitsByView[2] = seg->NEnabledHits(geo::kZ);
		segHits = nHitsByView[0] + nHitsByView[1] + nHitsByView[2];

		if (segHits < 15)
		{
			if ((nHitsByView[0] == 0) && ((nHitsByView[1] < 4) || (nHitsByView[2] < 4))) continue;
			if ((nHitsByView[1] == 0) && ((nHitsByView[0] < 4) || (nHitsByView[2] < 4))) continue;
			if ((nHitsByView[2] == 0) && ((nHitsByView[0] < 4) || (nHitsByView[1] < 4))) continue;
		}

		nHits = fNodes[i]->NEnabledHits() + seg->NEnabledHits() + fNodes[i-1]->NEnabledHits();
		if (nHits > maxHits)
		{
			maxHits = nHits;
			maxLength = seg->Length();
			maxSegHits = segHits;
			maxSeg = seg;
			vIndex = i;
		}
		else if (nHits == maxHits)
		{
			segLength = seg->Length();
			if (segLength > maxLength)
			{
				maxLength = segLength;
				maxSegHits = segHits;
				maxSeg = seg;
				vIndex = i;
			}
		}
	}

	if (maxSegHits > 1)
	{
		maxSeg->SortHits();

		nHitsByView[0] = maxSeg->NEnabledHits(geo::kU);
		nHitsByView[1] = maxSeg->NEnabledHits(geo::kV);
		nHitsByView[2] = maxSeg->NEnabledHits(geo::kZ);

		unsigned int maxViewIdx = 2, midViewIdx = 2;
		if ((nHitsByView[2] >= nHitsByView[1]) && (nHitsByView[1] >= nHitsByView[0])) { maxViewIdx = 2; midViewIdx = 1; }
		else if ((nHitsByView[1] >= nHitsByView[2]) && (nHitsByView[2] >= nHitsByView[0])) { maxViewIdx = 1; midViewIdx = 2; }
		else if ((nHitsByView[0] >= nHitsByView[2]) && (nHitsByView[2] >= nHitsByView[1])) { maxViewIdx = 0; midViewIdx = 2; }
		else if ((nHitsByView[2] >= nHitsByView[0]) && (nHitsByView[0] >= nHitsByView[1])) { maxViewIdx = 2; midViewIdx = 0; }
		else if ((nHitsByView[0] >= nHitsByView[1]) && (nHitsByView[1] >= nHitsByView[2])) { maxViewIdx = 0; midViewIdx = 1; }
		else if ((nHitsByView[1] >= nHitsByView[0]) && (nHitsByView[0] >= nHitsByView[2])) { maxViewIdx = 1; midViewIdx = 0; }
		if (nHitsByView[midViewIdx] < 2) midViewIdx = maxViewIdx;

		if (nHitsByView[midViewIdx] < 2) { std::cout << "AddNode(): too few hits." << std::endl; return false; }

		unsigned int mHits[3] = { 0, 0, 0 };
		unsigned int halfIndex = (nHitsByView[midViewIdx] >> 1) - 1;
		unsigned int n = 0, i = 0, i0 = 0, i1 = 0;
		while (i < maxSeg->NHits() - 1)
		{
			if (maxSeg->Hit(i).IsEnabled())
			{
				mHits[maxSeg->Hit(i).View2D()]++;
				if (maxSeg->Hit(i).View2D() == midViewIdx)
				{
					if (n == halfIndex) break;
					n++;
				}
			}
			i++;
		}
		i0 = i; i++;
		while ((i < maxSeg->NHits()) && !((maxSeg->Hit(i).View2D() == midViewIdx) && maxSeg->Hit(i).IsEnabled()))
		{
			i++;
		}
		i1 = i;

		if (!nHitsByView[0])
		{
			if (nHitsByView[1] && (mHits[1] < 2)) { std::cout << "AddNode(): low Ind2 hits." << std::endl; return false; }
			if (nHitsByView[2] && (mHits[2] < 2)) { std::cout << "AddNode(): low Coll hits." << std::endl; return false; }
		}

		maxSeg->SetProjection(maxSeg->Hit(i0));
		maxSeg->SetProjection(maxSeg->Hit(i1));

		unsigned int tpc = maxSeg->Hit(i0).TPC();
		unsigned int cryo = maxSeg->Hit(i0).Cryo();

		pma::Node3D* p = new pma::Node3D((maxSeg->Hit(i0).Point3D() + maxSeg->Hit(i1).Point3D()) * 0.5, tpc, cryo);
		fNodes.insert(fNodes.begin() + vIndex, p);

		maxSeg->AddNext(fNodes[vIndex]);

		seg = new pma::Segment3D(this, fNodes[vIndex], fNodes[vIndex + 1]);
		fSegments.insert(fSegments.begin() + vIndex, seg);

		return true;
	}
	else return false;
}

bool pma::Track3D::HasRefPoint(TVector3* p) const
{
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
		if (fAssignedPoints[i] == p) return true;
	return false;
}

double pma::Track3D::GetObjFunction(bool suppressPenalty) const
{
	double sum = 0.0;
	float p = fPenaltyValue;
	if (suppressPenalty) p *= 0.33F;
	for (size_t i = 0; i < fNodes.size(); i++)
	{
		sum += fNodes[i]->GetObjFunction(p, fEndSegWeight);
	}
	return sum / fNodes.size();
}

double pma::Track3D::Optimize(int newVertices, double eps, bool selAllHits )
{
	if (!fNodes.size()) { std::cout << "Track3D not initialized." << std:: endl; return 0.0; }

	UpdateParams();
	double g0 = GetObjFunction(), g1 = 0.0;
	if (g0 == 0.0) return g0;

	//std::cout << "objective function: " << g0 << std::endl;

	bool stop = false;
	fMinSegStop = fSegments.size();
	fMaxSegStop = (int)(fMaxSegStopFactor * size());
	do
	{
		bool stepDone = true;
		unsigned int stepIter = 0;
		do
		{
			double gstep = 1.0;
			unsigned int iter = 0;
			while ((gstep > eps) && (iter < 60))
			{
				MakeProjection();

				UpdateParams();

				for (unsigned int j = 0; j < fNodes.size(); j++)
				{
					fNodes[j]->Optimize(fPenaltyValue, fEndSegWeight);
				}

				g1 = g0;
				g0 = GetObjFunction();
				//std::cout << "obj fn: " << g0 << std::endl;
				if (g0 == 0.0F) { MakeProjection(); break; }
				gstep = fabs(g0 - g1) / g0;
				iter++;
			}

			stepIter++;
			if (fNodes.size() > 2)
			{
				stepDone = !(CheckEndSegment(pma::Track3D::kEnd) ||
				             CheckEndSegment(pma::Track3D::kBegin));
			}
		} while (!stepDone && (stepIter < 5));

		if (selAllHits && (size() / fNodes.size() < 300)) { SelectHits(); selAllHits = false; }
		switch (newVertices)
		{
			case 0: stop = true; break; // just optimize existing vertices

			case -1: // grow and optimize until automatic stop condition
				std::cout << "optimized segments: " << fSegments.size() << std::endl;
				if ((fSegments.size() >= fSegStopValue) ||
				    (fSegments.size() >= fMaxSegStop))
				{
					stop = true;
				}
				else
				{
					if (!AddNode()) stop = true;
				}
				break;

			default: // grow and optimize until fixed number of vertices is added
				//std::cout << "vertices to add:" << newVertices << std::endl;

				if (newVertices > 12)
				{
					if (AddNode()) { MakeProjection(); newVertices--; }
					else { std::cout << "stop (3)" << std::endl; stop = true; break; }

					if (AddNode())
					{
						MakeProjection(); newVertices--;
						if (AddNode()) newVertices--;
					}
				}
				else if (newVertices > 4)
				{
					if (AddNode()) { MakeProjection(); newVertices--; }
					else { std::cout << "stop (2)" << std::endl; stop = true; break; }

					if (AddNode()) newVertices--;
				}
				else
				{
					if (AddNode()) { newVertices--; }
					else { std::cout << "stop (1)" << std::endl; stop = true; break; }
				}
				break;
		}

	} while (!stop);
	//std::cout << "Done (optimized segments: " << fSegments.size() << ")." << std::endl << std::endl;

	MakeProjection();
	return GetObjFunction();
}

void pma::Track3D::RebuildSegments(void)
{
	for (size_t i = 0; i < fSegments.size(); i++) delete fSegments[i];
	fSegments.clear();

	for (size_t i = 1; i < fNodes.size(); i++)
	{
		pma::Segment3D* s = new pma::Segment3D(this, fNodes[i - 1], fNodes[i]);
		fSegments.push_back(s);
	}
}

pma::Element3D* pma::Track3D::GetNearestElement(
	const TVector2& p2d, unsigned int view) const
{
	pma::Element3D* pe_min = fNodes.front();
	double dist, min_dist = pe_min->GetDistance2To(p2d, view);
	for (size_t i = 1; i < fNodes.size(); i++)
	{
		dist = fNodes[i]->GetDistance2To(p2d, view);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fNodes[i];
		}
	}
	for (size_t i = 0; i < fSegments.size(); i++)
	{
		dist = fSegments[i]->GetDistance2To(p2d, view);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fSegments[i];
		}
	}
	return pe_min;
}

pma::Element3D* pma::Track3D::GetNearestElement(const TVector3& p3d) const
{
	pma::Element3D* pe_min = fNodes.front();
	double dist, min_dist = pe_min->GetDistance2To(p3d);
	for (size_t i = 1; i < fNodes.size(); i++)
	{
		dist = fNodes[i]->GetDistance2To(p3d);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fNodes[i];
		}
	}
	for (size_t i = 0; i < fSegments.size(); i++)
	{
		dist = fSegments[i]->GetDistance2To(p3d);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fSegments[i];
		}
	}
	return pe_min;
}

void pma::Track3D::SortHits(void)
{
	std::vector< pma::Hit3D* > hits_tmp;
	pma::Node3D* vtx = fNodes.front();
	pma::Segment3D* seg = NextSegment(vtx);
	while (vtx)
	{
		vtx->SortHits();
		for (size_t i = 0; i < vtx->NHits(); i++)
		{
			pma::Hit3D* h3d = &(vtx->Hit(i));
			for (size_t j = 0; j < size(); j++)
				if ((*this)[j] == h3d)
				{
					hits_tmp.push_back(h3d);
					break;
				}
		}

		if (seg)
		{
			seg->SortHits();
			for (size_t i = 0; i < seg->NHits(); i++)
			{
				pma::Hit3D* h3d = &(seg->Hit(i));
				for (size_t j = 0; j < size(); j++)
					if ((*this)[j] == h3d)
					{
						hits_tmp.push_back(h3d);
						break;
					}
			}

			vtx = static_cast< pma::Node3D* >(seg->Next());
			seg = NextSegment(vtx);
		}
		else break;
	}

	if (size() == hits_tmp.size())
	{
		for (size_t i = 0; i < size(); i++) (*this)[i] = hits_tmp[i];
	}
	else std::cout << "Hit sorting problem." << std::endl;
}

void pma::Track3D::SelectHits(float fraction)
{
	if (fraction < 0.0F) fraction = 0.0F;
	if (fraction > 1.0F) fraction = 1.0F;
	if (fraction < 1.0F) std::sort(fHits.begin(), fHits.end(), pma::bTrajectory3DDistLess());

	unsigned int nHitsColl = (unsigned int)(fraction * NHits(geo::kZ));
	unsigned int nHitsInd2 = (unsigned int)(fraction * NHits(geo::kV));
	unsigned int nHitsInd1 = (unsigned int)(fraction * NHits(geo::kU));
	unsigned int coll = 0, ind2 = 0, ind1 = 0;

	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		if (fraction < 1.0F)
		{
			hit->SetEnabled(false);
			switch (hit->View2D())
			{
				case geo::kZ: if (coll++ < nHitsColl) hit->SetEnabled(true); break;
				case geo::kV: if (ind2++ < nHitsInd2) hit->SetEnabled(true); break;
				case geo::kU: if (ind1++ < nHitsInd1) hit->SetEnabled(true); break;
			}
		}
		else hit->SetEnabled(true);
	}

	pma::Element3D* pe = FirstElement();
	for (unsigned int i = 0; i < pe->NHits(); i++)
		pe->Hit(i).SetEnabled(true);

	pe = LastElement();
	for (unsigned int i = 0; i < pe->NHits(); i++)
		pe->Hit(i).SetEnabled(true);
}

void pma::Track3D::MakeProjection(void)
{
	for (size_t i = 0; i < fNodes.size(); i++) fNodes[i]->ClearAssigned(this);
	for (size_t i = 0; i < fSegments.size(); i++) fSegments[i]->ClearAssigned(this);

	pma::Element3D* pe = 0;

	for (size_t i = 0; i < size(); i++) // assign hits to nodes/segments
	{
		pma::Hit3D* hit = (*this)[i];
		pe = GetNearestElement(hit->Point2D(), hit->View2D());
		pe->AddHit(hit);
	}

	for (size_t i = 0; i < fAssignedPoints.size(); i++) // assign ref points to nodes/segments
	{
		pe = GetNearestElement(*(fAssignedPoints[i]));
		pe->AddPoint(fAssignedPoints[i]);
	}

	// move hits to segments if not branching
	if (!(fNodes.front()->Prev()) && (fNodes.front()->NextCount() == 1))
	{
		for (size_t i = 0; i < fNodes.front()->NHits(); i++)
			fSegments.front()->AddHit(&(fNodes.front()->Hit(i)));
		fNodes.front()->ClearAssigned();
	}

	if (fNodes.back()->NextCount() == 0)
	{
		for (size_t i = 0; i < fNodes.back()->NHits(); i++)
			fSegments.back()->AddHit(&(fNodes.back()->Hit(i)));
		fNodes.back()->ClearAssigned();
	}

	for (unsigned int i = 0; i < fNodes.size(); i++) fNodes[i]->UpdateHitParams();
	for (unsigned int i = 0; i < fSegments.size(); i++) fSegments[i]->UpdateHitParams();
}

void pma::Track3D::UpdateProjection(void)
{
	for (size_t i = 0; i < fNodes.size(); i++) fNodes[i]->UpdateProjection();
	for (size_t i = 0; i < fSegments.size(); i++) fSegments[i]->UpdateProjection();
}

double pma::Track3D::AverageDist2(void) const
{
	double sum = 0.0;
	unsigned int count = 0;
	pma::Node3D* vtx = fNodes.front();
	pma::Segment3D* seg = NextSegment(vtx);
	while (vtx)
	{
		sum += vtx->SumDist2();
		count += vtx->NEnabledHits();
		if (seg)
		{
			sum += seg->SumDist2();
			count += seg->NEnabledHits();
			vtx = static_cast< pma::Node3D* >(seg->Next());
			seg = NextSegment(vtx);
		}
		else break;
	}
	return sum / count;
}

void pma::Track3D::UpdateParams(void)
{
	size_t n = size();
	if (!n) n = 1;

	float nCubeRoot = pow((double)n, 1.0/3.0);
	float avgDist2Root = sqrt(AverageDist2());	

	fPenaltyValue = fPenaltyFactor * pow((double)fSegments.size(), 1.8) * avgDist2Root / (fHitsRadius * nCubeRoot);

	fSegStopValue = (int)(fSegStopFactor * nCubeRoot * fHitsRadius / avgDist2Root);
	if (fSegStopValue < fMinSegStop) fSegStopValue = fMinSegStop;
}

bool pma::Track3D::SwapVertices(size_t v0, size_t v1)
{
	if (v0 == v1) return false;

	if (v0 > v1)
	{
		size_t vx = v0;
		v0 = v1; v1 = vx;
	}

	pma::Node3D* vtmp;
	if (v1 - v0 == 1)
	{
		pma::Segment3D* midSeg = NextSegment(fNodes[v0]);
		pma::Segment3D* prevSeg = PrevSegment(fNodes[v0]);
		pma::Segment3D* nextSeg = NextSegment(fNodes[v1]);

		fNodes[v1]->RemoveNext(nextSeg);
		midSeg->Disconnect();

		vtmp = fNodes[v0];
		fNodes[v0] = fNodes[v1];
		fNodes[v1] = vtmp;

		if (prevSeg) prevSeg->AddNext(fNodes[v0]);
		fNodes[v0]->AddNext(midSeg);
		midSeg->AddNext(fNodes[v1]);
		if (nextSeg) fNodes[v1]->AddNext(nextSeg);

		return false;
	}
	else
	{
		vtmp = fNodes[v0];
		fNodes[v0] = fNodes[v1];
		fNodes[v1] = vtmp;
		return true;
	}
}

bool pma::Track3D::CheckEndSegment(pma::Track3D::ETrackEnd endCode)
{
	unsigned int v1, v2;
	switch (endCode)
	{
		case pma::Track3D::kBegin:
			if (fSegments.front()->IsFrozen()) return false;
			if (fNodes.front()->NextCount() > 1) return false;
			v1 = 0; v2 = 1; break;
		case pma::Track3D::kEnd:
			if (fSegments.back()->IsFrozen()) return false;
			if (fNodes.back()->NextCount() > 1) return false;
			v1 = fNodes.size() - 1;
			v2 = fNodes.size() - 2;
			break;
		default: return false;
	}

	double g1, g0 = GetObjFunction();

	if (SwapVertices(v1, v2)) RebuildSegments();
	MakeProjection();
	g1 = GetObjFunction();

	if (g1 >= g0)
	{
		if (SwapVertices(v1, v2)) RebuildSegments();
		MakeProjection();
		return false;
	}
	else return true;
}

void pma::Track3D::UpdateHitsRadius(void)
{
	std::vector< pma::Hit3D* > hitsColl, hitsInd1, hitsInd2;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		switch (hit->View2D())
		{
			case geo::kZ: hitsColl.push_back(hit); break;
			case geo::kV: hitsInd2.push_back(hit); break;
			case geo::kU: hitsInd1.push_back(hit); break;
		}
	}
	fHitsRadius = pma::GetHitsRadius2D(hitsColl, true);
	float r = pma::GetHitsRadius2D(hitsInd2, true);
	if (r > fHitsRadius) fHitsRadius = r;
	r = pma::GetHitsRadius2D(hitsInd1, true);
	if (r > fHitsRadius) fHitsRadius = r;
}

