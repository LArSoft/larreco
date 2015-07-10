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
#include "RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::Track3D::Track3D(void) :
	fMaxHitsPerSeg(70),
	fPenaltyFactor(1.0F),
	fMaxSegStopFactor(8.0F),

	fSegStopValue(2), fMinSegStop(2), fMaxSegStop(2),
	fSegStopFactor(0.2F),
	fPenaltyValue(0.1F),
	fEndSegWeight(0.05F),
	fHitsRadius(1.0F)
{
}

pma::Track3D::Track3D(const Track3D& src) :
	fMaxHitsPerSeg(src.fMaxHitsPerSeg),
	fPenaltyFactor(src.fPenaltyFactor),
	fMaxSegStopFactor(src.fMaxSegStopFactor),

	fSegStopValue(src.fSegStopValue),
	fMinSegStop(src.fMinSegStop),
	fMaxSegStop(src.fMaxSegStop),
	fSegStopFactor(src.fSegStopFactor),
	fPenaltyValue(src.fPenaltyValue),
	fEndSegWeight(src.fEndSegWeight),
	fHitsRadius(src.fHitsRadius)
{
	for (auto const& hit : src.fHits) fHits.push_back(new pma::Hit3D(*hit));
	for (auto const& point : src.fAssignedPoints) fAssignedPoints.push_back(new TVector3(*point));
	for (auto const& node : src.fNodes) fNodes.push_back(new pma::Node3D(node->Point3D(), node->TPC(), node->Cryo()));

	RebuildSegments();
	MakeProjection();
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
	int tpc = TPCs().front(); // just the first tpc, many tpc's are ok, but need to generalize code
	
	if (Cryos().size() > 1)
	{
		mf::LogError("pma::Track3D") << "Only one cryostat for now, please.";
		return; // would need to generalize code even more if many cryostats
	}
	int cryo = Cryos().front(); // just the first cryo
	
	if (InitFromRefPoints(tpc, cryo)) mf::LogVerbatim("pma::Track3D") << "Track initialized with 3D reference points.";
	else
	{
		if (InitFromHits(tpc, cryo, initEndSegW)) mf::LogVerbatim("pma::Track3D") << "Track initialized with hit positions.";
		else { InitFromMiddle(tpc, cryo); mf::LogVerbatim("pma::Track3D") << "Track initialized in the module center."; }
	}
	UpdateHitsRadius();
}

bool pma::Track3D::PCEndpoints(TVector2 & start, TVector2 & stop,
	unsigned int view, double wpitch, double dpitch) const
{
	unsigned int nHits = 0;
	TVector2 mean(0., 0.), stdev(0., 0.), p(0., 0.);
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = (*this)[i];
		if (hit->View2D() != view) continue;

		p = pma::WireDriftToCm(hit->Wire(), hit->PeakTime(), hit->View2D(), hit->TPC(), hit->Cryo());
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

		p = pma::WireDriftToCm(hit->Wire(), hit->PeakTime(), hit->View2D(), hit->TPC(), hit->Cryo());
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
	start = pma::CmToWireDrift(start.X(), start.Y(), front()->View2D(), front()->TPC(), front()->Cryo());
	stop = pma::CmToWireDrift(stop.X(), stop.Y(), front()->View2D(), front()->TPC(), front()->Cryo());
	return true;
}

void pma::Track3D::ClearNodes(void)
{
	for (size_t i = 0; i < fNodes.size(); i++) delete fNodes[i];
	fNodes.clear();
}

bool pma::Track3D::InitFromHits(int tpc, int cryo, float initEndSegW)
{
	art::ServiceHandle<util::DetectorProperties> detprop;
	art::ServiceHandle<geo::Geometry> geom;

	float wtmp = fEndSegWeight;
	fEndSegWeight = initEndSegW;

	//TVector2 p00(0.0F, 0.0F);

	// endpoints for the first combination:
	TVector3 v3d_1(0., 0., 0.), v3d_2(0., 0., 0.);

	// endpoints for the inverted combination
	//TVector3 v3d_3(0., 0., 0.), v3d_4(0., 0., 0.);

	//unsigned int wireU_idx, wireV_idx, wireZ_idx;
	double x, y, z;

/*
//  another way of initialization ******** NEED TO CORRECT wire/drift <-> cm *********

	bool useColl = true, useInd2 = true, tryInd1 = false;

	const unsigned int nColl = NHits(geo::kZ);
	if (nColl < 3) { useColl = false; tryInd1 = true; }

	const unsigned int nInd2 = NHits(geo::kV);
	if (nInd2 < 3) { useInd2 = false; tryInd1 = true; }

	const unsigned int nInd1 = NHits(geo::kU);
	if (nInd1 < 3) tryInd1 = false;

	if (tryInd1 && !(useColl || useInd2)) return false;
	else if (!tryInd1 && !(useColl && useInd2)) return false;

	double wirePitch, driftPitch = detprop->GetXTicksCoefficient(tpc, cryo);

	ClearNodes();

	TVector2 ptCollStart(p00), ptCollStop(p00);
	wirePitch = geom->TPC(tpc, cryo).Plane(geo::kZ).WirePitch();
	if (useColl) useColl = PCEndpoints(ptCollStart, ptCollStop, geo::kZ, wirePitch, driftPitch);

	TVector2 ptInd2Start(p00), ptInd2Stop(p00);
	wirePitch = geom->TPC(tpc, cryo).Plane(geo::kV).WirePitch();
	if (useInd2) useInd2 = PCEndpoints(ptInd2Start, ptInd2Stop, geo::kV, wirePitch, driftPitch);

	TVector2 ptInd1Start(p00), ptInd1Stop(p00);
	wirePitch = geom->TPC(tpc, cryo).Plane(geo::kU).WirePitch();
	if (tryInd1 && !PCEndpoints(ptInd1Start, ptInd1Stop, geo::kU, wirePitch, driftPitch))
	{
		return false;
	}

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

		mf::LogVerbatim("pma::Track3D") << "Initialization from Coll - Ind2.";
	}
	else if (useColl && tryInd1)
	{
		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStart.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Start.X());
		wireZ_idx = (unsigned int)round(ptCollStart.X());
		geom->IntersectionPoint(wireU_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_1.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStop.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Stop.X());
		wireZ_idx = (unsigned int)round(ptCollStop.X());
		geom->IntersectionPoint(wireU_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_2.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStart.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Stop.X());
		wireZ_idx = (unsigned int)round(ptCollStart.X());
		geom->IntersectionPoint(wireU_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_3.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptCollStop.Y(), geo::kZ, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Start.X());
		wireZ_idx = (unsigned int)round(ptCollStop.X());
		geom->IntersectionPoint(wireU_idx, wireZ_idx, geo::kU, geo::kZ, cryo, tpc, y, z);
		v3d_4.SetXYZ(x, y, z);

		mf::LogVerbatim("pma::Track3D") << "Initialization from Coll - Ind1.";
	}
	else if (useInd2 && tryInd1)
	{
		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Start.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Start.X());
		wireV_idx = (unsigned int)round(ptInd2Start.X());
		geom->IntersectionPoint(wireU_idx, wireV_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_1.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Stop.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Stop.X());
		wireV_idx = (unsigned int)round(ptInd2Stop.X());
		geom->IntersectionPoint(wireU_idx, wireV_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_2.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Start.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Stop.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Stop.X());
		wireV_idx = (unsigned int)round(ptInd2Start.X());
		geom->IntersectionPoint(wireU_idx, wireV_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_3.SetXYZ(x, y, z);

		x = 0.5 * (
			detprop->ConvertTicksToX(ptInd2Stop.Y(), geo::kV, tpc, cryo) +
			detprop->ConvertTicksToX(ptInd1Start.Y(), geo::kU, tpc, cryo));
		wireU_idx = (unsigned int)round(ptInd1Start.X());
		wireV_idx = (unsigned int)round(ptInd2Stop.X());
		geom->IntersectionPoint(wireU_idx, wireV_idx, geo::kU, geo::kV, cryo, tpc, y, z);
		v3d_4.SetXYZ(x, y, z);

		mf::LogVerbatim("pma::Track3D") << "Initialization from Ind2 - Ind1.";
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
	
	//double g1 = GetObjFunction(true);
	double g1 = Optimize(0, 0.01F);
	mf::LogVerbatim("pma::Track3D") << "  g1 = " << g1;
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

	//double g2 = GetObjFunction(true);
	double g2 = Optimize(0, 0.01F);
	mf::LogVerbatim("pma::Track3D") << "  g2 = " << g2;
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
		mf::LogVerbatim("pma::Track3D") << "First combination good.";
	}
	else
	{
		mf::LogVerbatim("pma::Track3D") << "Inverted combination good.";
	}
	//--------------------------------------------------
*/

/*
	// try to correct if the optimization has converged to an extremely short segment
	if (((g1 > 1.0e+19) && (g2 > 1.0e+19)) ||
	    (fVertices[0]->GetDistanceTo(*(fVertices[1])) < 0.3) ||
	    ((nColl > 1) && (fSegments[0]->HitsRadius3D(T600::kViewColl) < 0.4F)) ||
	    ((nInd2 > 1) && (fSegments[0]->HitsRadius3D(T600::kViewInd2) < 0.4F)))
*/
	{
		pma::Hit3D* hit0_a = front();
		pma::Hit3D* hit0_b = 0;
		pma::Hit3D* hit = 0;
		float diff, minDiff, minX = fabs(detprop->ConvertTicksToX(hit0_a->PeakTime(), hit0_a->View2D(), tpc, cryo));
		for (size_t i = 1; i < size(); i++)
		{
			hit = (*this)[i];
			x = fabs(detprop->ConvertTicksToX(hit->PeakTime(), hit->View2D(), tpc, cryo));
			if (x < minX)
			{
				minX = x; hit0_a = hit;
			}
		}
		minDiff = 5000;
		for (size_t i = 0; i < size(); i++)
		{
			hit = (*this)[i];
			x = fabs(detprop->ConvertTicksToX(hit->PeakTime(), hit->View2D(), tpc, cryo));
			diff = fabs(x - minX);
			if ((diff < minDiff) && (hit->View2D() != hit0_a->View2D()))
			{
				minDiff = diff; hit0_b = hit;
			}
		}

		pma::Hit3D* hit1_a = front();
		pma::Hit3D* hit1_b = 0;
		float maxX = fabs(detprop->ConvertTicksToX(hit1_a->PeakTime(), hit1_a->View2D(), tpc, cryo));
		for (size_t i = 1; i < size(); i++)
		{
			hit = (*this)[i];
			x = fabs(detprop->ConvertTicksToX(hit->PeakTime(), hit->View2D(), tpc, cryo));
			if (x > maxX)
			{
				maxX = x; hit1_a = hit;
			}
		}
		minDiff = 5000;
		for (size_t i = 0; i < size(); i++)
		{
			hit = (*this)[i];
			x = fabs(detprop->ConvertTicksToX(hit->PeakTime(), hit->View2D(), tpc, cryo));
			diff = fabs(x - maxX);
			if ((diff < minDiff) && (hit->View2D() != hit1_a->View2D()))
			{
				minDiff = diff; hit1_b = hit;
			}
		}

		if (hit0_a && hit0_b && hit1_a && hit1_b)
		{
			x = 0.5 * (
				detprop->ConvertTicksToX(hit0_a->PeakTime(), hit0_a->View2D(), tpc, cryo) +
				detprop->ConvertTicksToX(hit0_b->PeakTime(), hit0_b->View2D(), tpc, cryo));
			geom->IntersectionPoint(hit0_a->Wire(), hit0_b->Wire(),
				hit0_a->View2D(), hit0_b->View2D(), cryo, tpc, y, z);
			v3d_1.SetXYZ(x, y, z);

			x = 0.5 * (
				detprop->ConvertTicksToX(hit1_a->PeakTime(), hit1_a->View2D(), tpc, cryo) +
				detprop->ConvertTicksToX(hit1_b->PeakTime(), hit1_b->View2D(), tpc, cryo));
			geom->IntersectionPoint(hit1_a->Wire(), hit1_b->Wire(),
				hit1_a->View2D(), hit1_b->View2D(), cryo, tpc, y, z);
			v3d_2.SetXYZ(x, y, z);

			ClearNodes();
			AddNode(v3d_1, tpc, cryo);
			AddNode(v3d_2, tpc, cryo);

			MakeProjection();
			UpdateHitsRadius();
			Optimize(0, 0.01F);
		}
		else
		{
			mf::LogVerbatim("pma::Track3D") << "Good hits not found.";
			fEndSegWeight = wtmp;
			return false;
		}
	}

	if (sqrt(pma::Dist2(fNodes.front()->Point3D(), fNodes.back()->Point3D())) < 0.3)
	{
		mf::LogVerbatim("pma::Track3D") << "Short initial segment.";
		fEndSegWeight = wtmp;
		return false;
	}

	fEndSegWeight = wtmp;
	return true;
}

bool pma::Track3D::InitFromRefPoints(int tpc, int cryo)
{
	if (fAssignedPoints.size() < 2) return false;

	ClearNodes();

	TVector3 mean(0., 0., 0.), stdev(0., 0., 0.), p(0., 0., 0.);
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
	{
		p = *(fAssignedPoints[i]);
		mean += p;
		p.SetXYZ( p.X()*p.X(), p.Y()*p.Y(), p.Z()*p.Z() );
		stdev += p;
	}
	stdev *= 1.0 / fAssignedPoints.size();
	mean *= 1.0 / fAssignedPoints.size();
	p = mean;
	p.SetXYZ( p.X()*p.X(), p.Y()*p.Y(), p.Z()*p.Z() );
	stdev -= p;

	double sx = stdev.X(), sy = stdev.Y(), sz = stdev.Z();
	if (sx >= 0.0) sx = sqrt(sx);
	else sx = 0.0;
	if (sy >= 0.0) sy = sqrt(sy);
	else sy = 0.0;
	if (sz >= 0.0) sz = sqrt(sz);
	else sz = 0.0;
	stdev.SetXYZ(sx, sy, sz);

	double scale = 2.0 * stdev.Mag();
	double iscale = 1.0 / scale;

	size_t max_index = 0;
	double norm2, max_norm2 = 0.0;
	std::vector< TVector3 > data;
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
	{
		p = *(fAssignedPoints[i]);
		p -= mean;
		p *= iscale;
		norm2 = p.Mag2();
		if (norm2 > max_norm2)
		{
			max_norm2 = norm2;
			max_index = i;
		}
		data.push_back(p);
	}

	double y = 0.0, kappa = 1.0, prev_kappa, kchg = 1.0;
	TVector3 w(data[max_index]);

	while (kchg > 0.0001)
		for (size_t i = 0; i < data.size(); i++)
		{
			y = (data[i] * w);
			w += (y/kappa) * (data[i] - y*w);

			prev_kappa = kappa;
			kappa += y*y;
			kchg = fabs((kappa - prev_kappa) / prev_kappa);
		}
	w *= 1.0 / w.Mag();

	TVector3 v1(w), v2(w);
	v1 *= scale; v1 += mean;
	v2 *= -scale; v2 += mean;
	std::sort(fAssignedPoints.begin(), fAssignedPoints.end(), pma::bSegmentProjLess(v1, v2));
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
	{
		AddNode(*(fAssignedPoints[i]), tpc, cryo);
	}

	RebuildSegments();
	MakeProjection();

	if (size()) UpdateHitsRadius();

	Optimize(0, 0.01F);

	return true;
}

void pma::Track3D::InitFromMiddle(int tpc, int cryo)
{
	art::ServiceHandle<geo::Geometry> geom;

	const auto& tpcGeo = geom->TPC(tpc, cryo);

	double minX = tpcGeo.MinX(), maxX = tpcGeo.MaxX();
	double minY = tpcGeo.MinY(), maxY = tpcGeo.MaxY();
	double minZ = tpcGeo.MinZ(), maxZ = tpcGeo.MaxZ();

	TVector3 v3d_1(0.5 * (minX + maxX), 0.5 * (minY + maxY), 0.5 * (minZ + maxZ));
	TVector3 v3d_2(v3d_1);

	TVector3 shift(5.0, 5.0, 5.0);
	v3d_1 += shift;
	v3d_2 -= shift;

	ClearNodes();
	AddNode(v3d_1, tpc, cryo);
	AddNode(v3d_2, tpc, cryo);

	MakeProjection();
	UpdateHitsRadius();

	Optimize(0, 0.01F);
}

int pma::Track3D::index_of(const pma::Hit3D* hit) const
{
	for (size_t i = 0; i < size(); i++)
		if (fHits[i] == hit) return (int)i;
	return -1;
}

bool pma::Track3D::push_back(art::Ptr< recob::Hit > hit)
{
	for (auto const& trk_hit : fHits)
	{
		if (trk_hit->fHit == hit) return false;
	}
	fHits.push_back(new pma::Hit3D(hit));
	return true;
}

void pma::Track3D::AddHits(const std::vector< art::Ptr<recob::Hit> >& hits)
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

std::vector< unsigned int > pma::Track3D::TPCs(void) const
{
	std::vector< unsigned int > tpc_idxs;
	for (size_t i = 0; i < size(); i++)
	{
		unsigned int tpc = (*this)[i]->TPC();

		bool found = false;
		for (size_t j = 0; j < tpc_idxs.size(); j++)
			if (tpc_idxs[j] == tpc) { found = true; break; }

		if (!found) tpc_idxs.push_back(tpc);
	}
	return tpc_idxs;
}

std::vector< unsigned int > pma::Track3D::Cryos(void) const
{
	std::vector< unsigned int > cryo_idxs;
	for (size_t i = 0; i < size(); i++)
	{
		unsigned int cryo = (*this)[i]->Cryo();

		bool found = false;
		for (size_t j = 0; j < cryo_idxs.size(); j++)
			if (cryo_idxs[j] == cryo) { found = true; break; }

		if (!found) cryo_idxs.push_back(cryo);
	}
	return cryo_idxs;
}

void pma::Track3D::InternalFlip(std::vector< pma::Track3D* >& toSort)
{
	bool branching = false;
	for (size_t i = 0; i < fNodes.size() - 1; i++)
	{
		if (fNodes[i]->NextCount() > 1)
		{
			for (size_t j = 0; j < fNodes[i]->NextCount(); j++)
			{
				pma::Segment3D* s = static_cast< pma::Segment3D* >(fNodes[i]->Next(j));
				if (s->Parent() != this) toSort.push_back(s->Parent());
			}
			branching = true;
		}
	}
	if (fNodes.back()->NextCount())
	{
		for (size_t j = 0; j < fNodes.back()->NextCount(); j++)
		{
			pma::Segment3D* s = static_cast< pma::Segment3D* >(fNodes.back()->Next(j));
			toSort.push_back(s->Parent());
		}
		branching = true;
	}

	if (fNodes.front()->Prev())
	{
		pma::Segment3D* s = static_cast< pma::Segment3D* >(fNodes.front()->Prev());
		toSort.push_back(s->Parent());
		s->Parent()->InternalFlip(toSort);
	}

	if (branching) mf::LogWarning("pma::Track3D") << "Branched track flipped.";

	std::reverse(fNodes.begin(), fNodes.end());
	toSort.push_back(this);
	RebuildSegments();
}

void pma::Track3D::Flip(void)
{
	std::vector< pma::Track3D* > toSort;
	InternalFlip(toSort);
	toSort.push_back(this);

	for (size_t t = 0; t < toSort.size(); t++)
	{
		bool sorted = false;
		for (size_t u = 0; u < t; u++)
			if (toSort[u] == toSort[t]) { sorted = true; break; }
		if (!sorted)
		{
			toSort[t]->MakeProjection();
			toSort[t]->SortHits();
		}
	}
}

double pma::Track3D::TestHitsMse(const std::vector< art::Ptr<recob::Hit> >& hits, bool normalized) const
{
	if (!hits.size())
	{
		mf::LogWarning("pma::Track3D") << "TestHitsMse(): Empty cluster.";
		return -1.0;
	}

	double mse = 0.0;
	for (auto const & h : hits)
	{
		unsigned int cryo = h->WireID().Cryostat;
		unsigned int tpc = h->WireID().TPC;
		unsigned int view = h->WireID().Plane;
		unsigned int wire = h->WireID().Wire;
		float drift = h->PeakTime();

		mse += Dist2(pma::WireDriftToCm(wire, drift, view, tpc, cryo), view);
	}
	if (normalized) return mse / hits.size();
	else return mse;
}

unsigned int pma::Track3D::TestHits(const std::vector< art::Ptr<recob::Hit> >& hits, double dist) const
{
	if (!hits.size())
	{
		mf::LogWarning("pma::Track3D") << "TestHitsMse(): Empty cluster.";
		return 0;
	}

	double d2 = dist * dist;
	unsigned int nhits = 0;
	for (auto const & h : hits)
	{
		unsigned int cryo = h->WireID().Cryostat;
		unsigned int tpc = h->WireID().TPC;
		unsigned int view = h->WireID().Plane;
		unsigned int wire = h->WireID().Wire;
		float drift = h->PeakTime();

		if (Dist2(pma::WireDriftToCm(wire, drift, view, tpc, cryo), view) < d2) nhits++;
	}
	return nhits;
}

double pma::Track3D::Length(size_t start, size_t stop, size_t step) const
{
	if (size() < 2) return 0.0;

	if (start > stop)
	{
		size_t tmp = stop; stop = start; start = tmp;
	}
	if (start >= size() - 1) return 0.0;
	if (stop >= size()) stop = size() - 1;

	double result = 0.0;
	pma::Hit3D *h1 = 0, *h0 = fHits[start];

	if (!step) step = 1;
	size_t i = start + step;
	while (i <= stop)
	{
		h1 = fHits[i];
		result += sqrt(pma::Dist2(h1->Point3D(), h0->Point3D()));
		h0 = h1;
		i += step;
	}
	if (i - step < stop) // last step jumped beyond the stop point
	{                    // so need to add the last short segment
		result += sqrt(pma::Dist2(h0->Point3D(), back()->Point3D()));
	}
	return result;
}

int pma::Track3D::NextHit(int index, unsigned int view, bool inclDisabled) const
{
	pma::Hit3D* hit = 0;
	if (index < -1) index = -1;
	while (++index < (int)size()) // look for the next index of hit from the view
	{
		hit = fHits[index];
		if (hit->View2D() == view)
		{
			if (inclDisabled) break;
			else if (hit->IsEnabled()) break;
		}
	}
	return index;
}

int pma::Track3D::PrevHit(int index, unsigned int view, bool inclDisabled) const
{
	pma::Hit3D* hit = 0;
	if (index > (int)size()) index = (int)size();
	while (--index >= 0) // look for the prev index of hit from the view
	{
		hit = fHits[index];
		if (hit->View2D() == view)
		{
			if (inclDisabled) break;
			else if (hit->IsEnabled()) break;
		}
	}
	return index;
}

double pma::Track3D::HitDxByView(size_t index, unsigned int view,
	pma::Track3D::EDirection dir, bool secondDir) const
{
	pma::Hit3D* nexthit = 0;
	pma::Hit3D* hit = fHits[index];

	if (hit->View2D() != view)
	{
		mf::LogWarning("pma::Track3D") << "Function used with the hit not matching specified view.";
	}

	double dx = 0.0; // [cm]
	bool hitFound = false;
	int i = index;
	switch (dir)
	{
		case pma::Track3D::kForward:
			while (!hitFound && (++i < (int)size()))
			{
				nexthit = fHits[i];
				dx += sqrt(pma::Dist2(hit->Point3D(), nexthit->Point3D()));
				
				if (nexthit->View2D() == view) hitFound = true;
				else hitFound = false;

				hit = nexthit;
			}
			if (!hitFound)
			{
				if (!secondDir) dx = 0.5 * HitDxByView(index, view, pma::Track3D::kBackward, true);
				else { dx = Length(); mf::LogWarning("pma::Track3D") << "Single hit in this view."; }
			}
			break;

		case pma::Track3D::kBackward:
			while (!hitFound && (--i >= 0))
			{
				nexthit = fHits[i];
				dx += sqrt(pma::Dist2(hit->Point3D(), nexthit->Point3D()));

				if (nexthit->View2D() == view) hitFound = true;
				else hitFound = false;

				hit = nexthit;
			}
			if (!hitFound)
			{
				if (!secondDir) dx = 0.5 * HitDxByView(index, view, pma::Track3D::kForward, true);
				else { dx = Length(); mf::LogWarning("pma::Track3D") << "Single hit in this view."; }
			}
			break;

		default:
			mf::LogError("pma::Track3D") << "Direction undefined.";
			break;
	}
	return dx;
}

double pma::Track3D::HitDxByView(size_t index, unsigned int view) const
{
	if (index < size())
	{
		return 0.5 * (HitDxByView(index, view, pma::Track3D::kForward)
			+ HitDxByView(index, view, pma::Track3D::kBackward));
	}
	else
	{
		mf::LogError("pma::Track3D") << "Hit index out of range.";
		return 0.0;
	}
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

double pma::Track3D::GetRawdEdxSequence(
		std::map< size_t, std::vector<double> >& dedx,
		unsigned int view, unsigned int skip,
		bool inclDisabled) const
{
	dedx.clear();

	if (!size()) return 0.0;

	size_t step = 1;

	pma::Hit3D* hit = 0;

	double dr, dR, dq, dEq, qSkipped = 0.0;

	size_t j = NextHit(-1, view, inclDisabled), s = skip;
	if (j >= size()) return 0.0F; // no charged hits at all
	while (j < size()) // look for the first hit index
	{
		hit = fHits[j];
		dq = hit->SummedADC();
		if (s) { qSkipped += dq; s--; }
		else break;

		j = NextHit(j, view, inclDisabled);
	}

	size_t jmax = PrevHit(size(), view, inclDisabled);

	std::vector< size_t > indexes;
	TVector3 p0(0., 0., 0.), p1(0., 0., 0.);
	TVector2 c0(0., 0.), c1(0., 0.);
	while (j <= jmax)
	{
		indexes.clear(); // prepare to collect hit indexes for used for this dE/dx entry

		indexes.push_back(j);
		hit = fHits[j];

		p0 = hit->Point3D();
		p1 = hit->Point3D();

		c0.Set(hit->Wire(), hit->PeakTime());
		c1.Set(hit->Wire(), hit->PeakTime());

		dEq = hit->SummedADC(); // [now it is ADC sum]

		dr = HitDxByView(j, view, pma::Track3D::kForward); // protection against hits on the same position
		dR = HitDxByView(j, view); // dx seen by j-th hit

		size_t m = 1; // number of hits with charge > 0
		while (((m < step) || (dR < 0.1) || (dr == 0.0)) && (j <= jmax))
		{
			j = NextHit(j, view); // just next, even if tagged as outlier
			if (j > jmax) break; // no more hits in this view

			hit = fHits[j];
			if (!inclDisabled && !hit->IsEnabled())
			{
				if (dr == 0.0) continue;
				else break;
			}
			indexes.push_back(j);

			p1 = hit->Point3D();

			c1.Set(hit->Wire(), hit->PeakTime());

			dq = hit->SummedADC();

			dEq += dq;

			dr = HitDxByView(j, view, pma::Track3D::kForward);
			dR += HitDxByView(j, view);
			m++;
		}
		p0 += p1; p0 *= 0.5;
		c0 += c1; c0 *= 0.5;

		double range = Length(0, j);

		std::vector<double> trk_section;
		trk_section.push_back(c0.X());
		trk_section.push_back(c0.Y());
		trk_section.push_back(p0.X());
		trk_section.push_back(p0.Y());
		trk_section.push_back(p0.Z());
		trk_section.push_back(dEq);
		trk_section.push_back(dR);
		trk_section.push_back(range);

		for (auto const idx : indexes) dedx[idx] = trk_section;

		j = NextHit(j, view, inclDisabled);
	}

	return qSkipped;
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
	pma::Segment3D* maxSeg = 0;

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

		if (nHitsByView[midViewIdx] < 2) { mf::LogVerbatim("pma::Track3D") << "AddNode(): too few hits."; return false; }

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
			if (nHitsByView[1] && (mHits[1] < 2)) { mf::LogVerbatim("pma::Track3D") << "AddNode(): low Ind2 hits."; return false; }
			if (nHitsByView[2] && (mHits[2] < 2)) { mf::LogVerbatim("pma::Track3D") << "AddNode(): low Coll hits."; return false; }
		}

		maxSeg->SetProjection(maxSeg->Hit(i0));
		maxSeg->SetProjection(maxSeg->Hit(i1));

		unsigned int tpc = maxSeg->Hit(i0).TPC();
		unsigned int cryo = maxSeg->Hit(i0).Cryo();

		pma::Node3D* p = new pma::Node3D((maxSeg->Hit(i0).Point3D() + maxSeg->Hit(i1).Point3D()) * 0.5, tpc, cryo);

		//mf::LogVerbatim("pma::Track3D") << "add node x:" << p->Point3D().X()
		//	<< " y:" << p->Point3D().Y() << " z:" << p->Point3D().Z();
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

double pma::Track3D::GetMse(void) const
{
	double sumMse = 0.0;
	unsigned int nEnabledHits = 0;
	for (size_t i = 0; i < fNodes.size(); i++)
	{
		sumMse += fNodes[i]->SumDist2();
		nEnabledHits += fNodes[i]->NEnabledHits();
	}
	for (size_t i = 0; i < fSegments.size(); i++)
	{
		sumMse += fSegments[i]->SumDist2();
		nEnabledHits += fSegments[i]->NEnabledHits();
	}

	if (nEnabledHits) return sumMse / nEnabledHits;
	else return 0.0;
}

double pma::Track3D::GetObjFunction(float penaltyFactor) const
{
	double sum = 0.0;
	float p = penaltyFactor * fPenaltyValue;
	for (size_t i = 0; i < fNodes.size(); i++)
	{
		sum += fNodes[i]->GetObjFunction(p, fEndSegWeight);
	}
	return sum / fNodes.size();
}

double pma::Track3D::Optimize(int nNodes, double eps, bool selAllHits)
{
	if (!fNodes.size()) { mf::LogError("pma::Track3D") << "Track3D not initialized."; return 0.0; }

	UpdateParams();
	double g0 = GetObjFunction(), g1 = 0.0;
	if (g0 == 0.0) return g0;

	//mf::LogVerbatim("pma::Track3D") << "objective function at opt start: " << g0;

	bool stop = false;
	fMinSegStop = fSegments.size();
	fMaxSegStop = (int)(size() / fMaxSegStopFactor) + 1;
	do
	{
		bool stepDone = true;
		unsigned int stepIter = 0;
		do
		{
			double gstep = 1.0;
			unsigned int iter = 0;
			while ((gstep > eps) && (iter < 1000))
			{
				MakeProjection();
				UpdateParams();

				for (size_t j = 0; j < fNodes.size(); j++)
				{
					fNodes[j]->Optimize(fPenaltyValue, fEndSegWeight);
				}

				g1 = g0; g0 = GetObjFunction();

				//mf::LogVerbatim("pma::Track3D") << "obj fn: " << g0;
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
		switch (nNodes)
		{
			case 0: stop = true; break; // just optimize existing vertices

			case -1: // grow and optimize until automatic stop condition
				mf::LogVerbatim("pma::Track3D") << "optimized segments: " << fSegments.size();
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
				if (nNodes > 12)
				{
					if (AddNode()) { MakeProjection(); nNodes--; }
					else { mf::LogVerbatim("pma::Track3D") << "stop (3)"; stop = true; break; }

					if (AddNode())
					{
						MakeProjection(); nNodes--;
						if (AddNode()) nNodes--;
					}
				}
				else if (nNodes > 4)
				{
					if (AddNode()) { MakeProjection(); nNodes--; }
					else { mf::LogVerbatim("pma::Track3D") << "stop (2)"; stop = true; break; }

					if (AddNode()) nNodes--;
				}
				else
				{
					if (AddNode()) { nNodes--; }
					else { mf::LogVerbatim("pma::Track3D") << "stop (1)"; stop = true; break; }
				}
				break;
		}

	} while (!stop);
	//mf::LogVerbatim("pma::Track3D") << "Done (optimized segments: " << fSegments.size() << ").";

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

bool pma::Track3D::ShiftEndsToHits(void)
{
	pma::Element3D* el;
	pma::Node3D* vtx;

	if (!(fNodes.front()->Prev()))
	{
		el = GetNearestElement(front()->Point3D());
		vtx = dynamic_cast< pma::Node3D* >(el);
		if (vtx)
		{
			if (vtx == fNodes.front()) fNodes.front()->SetPoint3D(front()->Point3D());
			else
			{
				mf::LogWarning("pma::Track3D") << "First hit is projected to inner node.";
				return false;
			}
		}
		else
		{
			pma::Segment3D* seg = dynamic_cast< pma::Segment3D* >(el);
			if (seg)
			{
				if (seg->Prev() == fNodes.front())
				{
					double l0 = seg->Length();
					fNodes.front()->SetPoint3D(front()->Point3D());
					if ((seg->Length() < 0.2 * l0) && (fNodes.size() > 2))
					{
						mf::LogWarning("pma::Track3D") << "ShiftEndsToHits(): Short segment, node removed.";
						fNodes.erase(fNodes.begin() + 1);
					}
				}
				else
				{
					mf::LogWarning("pma::Track3D") << "First hit is projected to inner segment.";
					return false;
				}
			}
		}
	}

	if (!(fNodes.back()->NextCount()))
	{
		el = GetNearestElement(back()->Point3D());
		vtx = dynamic_cast< pma::Node3D* >(el);
		if (vtx)
		{
			if (vtx == fNodes.back()) fNodes.back()->SetPoint3D(back()->Point3D());
			else
			{
				mf::LogWarning("pma::Track3D") << "First hit is projected to inner node.";
				return false;
			}
		}
		else
		{
			pma::Segment3D* seg = dynamic_cast< pma::Segment3D* >(el);
			if (seg)
			{
				if (seg->Next() == fNodes.back())
				{
					double l0 = seg->Length();
					fNodes.back()->SetPoint3D(back()->Point3D());
					if ((seg->Length() < 0.2 * l0) && (fNodes.size() > 2))
					{
						mf::LogWarning("pma::Track3D") << "ShiftEndsToHits(): Short segment, node removed.";
						fNodes.erase(fNodes.end() - 1);
					}
				}
				else
				{
					mf::LogWarning("pma::Track3D") << "First hit is projected to inner segment.";
					return false;
				}
			}
		}
	}

	return true;
}

double pma::Track3D::Dist2(const TVector2& p2d, unsigned int view) const
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
	return min_dist;
}

double pma::Track3D::Dist2(const TVector3& p3d) const
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
	return min_dist;
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
		if (fSegments[i]->TPC() < 0) continue; // segment between TPC's

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
		for (size_t i = 0; i < size(); i++)
		{
			(*this)[i] = hits_tmp[i];
			/*mf::LogVerbatim("pma::Track3D")
				<< " w:" << hits_tmp[i]->Wire()
				<< " d:["
					<< hits_tmp[i]->Hit2DPtr()->StartTick()
					<< "; " << hits_tmp[i]->PeakTime()
					<< "; " << hits_tmp[i]->Hit2DPtr()->EndTick() << "]"
				<< " f:" << hits_tmp[i]->GetSegFraction()
				<< " p:" << hits_tmp[i]->View2D();*/
		}
	}
	else mf::LogError("pma::Track3D") << "Hit sorting problem.";
}

unsigned int pma::Track3D::DisableSingleViewEnds(void)
{
	SortHits();

	unsigned int nDisabled = 0;

	int hasHits[4];

	pma::Hit3D* nextHit = 0;
	int hitIndex = -1;

	bool stop = false;
	int nViews = 0;
	hasHits[0] = hasHits[1] = hasHits[2] = hasHits[3] = 0;
	do
	{
		pma::Node3D* vtx = fNodes.front();
		pma::Segment3D* seg = NextSegment(vtx);
		if (!seg) break;

		if (vtx->NPoints() + seg->NPoints() > 0) hasHits[3] = 1;

		for (size_t i = 0; i < vtx->NHits(); i++)
		{
			hitIndex = index_of(&(vtx->Hit(i)));
			if ((hitIndex >= 0) && (hitIndex + 1 < (int)size())) nextHit = fHits[hitIndex + 1];
			else nextHit = 0;

			if (vtx->Hit(i).IsEnabled()) hasHits[vtx->Hit(i).View2D()] = 1;
			if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
			nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
			if (nViews < 2)
			{
				if (vtx->Hit(i).IsEnabled())
				{
					vtx->Hit(i).SetEnabled(false);
					nDisabled++;
				}
			}
		}
		for (size_t i = 0; i < seg->NHits(); i++)
		{
			hitIndex = index_of(&(seg->Hit(i)));
			if ((hitIndex >= 0) && (hitIndex + 1 < (int)size())) nextHit = fHits[hitIndex + 1];
			else nextHit = 0;

			if (seg->Hit(i).IsEnabled()) hasHits[seg->Hit(i).View2D()] = 1;
			if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
			nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
			if (nViews < 2)
			{
				if (seg->Hit(i).IsEnabled())
				{
					seg->Hit(i).SetEnabled(false);
					nDisabled++;
				}
			}
		}

		if (fNodes.size() < 3) break;

		nViews = hasHits[1] + hasHits[2] + hasHits[3];
		if (hasHits[0] || (nViews > 1)) stop = true;
		else
		{
			pma::Node3D* vtx_front = fNodes.front();
			fNodes.erase(fNodes.begin());
			delete vtx_front;
		}

	} while (!stop);

	stop = false;
	nViews = 0;
	hasHits[0] = hasHits[1] = hasHits[2] = hasHits[3] = 0;
	do
	{
		pma::Node3D* vtx = fNodes.back();
		pma::Segment3D* seg = PrevSegment(vtx);
		if (!seg) break;

		if (vtx->NPoints() || seg->NPoints()) hasHits[3] = 1;

		for (int i = vtx->NHits() - 1; i >= 0; i--)
		{
			hitIndex = index_of(&(vtx->Hit(i)));
			if ((hitIndex >= 0) && (hitIndex - 1 >= 0)) nextHit = fHits[hitIndex - 1];
			else nextHit = 0;

			if (vtx->Hit(i).IsEnabled()) hasHits[vtx->Hit(i).View2D()] = 1;
			if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
			nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
			if (nViews < 2)
			{
				if (vtx->Hit(i).IsEnabled())
				{
					vtx->Hit(i).SetEnabled(false);
					nDisabled++;
				}
			}
		}
		for (int i = seg->NHits() - 1; i >= 0; i--)
		{
			hitIndex = index_of(&(seg->Hit(i)));
			if ((hitIndex >= 0) && (hitIndex - 1 >= 0)) nextHit = fHits[hitIndex - 1];
			else nextHit = 0;

			if (seg->Hit(i).IsEnabled()) hasHits[seg->Hit(i).View2D()] = 1;
			if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
			nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
			if (nViews < 2)
			{
				if (seg->Hit(i).IsEnabled())
				{
					seg->Hit(i).SetEnabled(false);
					nDisabled++;
				}
			}
		}

		if (fNodes.size() < 3) break;

		nViews = hasHits[1] + hasHits[2] + hasHits[3];
		if (hasHits[0] || (nViews > 1)) stop = true;
		else
		{
			pma::Node3D* vtx_back = fNodes.back();
			fNodes.pop_back();
			delete vtx_back;
		}

	} while (!stop);

	RebuildSegments();
	MakeProjection();

	return nDisabled;
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
		pma::Hit3D* hit = fHits[i];
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

	if (fraction < 1.0F)
	{
		pma::Element3D* pe = FirstElement();
		for (unsigned int i = 0; i < pe->NHits(); i++)
			pe->Hit(i).SetEnabled(true);

		pe = LastElement();
		for (unsigned int i = 0; i < pe->NHits(); i++)
			pe->Hit(i).SetEnabled(true);
	}
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

