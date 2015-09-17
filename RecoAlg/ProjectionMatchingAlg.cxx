////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ProjectionMatchingAlg.h"

#include "RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::ProjectionMatchingAlg::ProjectionMatchingAlg(const fhicl::ParameterSet& pset)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

pma::ProjectionMatchingAlg::~ProjectionMatchingAlg(void)
{
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fOptimizationEps = p.get< double >("OptimizationEps");
	fFineTuningEps = p.get< double >("FineTuningEps");

	fTrkValidationDist2D = p.get< double >("TrkValidationDist2D");
	fHitTestingDist2D = p.get< double >("HitTestingDist2D");

	fMinTwoViewFraction = p.get< double >("MinTwoViewFraction");

	pma::Node3D::SetMargin(p.get< double >("NodeMargin3D"));

	pma::Element3D::SetOptFactor(geo::kZ, p.get< double >("HitWeightZ"));
	pma::Element3D::SetOptFactor(geo::kV, p.get< double >("HitWeightV"));
	pma::Element3D::SetOptFactor(geo::kU, p.get< double >("HitWeightU"));
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView) const
{
	double step = 0.3;
	double max_d = fTrkValidationDist2D;
	double d2, max_d2 = max_d * max_d;
	unsigned int nAll = 0, nPassed = 0;

	std::vector< unsigned int > trkTPCs = trk.TPCs();
	std::vector< unsigned int > trkCryos = trk.Cryos();
	std::map< std::pair< unsigned int, unsigned int >, std::pair< TVector2, TVector2 > > ranges;
	for (auto c : trkCryos)
		for (auto t : trkTPCs)
		{
			ranges[std::pair< unsigned int, unsigned int >(t, c)] = trk.WireDriftRange(testView, t, c);
		}

	unsigned int tpc, cryo;
	std::map< std::pair< unsigned int, unsigned int >, std::vector< TVector2 > > all_close_points;

	for (const auto h : hits)
		if (h->WireID().Plane == testView)
	{
		std::pair< unsigned int, unsigned int > tpc_cryo(h->WireID().TPC, h->WireID().Cryostat);
		std::pair< TVector2, TVector2 > rect = ranges[tpc_cryo];

		if ((h->WireID().Wire > rect.first.X() - 10) &&  // chceck only hits in the rectangle around
		    (h->WireID().Wire < rect.second.X() + 10) && // the track projection, it is faster than
		    (h->PeakTime() > rect.first.Y() - 100) &&    // calculation of trk.Dist2(p2d, testView)
		    (h->PeakTime() < rect.second.Y() + 100))
		{
			TVector2 p2d = pma::WireDriftToCm(h->WireID().Wire, h->PeakTime(), testView, tpc_cryo.first, tpc_cryo.second);

			d2 = trk.Dist2(p2d, testView);

			if (d2 < max_d2) all_close_points[tpc_cryo].push_back(p2d);
		}
	}

	// then check how points close to the track projection are distributed along the
	// track, namely: are there track sections crossing empty spaces?
	TVector3 p(trk.front()->Point3D());
	for (size_t i = 0; i < trk.Nodes().size() - 1; i++)
	{
		tpc = trk.Nodes()[i]->TPC();
		cryo = trk.Nodes()[i]->Cryo();

		TVector3 vNext(trk.Nodes()[i + 1]->Point3D());
		TVector3 vThis(trk.Nodes()[i]->Point3D());

		const std::vector< TVector2 >& points = all_close_points[std::pair< unsigned int, unsigned int >(tpc, cryo)];
		if (trk.Nodes()[i + 1]->TPC() == (int)tpc) // skip segments between tpc's
		{
			TVector3 dc(vNext); dc -= vThis;
			dc *= step / dc.Mag();

			double f = pma::GetSegmentProjVector(p, vThis, vNext);
			while ((f < 1.0) && trk.Nodes()[i]->SameTPC(p))
			{
				if (points.size())
				{
					TVector2 p2d = pma::GetProjectionToPlane(p, testView, tpc, cryo);

					for (const auto & h : points)
					{
						d2 = pma::Dist2(p2d, h);
						if (d2 < max_d2) { nPassed++; break; }
					}
				}
				nAll++;

				p += dc; f = pma::GetSegmentProjVector(p, vThis, vNext);
			}
		}
		p = vNext;
	}

	if (nAll > 3) // validate actually only if 2D projection in testView has some minimum length
	{
		double v = nPassed / (double)nAll;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  trk fraction ok: " << v;
		return v;
	}
	else return 1.0;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(
	const TVector3& p0, const TVector3& p1,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, unsigned int tpc, unsigned int cryo) const
{
	double step = 0.3;
	double max_d = fTrkValidationDist2D;
	double d2, max_d2 = max_d * max_d;
	unsigned int nAll = 0, nPassed = 0;

	TVector3 p(p0);
	TVector3 dc(p1); dc -= p;
	dc *= step / dc.Mag();

	double f = pma::GetSegmentProjVector(p, p0, p1);
	while (f < 1.0)
	{
		TVector2 p2d = pma::GetProjectionToPlane(p, testView, tpc, cryo);

		for (const auto & h : hits)
			if (h->WireID().Plane == testView)
		{
			d2 = pma::Dist2(p2d, pma::WireDriftToCm(h->WireID().Wire, h->PeakTime(), testView, tpc, cryo));
			if (d2 < max_d2) { nPassed++; break; }
		}
		nAll++;

		p += dc; f = pma::GetSegmentProjVector(p, p0, p1);
	}

	if (nAll > 3) // validate actually only if 2D projection in testView has some minimum length
	{
		double v = nPassed / (double)nAll;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  segment fraction ok: " << v;
		return v;
	}
	else return 1.0;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::twoViewFraction(pma::Track3D& trk) const
{
	trk.SelectHits();
	trk.DisableSingleViewEnds();

	size_t idx = 0;
	while ((idx < trk.size() - 1) && !trk[idx]->IsEnabled()) idx++;
	double l0 = trk.Length(0, idx + 1);

	idx = trk.size() - 1;
	while ((idx > 1) && !trk[idx]->IsEnabled()) idx--;
	double l1 = trk.Length(idx - 1, trk.size() - 1);

	trk.SelectHits();

	return 1.0 - (l0 + l1) / trk.Length();
}
// ------------------------------------------------------

size_t pma::ProjectionMatchingAlg::getSegCount(size_t trk_size)
{
	int nSegments = (int)( 0.8 * trk_size / sqrt(trk_size) );

	if (nSegments > 1) return (size_t)nSegments;
	else return 1;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildTrack(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2) const
{
	pma::Track3D* trk = new pma::Track3D(); // track candidate
	trk->AddHits(hits_1);
	trk->AddHits(hits_2);

	mf::LogVerbatim("ProjectionMatchingAlg") << "track size: " << trk->size();
	std::vector< unsigned int > tpcs = trk->TPCs();
	for (size_t t = 0; t < tpcs.size(); ++t)
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "  tpc:" << tpcs[t];
	}
	mf::LogVerbatim("ProjectionMatchingAlg")
		<< "  #coll:" << trk->NHits(geo::kZ)
		<< "  #ind2:" << trk->NHits(geo::kV)
		<< "  #ind1:" << trk->NHits(geo::kU);

	size_t nSegments = getSegCount(trk->size());
	size_t nNodes = (size_t)( nSegments - 1 ); // n nodes to add

	mf::LogVerbatim("ProjectionMatchingAlg") << "  initialize trk";
	trk->Initialize();

	double f = twoViewFraction(*trk);
	if (f > fMinTwoViewFraction)
	{
		double g = 0.0;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  optimize trk (" << nSegments << " seg)";
		if (nNodes)
		{
			g = trk->Optimize(nNodes, fOptimizationEps);   // build nodes
			mf::LogVerbatim("ProjectionMatchingAlg") << "  nodes done, g = " << g;
		}
		g = trk->Optimize(0, fFineTuningEps);              // final tuning
		mf::LogVerbatim("ProjectionMatchingAlg") << "  tune done, g = " << g;

		trk->SortHits();
		// trk->ShiftEndsToHits(); // not sure if useful already here
		return trk;
	}
	else
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "  clusters do not match, f = " << f;
		delete trk;
		return 0;
	}
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2) const
{
	pma::Track3D* trk = new pma::Track3D();
	trk->SetEndSegWeight(0.001F);
	trk->AddHits(hits_1);
	trk->AddHits(hits_2);

	trk->Initialize(0.001F);
	trk->Optimize(0, fFineTuningEps);

	trk->SortHits();
	return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2,
	const TVector3& point) const
{
	pma::Track3D* trk = buildSegment(hits_1, hits_2);

	double dfront = pma::Dist2(trk->front()->Point3D(), point);
	double dback = pma::Dist2(trk->back()->Point3D(), point);
	if (dfront > dback) trk->Flip();

	trk->Nodes().front()->SetPoint3D(point);
	trk->Nodes().front()->SetFrozen(true);
	trk->Optimize(0, fFineTuningEps);

	trk->SortHits();
	return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	const TVector3& point) const
{
	pma::Track3D* trk = buildSegment(hits);

	double dfront = pma::Dist2(trk->front()->Point3D(), point);
	double dback = pma::Dist2(trk->back()->Point3D(), point);
	if (dfront > dback) trk->Flip();

	trk->Nodes().front()->SetPoint3D(point);
	trk->Nodes().front()->SetFrozen(true);
	trk->Optimize(0, fFineTuningEps);

	trk->SortHits();
	return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::extendTrack(
	const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	bool add_nodes) const
{
	pma::Track3D* copy = new pma::Track3D(trk);
	copy->AddHits(hits);

	mf::LogVerbatim("ProjectionMatchingAlg") << "ext. track size: " << copy->size()
		<< "  #coll:" << copy->NHits(geo::kZ)
		<< "  #ind2:" << copy->NHits(geo::kV)
		<< "  #ind1:" << copy->NHits(geo::kU);

	if (add_nodes)
	{
		size_t nSegments = getSegCount(copy->size());
		int nNodes = nSegments - copy->Nodes().size() + 1; // n nodes to add
		if (nNodes < 0) nNodes = 0;

		if (nNodes)
		{
			mf::LogVerbatim("ProjectionMatchingAlg") << "  add " << nNodes << " nodes";
			copy->Optimize(nNodes, fOptimizationEps);
		}
	}
	double g = copy->Optimize(0, fFineTuningEps);
	mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt done, g = " << g;

	return copy;
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::fixReadoutPlanesParallel(
		pma::Track3D& trk, const std::vector< art::Ptr<recob::Hit> >& hits) const
{
}
// ------------------------------------------------------

std::vector< pma::Hit3D* > pma::ProjectionMatchingAlg::trimTrackToVolume(
	pma::Track3D& trk, TVector3 p0, TVector3 p1) const
{
	std::vector< pma::Hit3D* > trimmedHits;



	return trimmedHits;
}
// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::alignTracks(pma::Track3D& first, pma::Track3D& second) const
{
	unsigned int k = 0;
	double distFF = pma::Dist2(first.front()->Point3D(), second.front()->Point3D());
	double dist = distFF;

	double distFB = pma::Dist2(first.front()->Point3D(), second.back()->Point3D());
	if (distFB < dist) { k = 1; dist = distFB; }

	double distBF = pma::Dist2(first.back()->Point3D(), second.front()->Point3D());
	if (distBF < dist) { k = 2; dist = distBF; }

	double distBB = pma::Dist2(first.back()->Point3D(), second.back()->Point3D());
	if (distBB < dist) { k = 3; dist = distBB; }

	switch (k) // flip to get dst end before src start, do not merge if track's order reversed
	{
		case 0:	first.Flip(); break;
		case 1: mf::LogError("PMAlgTrackMaker") << "Tracks in reversed order."; return false;
		case 2: break;
		case 3: second.Flip(); break;
		default: mf::LogError("PMAlgTrackMaker") << "Should never happen."; return false;
	}
	return true;
}

void pma::ProjectionMatchingAlg::mergeTracks(pma::Track3D& dst, pma::Track3D& src, bool reopt) const
{
	if (!alignTracks(dst, src)) return;

	unsigned int tpc = src.FrontTPC();
	unsigned int cryo = src.FrontCryo();
	double lmean = dst.Length() / (dst.Nodes().size() - 1);
	if ((pma::Dist2(
			dst.Nodes().back()->Point3D(),
			src.Nodes().front()->Point3D()) > 0.5 * lmean) ||
	    (tpc != dst.BackTPC()) || (cryo != dst.BackCryo()))
	{
		dst.AddNode(src.Nodes().front()->Point3D(), tpc, cryo);
		if (src.Nodes().front()->IsFrozen())
			dst.Nodes().back()->SetFrozen(true);
	}
	for (size_t n = 1; n < src.Nodes().size(); n++)
	{
		pma::Node3D* node = src.Nodes()[n];

		dst.AddNode(src.Nodes()[n]->Point3D(),
		            node->TPC(), node->Cryo());

		if (node->IsFrozen())
			dst.Nodes().back()->SetFrozen(true);
	}
	for (size_t h = 0; h < src.size(); h++)
	{
		dst.push_back(src[h]->Hit2DPtr());
	}
	if (reopt)
	{
		double g = dst.Optimize(0, fFineTuningEps);
		mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt after merging done, g = " << g;
	}
	else
	{
		dst.MakeProjection();
	}

	dst.SortHits();
	dst.ShiftEndsToHits();
}

void pma::ProjectionMatchingAlg::autoFlip(pma::Track3D& trk,
	pma::Track3D::EDirection dir, double thr, unsigned int n) const
{
	unsigned int nViews = 3;
	std::map< size_t, std::vector<double> > dedx_map[3];
	for (unsigned int i = 0; i < nViews; i++)
	{
		trk.GetRawdEdxSequence(dedx_map[i], i, 1);
	}
	unsigned int bestView = 2;
	if (dedx_map[0].size() > 2 * dedx_map[2].size()) bestView = 0;
	if (dedx_map[1].size() > 2 * dedx_map[2].size()) bestView = 1;

	std::vector< std::vector<double> > dedx;
	for (size_t i = 0; i < trk.size(); i++)
	{
		auto it = dedx_map[bestView].find(i);
		if (it != dedx_map[bestView].end())
		{
			dedx.push_back(it->second);
		}
	}

	float dEdxStart = 0.0F, dEdxStop = 0.0F;
	float dEStart = 0.0F, dxStart = 0.0F;
	float dEStop = 0.0F, dxStop = 0.0F;
	if (dedx.size() > 4)
	{
		if (!n) // use default options
		{
			if (dedx.size() > 30) n = 12;
			else if (dedx.size() > 20) n = 8;
			else if (dedx.size() > 10) n = 4;
			else n = 3;
		}

		size_t k = (dedx.size() - 2) >> 1;
		if (n > k) n = k;

		for (size_t i = 1, j = 0; j < n; i++, j++)
		{
			dEStart += dedx[i][5]; dxStart += dedx[i][6];
		}
		if (dxStart > 0.0F) dEdxStart = dEStart / dxStart;

		for (size_t i = dedx.size() - 2, j = 0; j < n; i--, j++)
		{
			dEStop += dedx[i][5]; dxStop += dedx[i][6];
		}
		if (dxStop > 0.0F) dEdxStop = dEStop / dxStop;
	}
	else if (dedx.size() == 4)
	{
		dEStart = dedx[0][5] + dedx[1][5]; dxStart = dedx[0][6] + dedx[1][6];
		dEStop = dedx[2][5] + dedx[3][5]; dxStop = dedx[2][6] + dedx[3][6];
		if (dxStart > 0.0F) dEdxStart = dEStart / dxStart;
		if (dxStop > 0.0F) dEdxStop = dEStop / dxStop;

	}
	else if (dedx.size() > 1)
	{
		if (dedx.front()[2] > 0.0F) dEdxStart = dedx.front()[5] / dedx.front()[6];
		if (dedx.back()[2] > 0.0F) dEdxStop = dedx.back()[5] / dedx.back()[6];
	}
	else return;

	if ((dir == pma::Track3D::kForward) && ((1.0 + thr) * dEdxStop < dEdxStart)) trk.Flip();  // particle stops at the end of the track
	if ((dir == pma::Track3D::kBackward) && (dEdxStop > (1.0 + thr) * dEdxStart)) trk.Flip(); // particle stops at the front of the track
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::selectInitialHits(pma::Track3D& trk, unsigned int view) const
{
	for (size_t i = 0; i < trk.size(); i++)
	{
		pma::Hit3D* hit = trk[i];
		if (hit->View2D() == view)
		{
			if ((hit->GetDistToProj() > 0.5) || // more than 0.5cm away away from the segment
			    (hit->GetSegFraction() < -1.0)) // projects before segment start (to check!!!)
				hit->TagOutlier(true);
			else hit->TagOutlier(false);
		}
	}

	unsigned int nhits = 0;
	double last_x, dx = 0.0, last_q, dq = 0.0, dqdx = 0.0;
	int ih = trk.NextHit(-1, view);

	pma::Hit3D* hit = trk[ih];
	pma::Hit3D* lastHit = hit;

	if ((ih >= 0) && (ih < (int)trk.size()))
	{
		hit->TagOutlier(true);

		ih = trk.NextHit(ih, view);
		while ((dx < 2.5) && (ih >= 0) && (ih < (int)trk.size()))
		{
			hit = trk[ih];

			if (abs(hit->Wire() - lastHit->Wire()) > 2)
				break; // break on gap in wire direction

			last_x = trk.HitDxByView(ih, view);
			last_q = hit->SummedADC();
			if (dx + last_x < 3.0)
			{
				dq += last_q;
				dx += last_x;
				nhits++;
			}
			else break;

			lastHit = hit;
			ih = trk.NextHit(ih, view);
		}
		while ((ih >= 0) && (ih < (int)trk.size()))
		{
			hit = trk[ih];
			hit->TagOutlier(true);
			ih = trk.NextHit(ih, view);
		}
	}
	else { mf::LogError("ProjectionMatchingAlg") << "Initial part selection failed."; }

	if (!nhits) { mf::LogError("ProjectionMatchingAlg") << "Initial part too short to select useful hits."; }

	if (dx > 0.0) dqdx = dq / dx;

	//std::cout << "nhits=" << nhits << ", dq=" << dq << ", dx=" << dx << std::endl;
	return dqdx;
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::setTrackTag(pma::Track3D& trk) const
{
	double length = trk.Length();
	double mse = trk.GetMse();
	double meanAngle = trk.GetMeanAng();

	if ( // (length < 80.0) &&  // tag only short tracks as EM shower-like
	    ((mse > 0.0001 * length) || (meanAngle < 3.0)))
	{
		trk.SetTag(pma::Track3D::kEmLike);
	}
	else
	{
		trk.SetTag(pma::Track3D::kTrackLike);
	}
}
// ------------------------------------------------------


