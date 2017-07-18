////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTracking
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl),
//              R.Sulej (Robert.Sulej@cern.ch),
//              L.Whitehead (leigh.howard.whitehead@cern.ch), June 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgTracking.h"
#include "larreco/RecoAlg/PMAlgStitching.h"

#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

using Point_t      = recob::tracking::Point_t;
using Vector_t     = recob::tracking::Vector_t;
using SMatrixSym55 = recob::tracking::SMatrixSym55;

recob::Track pma::convertFrom(const pma::Track3D& src, unsigned int tidx, int pdg)
{
    std::vector<Point_t> positions;                    positions.reserve(src.size());
    std::vector<Vector_t> momenta;                     momenta.reserve(src.size());
    std::vector<recob::TrajectoryPointFlags> outFlags; outFlags.reserve(src.size());

	for (size_t i = 0, h = 0; i < src.size(); i++)
		if (src[i]->IsEnabled())
	{
        auto const & point3d = src[i]->Point3D();
        positions.emplace_back(point3d.X(), point3d.Y(), point3d.Z());
        momenta.push_back(src.GetDirection3D(i));
        outFlags.emplace_back(h++, recob::TrajectoryPointFlags::makeMask());
	}

    int ndof = 0;
    float totChi2 = 0;

    SMatrixSym55 covStart, covEnd;
    return recob::Track(
        recob::TrackTrajectory(std::move(positions), std::move(momenta), std::move(outFlags), false),
        pdg, totChi2, ndof, covStart, covEnd, tidx);
}
// ------------------------------------------------------

pma::PMAlgTrackingBase::PMAlgTrackingBase(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
                                          const pma::ProjectionMatchingAlg::Config& pmalgConfig,
                                          const pma::PMAlgVertexing::Config& pmvtxConfig) :
	fProjectionMatchingAlg(pmalgConfig),
	fPMAlgVertexing(pmvtxConfig)
{
	unsigned int cryo, tpc, view;
	for (auto const& h : allhitlist)
	{
		cryo = h->WireID().Cryostat;
		tpc = h->WireID().TPC;
		view = h->WireID().Plane;

		fHitMap[cryo][tpc][view].push_back(h);
	}
}
// ------------------------------------------------------

pma::PMAlgTrackingBase::~PMAlgTrackingBase(void)
{
	for (auto t : fResult.tracks()) t.DeleteTrack();
}
// ------------------------------------------------------

void pma::PMAlgTrackingBase::guideEndpoints(pma::TrkCandidateColl & tracks)
{
	for (auto const & t : tracks.tracks())
	{
		auto & trk = *(t.Track());

		unsigned int tpc = trk.FrontTPC(), cryo = trk.FrontCryo();
		if ((tpc == trk.BackTPC()) && (cryo == trk.BackCryo()))
		{
			fProjectionMatchingAlg.guideEndpoints(trk, fHitMap[cryo][tpc]);
		}
		else
		{
			fProjectionMatchingAlg.guideEndpoints(trk, pma::Track3D::kBegin,
				fHitMap[trk.FrontCryo()][trk.FrontTPC()]);
			fProjectionMatchingAlg.guideEndpoints(trk, pma::Track3D::kEnd,
				fHitMap[trk.BackCryo()][trk.BackTPC()]);
		}
	}
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

pma::PMAlgFitter::PMAlgFitter(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
                              const std::vector< recob::Cluster > & clusters,
                              const std::vector< recob::PFParticle > & pfparticles,
                              const art::FindManyP< recob::Hit > & hitsFromClusters,
                              const art::FindManyP< recob::Cluster > & clusFromPfps,
                              const art::FindManyP< recob::Vertex > & vtxFromPfps,
                              const pma::ProjectionMatchingAlg::Config& pmalgConfig,
                              const pma::PMAlgFitter::Config& pmalgFitterConfig,
                              const pma::PMAlgVertexing::Config& pmvtxConfig) :
	PMAlgTrackingBase(allhitlist, pmalgConfig, pmvtxConfig),
	fTrackingOnlyPdg(pmalgFitterConfig.TrackingOnlyPdg()),
	fTrackingSkipPdg(pmalgFitterConfig.TrackingSkipPdg()),
	fRunVertexing(pmalgFitterConfig.RunVertexing())
{
	mf::LogVerbatim("PMAlgFitter") << "Found " << allhitlist.size() << "hits in the event.";
	mf::LogVerbatim("PMAlgFitter") << "Sort hits by clusters assigned to PFParticles...";

	fCluHits.resize(clusters.size());
	for (size_t i = 0; i < pfparticles.size(); ++i)
	{
		fPfpPdgCodes[i] = pfparticles[i].PdgCode();

		auto cv = clusFromPfps.at(i);
		for (const auto & c : cv)
		{
			fPfpClusters[i].push_back(c);
			if (fCluHits[c.key()].empty())
			{
				auto hv = hitsFromClusters.at(c.key());
				fCluHits[c.key()].reserve(hv.size());
				for (auto const & h : hv) fCluHits[c.key()].push_back(h);
			}
		}

		if (vtxFromPfps.isValid() && vtxFromPfps.at(i).size())
		{
			double xyz[3];
			vtxFromPfps.at(i).front()->XYZ(xyz);
			fPfpVtx[i] = pma::Vector3D(xyz[0], xyz[1], xyz[2]);
		}
	}

	mf::LogVerbatim("PMAlgFitter") << "...done, "
		<< fCluHits.size() << " clusters from "
		<< fPfpClusters.size() << " pfparticles for 3D tracking.";
}
// ------------------------------------------------------

// ------------------------------------------------------
// ------------------------------------------------------
int pma::PMAlgFitter::build(void)
{
    if (!fPfpClusters.empty() && !fCluHits.empty())
    {
			// build pm tracks
			buildTracks();

			// add 3D ref.points for clean endpoints of wire-plae parallel tracks
			guideEndpoints(fResult);

			if (fRunVertexing) fPMAlgVertexing.run(fResult);

			// build segment of shower
			buildShowers();
    }
    else
    {
        mf::LogWarning("PMAlgFitter") << "no clusters, no pfparticles";
        return -1;
    }
    
    return fResult.size();
}
// ------------------------------------------------------
// ------------------------------------------------------

void pma::PMAlgFitter::buildTracks(void)
{
		bool skipPdg = true;
		if (!fTrackingSkipPdg.empty() && (fTrackingSkipPdg.front() == 0)) skipPdg = false;

		bool selectPdg = true;
		if (!fTrackingOnlyPdg.empty() && (fTrackingOnlyPdg.front() == 0)) selectPdg = false;

		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (pdg == 11) continue;
			if (skipPdg && has(fTrackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(fTrackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgFitter") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

			std::vector< art::Ptr<recob::Hit> > allHits;

			pma::TrkCandidate candidate;
			for (const auto & c : pfpCluEntry.second)
			{
				candidate.Clusters().push_back(c.key());

				allHits.reserve(allHits.size() + fCluHits.at(c.key()).size());
				for (const auto & h : fCluHits.at(c.key()))
				{
					allHits.push_back(h);
				}
			}
			candidate.SetKey(pfpCluEntry.first);

			candidate.SetTrack(fProjectionMatchingAlg.buildMultiTPCTrack(allHits));

			if (candidate.IsValid() &&
			    candidate.Track()->HasTwoViews() &&
			    (candidate.Track()->Nodes().size() > 1))
			{
			    if (!std::isnan(candidate.Track()->Length())) { fResult.push_back(candidate); }
			    else
			    {
			        mf::LogError("PMAlgFitter") << "Trajectory fit lenght is nan.";
				    candidate.DeleteTrack();
			    }
			}
			else
			{
				candidate.DeleteTrack();
			}
		}
}
// ------------------------------------------------------

void pma::PMAlgFitter::buildShowers(void)
{
		bool skipPdg = true;
		if (!fTrackingSkipPdg.empty() && (fTrackingSkipPdg.front() == 0))
			skipPdg = false;

		bool selectPdg = true;
		if (!fTrackingOnlyPdg.empty() && (fTrackingOnlyPdg.front() == 0))
			selectPdg = false;

		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (pdg != 11) continue;
			if (skipPdg && has(fTrackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(fTrackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgFitter") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

			std::vector< art::Ptr<recob::Hit> > allHits;

			pma::TrkCandidate candidate;
			for (const auto & c : pfpCluEntry.second)
			{
				candidate.Clusters().push_back(c.key());

				allHits.reserve(allHits.size() + fCluHits.at(c.key()).size());
				for (const auto & h : fCluHits.at(c.key()))
					allHits.push_back(h);
			}

			candidate.SetKey(pfpCluEntry.first);

			mf::LogVerbatim("PMAlgFitter") << "building..." << ", pdg:" << pdg;

			auto search = fPfpVtx.find(pfPartIdx);
			if (search != fPfpVtx.end())
			{
				candidate.SetTrack(fProjectionMatchingAlg.buildShowerSeg(allHits, search->second));

				if (candidate.IsValid()
						&& candidate.Track()->HasTwoViews()
						&& (candidate.Track()->Nodes().size() > 1)
						&& !std::isnan(candidate.Track()->Length()))
				{
					fResult.push_back(candidate);
				}
				else
				{
					candidate.DeleteTrack();
				}
			}
		}
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

pma::PMAlgTracker::PMAlgTracker(const std::vector< art::Ptr<recob::Hit> > & allhitlist, const std::vector<recob::Wire> & wires,
	const pma::ProjectionMatchingAlg::Config& pmalgConfig,
	const pma::PMAlgTracker::Config& pmalgTrackerConfig,
	const pma::PMAlgVertexing::Config& pmvtxConfig,
	const pma::PMAlgStitching::Config& pmstitchConfig,
	const pma::PMAlgCosmicTagger::Config& pmtaggerConfig,
	
	const std::vector< TH1F* > & hclose, const std::vector< TH1F* > & hdist) :

	PMAlgTrackingBase(allhitlist, pmalgConfig, pmvtxConfig),

    fWires(wires),

	fMinSeedSize1stPass(pmalgTrackerConfig.MinSeedSize1stPass()),
	fMinSeedSize2ndPass(pmalgTrackerConfig.MinSeedSize2ndPass()),
	fTrackLikeThreshold(pmalgTrackerConfig.TrackLikeThreshold()),

	fMinTwoViewFraction(pmalgConfig.MinTwoViewFraction()),

	fFlipToBeam(pmalgTrackerConfig.FlipToBeam()),
	fFlipDownward(pmalgTrackerConfig.FlipDownward()),
	fAutoFlip_dQdx(pmalgTrackerConfig.AutoFlip_dQdx()),

	fMergeWithinTPC(pmalgTrackerConfig.MergeWithinTPC()),
	fMergeTransverseShift(pmalgTrackerConfig.MergeTransverseShift()),
	fMergeAngle(pmalgTrackerConfig.MergeAngle()),

    fCosmicTagger(pmtaggerConfig),
    fTagCosmicTracks(fCosmicTagger.tagAny()),

	fStitchBetweenTPCs(pmalgTrackerConfig.StitchBetweenTPCs()),
	fStitchDistToWall(pmalgTrackerConfig.StitchDistToWall()),
	fStitchTransverseShift(pmalgTrackerConfig.StitchTransverseShift()),
	fStitchAngle(pmalgTrackerConfig.StitchAngle()),

	fMatchT0inAPACrossing(pmalgTrackerConfig.MatchT0inAPACrossing()),
	fMatchT0inCPACrossing(pmalgTrackerConfig.MatchT0inCPACrossing()),
    fStitcher(pmstitchConfig),

	fRunVertexing(pmalgTrackerConfig.RunVertexing()),

	fAdcInPassingPoints(hclose), fAdcInRejectedPoints(hdist),

    fGeom(&*(art::ServiceHandle<geo::Geometry>())),
	fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    for (const auto v : fGeom->Views()) { fAvailableViews.push_back(v); }
    std::reverse(fAvailableViews.begin(), fAvailableViews.end());

    mf::LogVerbatim("PMAlgTracker") << "Using views in the following order:";
    for (const auto v : fAvailableViews) {  mf::LogInfo("PMAlgTracker") << " " << v; }

    // track validation settings:
    mf::LogVerbatim("PMAlgTracker") << "Validation mode: " << pmalgTrackerConfig.Validation();

    size_t nplanes = fGeom->MaxPlanes();
    for (size_t p = 0; p < nplanes; ++p) { fAdcImages.emplace_back(pmalgTrackerConfig.AdcImageAlg()); }

    if (pmalgTrackerConfig.Validation() == "hits")         { fValidation = pma::PMAlgTracker::kHits;  }
    else if (pmalgTrackerConfig.Validation() == "adc")     { fValidation = pma::PMAlgTracker::kAdc;   }
    else if (pmalgTrackerConfig.Validation() == "calib")   { fValidation = pma::PMAlgTracker::kCalib; }
    else { throw cet::exception("pma::PMAlgTracker") << "validation name not supported" << std::endl; }

    fAdcValidationThr = pmalgTrackerConfig.AdcValidationThr();
    if (fValidation == pma::PMAlgTracker::kAdc) { mf::LogVerbatim("PMAlgTracker") << "Validation ADC threshold: " << fAdcValidationThr; }
}
// ------------------------------------------------------

void pma::PMAlgTracker::init(const art::FindManyP< recob::Hit > & hitsFromClusters)
{
	mf::LogVerbatim("PMAlgTracker") << "Sort hits by clusters...";
	fCluHits.clear(); fCluHits.reserve(hitsFromClusters.size());
	fCluWeights.clear(); fCluWeights.reserve(fCluHits.size());
	for (size_t i = 0; i < hitsFromClusters.size(); ++i)
	{
		auto v = hitsFromClusters.at(i);
		fCluHits.push_back(std::vector< art::Ptr<recob::Hit> >());
		for (auto const & h : v) fCluHits.back().push_back(h);
		fCluWeights.push_back(1);
	}
	mf::LogVerbatim("PMAlgTracker") << "...done, " << fCluHits.size() << " clusters for 3D tracking.";
}
// ------------------------------------------------------

void pma::PMAlgTracker::init(const art::FindManyP< recob::Hit > & hitsFromClusters,
                             const std::vector< float > & trackLike)
{
	mf::LogVerbatim("PMAlgTracker") << "Filter track-like clusters using likelihood values...";
	fCluHits.clear(); fCluHits.reserve(hitsFromClusters.size());
	fCluWeights.clear(); fCluWeights.reserve(fCluHits.size());
	for (size_t i = 0; i < hitsFromClusters.size(); ++i)
	{
		auto v = hitsFromClusters.at(i);
		fCluHits.push_back(std::vector< art::Ptr<recob::Hit> >());
		for (auto const & h : v) fCluHits.back().push_back(h);
		fCluWeights.push_back(trackLike[i]);
	}
}
// ------------------------------------------------------

void pma::PMAlgTracker::init(const art::FindManyP< recob::Hit > & hitsFromClusters,
                             const art::FindManyP< recob::Hit > & hitsFromEmParts)
{
	mf::LogVerbatim("PMAlgTracker") << "Filter track-like clusters...";
	fCluHits.clear(); fCluHits.reserve(hitsFromClusters.size());
	fCluWeights.clear(); fCluWeights.reserve(fCluHits.size());
	size_t n = 0; // count not-empty clusters
	for (size_t i = 0; i < hitsFromClusters.size(); ++i)
	{
		auto v = hitsFromClusters.at(i);
		fCluHits.push_back(std::vector< art::Ptr<recob::Hit> >());
		for (auto const & h : v)
		{
			bool trkLike = true;
			for (size_t j = 0; j < hitsFromEmParts.size(); ++j)
			{
				auto u = hitsFromEmParts.at(j);
				for (auto const & g : u) // is hit clustered in one of em-like?
				{
					if (g.key() == h.key()) { trkLike = false; break; }
				}
				if (!trkLike) break;
			}
			if (trkLike) fCluHits.back().push_back(h);
		}
		if (!fCluHits.back().empty()) { fCluWeights.push_back(1); ++n; }
		else { fCluWeights.push_back(0); }
	}
	mf::LogVerbatim("PMAlgTracker") << "...done, " << n << " clusters for 3D tracking.";
}
// ------------------------------------------------------

double pma::PMAlgTracker::validate(pma::Track3D& trk, unsigned int testView)
{
	if ((trk.FirstElement()->GetDistToWall() < -3.0) ||
	    (trk.LastElement()->GetDistToWall() < -3.0))
	{
		mf::LogVerbatim("PMAlgTracker") << "first or last node too far out of its initial TPC";
		return 0.0;
	}

	if (testView != geo::kUnknown) { mf::LogVerbatim("PMAlgTracker") << "validation in plane: " << testView; }
	else { return 1.0; }

    double v = 0;
    switch (fValidation)
    {
        case pma::PMAlgTracker::kAdc:
            v = fProjectionMatchingAlg.validate_on_adc(trk, fAdcImages[testView], fAdcValidationThr);
            break;

        case pma::PMAlgTracker::kHits:
            v = fProjectionMatchingAlg.validate(trk, fHitMap[trk.FrontCryo()][trk.FrontTPC()][testView], testView);
            break;

        case pma::PMAlgTracker::kCalib:
            v = fProjectionMatchingAlg.validate_on_adc_test(
                trk, fAdcImages[testView], fHitMap[trk.FrontCryo()][trk.FrontTPC()][testView],
                fAdcInPassingPoints[testView], fAdcInRejectedPoints[testView]);
            break;
            
        default: throw cet::exception("pma::PMAlgTracker") << "validation mode not supported" << std::endl; break;
    }

	return v;
}
// ------------------------------------------------------

bool pma::PMAlgTracker::reassignHits_1(const std::vector< art::Ptr<recob::Hit> > & hits,
	pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2)
{
	pma::Track3D* trk1 = tracks[trk_idx].Track();

	bool result = false;
	if ((hits.size() > 1) || (dist2 > 1.0)) // min. 2 hits or single hit separated from the rest
	{
		pma::Track3D* best_trk = 0;

		size_t best_u = 0, n_max = 0;
		for (size_t u = 0; u < tracks.size(); u++)
			if (trk_idx != u)
		{
			pma::Track3D* trk2 = tracks[u].Track();
			size_t n = fProjectionMatchingAlg.testHits(*trk2, hits);
			if (n > n_max) { n_max = n; best_u = u; best_trk = trk2; }
		}

		if (best_trk && (n_max >= hits.size() / 3)) // /2
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "  Reassign(v1) " << n_max << " hits." << std::endl;

			trk1->RemoveHits(hits);
			trk1->CleanupTails();
			trk1->ShiftEndsToHits();

			pma::Track3D* ext = fProjectionMatchingAlg.extendTrack(*best_trk, hits,	false);
			ext->SortHits(); ext->ShiftEndsToHits();
			if (fProjectionMatchingAlg.isContained(*ext))
			{
				tracks[best_u].SetTrack(ext); // and this deletes best_trk stored at best_u
				result = true;
			}
			else delete ext;
		}
		else if (hits.size() >= fMinSeedSize2ndPass)
		{
			size_t minSizeCompl = hits.size() / 8;  // much smaller minimum required in complementary views
			if (minSizeCompl < 3) minSizeCompl = 3; // but at least three hits!

			geo::View_t first_view = (geo::View_t)hits.front()->WireID().Plane;
			unsigned int tpc = hits.front()->WireID().TPC;
			unsigned int cryo = hits.front()->WireID().Cryostat;

			pma::TrkCandidate candidate = matchCluster(-1, hits, minSizeCompl, tpc, cryo, first_view);

			if (candidate.IsGood())
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "  Add new track, cut hits from source track." << std::endl;
				tracks.push_back(candidate);

				trk1->RemoveHits(hits);
				trk1->CleanupTails();
				trk1->ShiftEndsToHits();
			}
		}
	}
	else if ((hits.size() == 1) || (dist2 > 2.25)) // dist > 1.5cm
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  Cut single-view isolated hit." << std::endl;
		trk1->RemoveHits(hits);
		trk1->CleanupTails();
		trk1->ShiftEndsToHits();
	}
	return result;
}

double pma::PMAlgTracker::collectSingleViewEnd(pma::Track3D & trk,
	std::vector< art::Ptr<recob::Hit> > & hits)
{
	size_t idx = 0;
	while ((idx < trk.size() - 1) && !(trk[idx]->IsEnabled()))
	{
		hits.push_back(trk[idx++]->Hit2DPtr());
	}

	double d2 = 0.0;
	if (idx > 0)
	{
		if ((idx < trk.size() - 1) &&
		    (trk[idx]->View2D() == trk[idx - 1]->View2D()))
		{
			double dprev = pma::Dist2(trk[idx]->Point3D(), trk[idx - 1]->Point3D());
			double dnext = pma::Dist2(trk[idx]->Point3D(), trk[idx + 1]->Point3D());
			if (dprev < dnext)
			{
				hits.push_back(trk[idx++]->Hit2DPtr());
			}
		}
		d2 = pma::Dist2(trk[idx]->Point3D(), trk[idx - 1]->Point3D());
	}
	return d2;
}

double pma::PMAlgTracker::collectSingleViewFront(pma::Track3D & trk,
	std::vector< art::Ptr<recob::Hit> > & hits)
{
	size_t idx = trk.size() - 1;
	while ((idx > 0) && !(trk[idx]->IsEnabled()))
	{
		hits.push_back(trk[idx--]->Hit2DPtr());
	}

	double d2 = 0.0;
	if (idx < trk.size() - 1)
	{
		if ((idx > 0) &&
		    (trk[idx]->View2D() == trk[idx + 1]->View2D()))
		{
			double dprev = pma::Dist2(trk[idx]->Point3D(), trk[idx + 1]->Point3D());
			double dnext = pma::Dist2(trk[idx]->Point3D(), trk[idx - 1]->Point3D());
			if (dprev < dnext)
			{
				hits.push_back(trk[idx--]->Hit2DPtr());
			}
		}
		d2 = pma::Dist2(trk[idx]->Point3D(), trk[idx + 1]->Point3D());
	}
	return d2;
}

bool pma::PMAlgTracker::reassignSingleViewEnds_1(pma::TrkCandidateColl & tracks)
{
	bool result = false;
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D & trk = *(tracks[t].Track());
		if (trk.size() < 6) continue;

		trk.DisableSingleViewEnds();

		std::vector< art::Ptr<recob::Hit> > hits;

		double d2 = collectSingleViewEnd(trk, hits);
		result |= reassignHits_1(hits, tracks, t, d2);

		hits.clear();

		d2 = collectSingleViewFront(trk, hits);
		result |= reassignHits_1(hits, tracks, t, d2);

		trk.SelectHits();
	}
	return result;
}
// ------------------------------------------------------

bool pma::PMAlgTracker::reassignHits_2(const std::vector< art::Ptr<recob::Hit> > & hits,
	pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2)
{
	pma::Track3D* trk1 = tracks[trk_idx].Track();

	bool result = false;
	if ((hits.size() > 1) || (dist2 > 1.0)) // min. 2 hits or single hit separated from the rest
	{
		pma::Track3D* best_trk = 0;

		size_t n_max = 0;
		for (size_t u = 0; u < tracks.size(); u++)
			if (trk_idx != u)
		{
			pma::Track3D* trk2 = tracks[u].Track();
			size_t n = fProjectionMatchingAlg.testHits(*trk2, hits, 0.5);
			if (n > n_max) { n_max = n; best_trk = trk2; }
		}

		if (best_trk && (n_max >= (size_t)(0.8 * hits.size()))) // most hits!
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "  Reassign(v2) " << n_max << " hits." << std::endl;

			trk1->RemoveHits(hits);
			trk1->CleanupTails();
			trk1->ShiftEndsToHits();

			best_trk->AddHits(hits);

			result = true;
		}
	}
	else if ((hits.size() == 1) || (dist2 > 2.25)) // dist > 1.5cm
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  Cut single-view isolated hit." << std::endl;
		trk1->RemoveHits(hits);
		trk1->CleanupTails();
		trk1->ShiftEndsToHits();

		result = true;
	}

	if (result)
	{
		// reopt trees
	}

	return result;
}

bool pma::PMAlgTracker::reassignSingleViewEnds_2(pma::TrkCandidateColl & tracks)
{
	bool result = false;
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D & trk = *(tracks[t].Track());
		if (trk.size() < 6) continue;

		trk.DisableSingleViewEnds();

		std::vector< art::Ptr<recob::Hit> > hits;

		double d2 = collectSingleViewEnd(trk, hits);
		result |= reassignHits_2(hits, tracks, t, d2);

		hits.clear();

		d2 = collectSingleViewFront(trk, hits);
		result |= reassignHits_2(hits, tracks, t, d2);

		trk.SelectHits();
	}
	return result;
}
// ------------------------------------------------------

bool pma::PMAlgTracker::areCoLinear(pma::Track3D* trk1, pma::Track3D* trk2,
	double& dist, double& cos3d, bool& reverseOrder,
	double distThr, double distThrMin,
	double distProjThr,
	double cosThr)
{
	double lmax;
	double l1 = trk1->Length();
	double l2 = trk2->Length();

	if (l1 > l2) lmax = l1;
	else lmax = l2;

	double d = lmax * distThr;
	if (d < distThrMin) d = distThrMin;

	unsigned int k = 0;
	double distFF = pma::Dist2(trk1->front()->Point3D(), trk2->front()->Point3D());
	dist = distFF;

	double distFB = pma::Dist2(trk1->front()->Point3D(), trk2->back()->Point3D());
	if (distFB < dist) { k = 1; dist = distFB; }

	double distBF = pma::Dist2(trk1->back()->Point3D(), trk2->front()->Point3D());
	if (distBF < dist) { k = 2; dist = distBF; }

	double distBB = pma::Dist2(trk1->back()->Point3D(), trk2->back()->Point3D());
	if (distBB < dist) { k = 3; dist = distBB; }

	dist = sqrt(dist);
	cos3d = 0.0;

	if (dist < d)
	{
		pma::Track3D* tmp = 0;
		switch (k) // swap or flip to get trk1 end before trk2 start
		{
			case 0:	trk1->Flip(); break;
			case 1: tmp = trk1;	trk1 = trk2; trk2 = tmp; break;
			case 2: break;
			case 3: trk2->Flip(); break;
			default: mf::LogError("PMAlgTracker") << "Should never happen.";
		}
		if (k == 1) reverseOrder = true;
		else reverseOrder = false;

		size_t nodeEndIdx = trk1->Nodes().size() - 1;

		TVector3 endpoint1 = trk1->back()->Point3D();
		TVector3 trk2front0 = trk2->Nodes()[0]->Point3D();
		TVector3 trk2front1 = trk2->Nodes()[1]->Point3D();
		TVector3 proj1 = pma::GetProjectionToSegment(endpoint1, trk2front0, trk2front1);
		double distProj1 = sqrt( pma::Dist2(endpoint1, proj1) );

		TVector3 endpoint2 = trk2->front()->Point3D();
		TVector3 trk1back0 = trk1->Nodes()[nodeEndIdx]->Point3D();
		TVector3 trk1back1 = trk1->Nodes()[nodeEndIdx - 1]->Point3D();
		TVector3 proj2 = pma::GetProjectionToSegment(endpoint2, trk1back1, trk1back0);
		double distProj2 = sqrt( pma::Dist2(endpoint2, proj2) );

		pma::Vector3D dir1 = trk1->Segments().back()->GetDirection3D();
		pma::Vector3D dir2 = trk2->Segments().front()->GetDirection3D();

		cos3d = dir1.Dot(dir2);

		if ((cos3d > cosThr) && (distProj1 < distProjThr) && (distProj2 < distProjThr))
			return true;
		else // check if parallel to wires & colinear in 2D
		{
			const double maxCosXZ = 0.996195; // 5 deg

			pma::Vector3D dir1_xz(dir1.X(), 0., dir1.Z());
			dir1_xz *= 1.0 / dir1_xz.R();

			pma::Vector3D dir2_xz(dir2.X(), 0., dir2.Z());
			dir2_xz *= 1.0 / dir2_xz.R();

			if ((fabs(dir1_xz.Z()) > maxCosXZ) && (fabs(dir2_xz.Z()) > maxCosXZ))
			{
				endpoint1.SetY(0.);
				trk2front0.SetY(0.);
				trk2front1.SetY(0.);
				proj1 = pma::GetProjectionToSegment(endpoint1, trk2front0, trk2front1);
				distProj1 = sqrt( pma::Dist2(endpoint1, proj1) );

				endpoint2.SetY(0.);
				trk1back0.SetY(0.);
				trk1back1.SetY(0.);
				proj2 = pma::GetProjectionToSegment(endpoint2, trk1back1, trk1back0);
				distProj2 = sqrt( pma::Dist2(endpoint2, proj2) );
			
				double cosThrXZ = cos(0.5 * acos(cosThr));
				double distProjThrXZ = 0.5 * distProjThr;
				double cosXZ = dir1_xz.Dot(dir2_xz);
				if ((cosXZ > cosThrXZ) && (distProj1 < distProjThrXZ) && (distProj2 < distProjThrXZ))
					return true;
			}
		}
	}
	return false;
}

bool pma::PMAlgTracker::mergeCoLinear(pma::TrkCandidateColl & tracks)
{
	double distThr = 0.25;    // max gap as a fraction of the longer track length
	double distThrMin = 0.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fMergeTransverseShift;
	double cosThr = cos(TMath::Pi() * fMergeAngle / 180.0);

	bool foundMerge = false;

	std::sort(tracks.tracks().begin(), tracks.tracks().end(), pma::bTrack3DLonger());

	bool r;
	double d, dmin, c, cmax, l, lbest;
	size_t t = 0, u = 0;
	while (t < tracks.size())
	{
		pma::Track3D* trk1 = tracks[t].Track();

		pma::Track3D* trk2 = 0;
		pma::Track3D* best_trk2 = 0;
		dmin = 1.0e12; cmax = 0; lbest = 0;
		for (u = t + 1; u < tracks.size(); u++)
		{
			trk2 = tracks[u].Track();
			if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr))
			{
				l = std::sqrt(pma::Dist2(trk2->front()->Point3D(), trk2->back()->Point3D()));
				if (((c > cmax) && (d < dmin + 0.5 * lbest)) ||
				    ((d < dmin) && (l > 1.5 * lbest)))
				{
					cmax = c; dmin = d;
					best_trk2 = trk2;
					lbest = l;
				}
			}
			trk2 = 0;
		}
		trk2 = best_trk2;

		if (trk2)
		{
			mf::LogVerbatim("PMAlgTracker") << "Merge track ("
				<< trk1->size() << ") with track (" << trk2->size() << ")";
			if (r)
			{
				fProjectionMatchingAlg.mergeTracks(*trk2, *trk1, true);
				tracks[t].SetTrack(trk2); // deletes old trk1
			}
			else
			{
				fProjectionMatchingAlg.mergeTracks(*trk1, *trk2, true);
				tracks[u].DeleteTrack();
			}
			tracks.erase_at(u);
			foundMerge = true;
		}
		else t++;
	}

	return foundMerge;
}
// ------------------------------------------------------


void pma::PMAlgTracker::freezeBranchingNodes(pma::TrkCandidateColl & tracks)
{
	for (auto const & trk : tracks.tracks())
		for (auto node : trk.Track()->Nodes())
			if (node->IsBranching()) node->SetFrozen(true);
}
void pma::PMAlgTracker::releaseAllNodes(pma::TrkCandidateColl & tracks)
{
	for (auto const & trk : tracks.tracks())
		for (auto node : trk.Track()->Nodes())
			node->SetFrozen(false);
}

void pma::PMAlgTracker::mergeCoLinear(pma::tpc_track_map& tracks)
{
	double distThr = 0.25;    // max gap as a fraction of the longer track length
	double distThrMin = 2.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fStitchTransverseShift;
	double cosThr = cos(TMath::Pi() * fStitchAngle / 180.0);

	double wallDistThr = fStitchDistToWall;
	double dfront1, dback1, dfront2, dback2;

	//for (auto & tpc_entry : tracks) freezeBranchingNodes(tpc_entry.second);

	for (auto & tpc_entry1 : tracks)
	{
		unsigned int tpc1 = tpc_entry1.first;
		pma::TrkCandidateColl & tracks1 = tpc_entry1.second;

		size_t t = 0;
		while (t < tracks1.size())
		{
			bool r, reverse = false;
			double l, lbest = 0, d, dmin = 1.0e12, c, cmax = 0.0;
			pma::Track3D* best_trk2 = 0;
			unsigned int best_tpc = 0;
			size_t best_idx = 0;

			pma::Track3D* trk1 = tracks1[t].Track();
			dfront1 = trk1->Nodes().front()->GetDistToWall();
			dback1 = trk1->Nodes().back()->GetDistToWall();
			if ((dfront1 < wallDistThr) || (dback1 < wallDistThr))
			{
				for (auto & tpc_entry2 : tracks)
				{
					unsigned int tpc2 = tpc_entry2.first;
					if (tpc2 == tpc1) continue;

					pma::TrkCandidateColl & tracks2 = tpc_entry2.second;

					for (size_t u = 0; u < tracks2.size(); u++)
					{
						pma::Track3D* trk2 = tracks2[u].Track();
						dfront2 = trk2->Nodes().front()->GetDistToWall();
						dback2 = trk2->Nodes().back()->GetDistToWall();
						if ((dfront2 < wallDistThr) || (dback2 < wallDistThr))
						{
							if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr))
							{
								l = std::sqrt(pma::Dist2(trk2->front()->Point3D(), trk2->back()->Point3D()));
								if (((c > cmax) && (d < dmin + 0.5 * lbest)) || (0.75 * l < dmin))
								{
									cmax = c; dmin = d; lbest = l;
									best_trk2 = trk2;
									best_tpc = tpc2;
									best_idx = u;
									reverse = r;
								}
							}
						}
					}
				}
			}

			if (best_trk2)
			{
				mf::LogVerbatim("PMAlgTracker") << "Merge track ("
					<< tpc1 << ":" << tracks1.size() << ":" << trk1->size() << ") with track ("
					<< best_tpc  << ":" << tracks[best_tpc].size() << ":" << best_trk2->size() << ")";
				if (reverse)
				{
					fProjectionMatchingAlg.mergeTracks(*best_trk2, *trk1, true);
					tracks1[t].SetTrack(best_trk2);
				}
				else
				{
					fProjectionMatchingAlg.mergeTracks(*trk1, *best_trk2, true);
					tracks[best_tpc][best_idx].DeleteTrack();
				}
				tracks[best_tpc].erase_at(best_idx);
			}
			else t++;
		}
	}

	//for (auto & tpc_entry : tracks) releaseAllNodes(tpc_entry.second);
}
// ------------------------------------------------------

// ------------------------------------------------------
// ------------------------------------------------------
int pma::PMAlgTracker::build(void)
{
	fInitialClusters.clear();
	fTriedClusters.clear();
	fUsedClusters.clear();

    size_t nplanes = fGeom->MaxPlanes();

	pma::tpc_track_map tracks; // track parts in tpc's

	for (auto tpc_iter = fGeom->begin_TPC_id();
	          tpc_iter != fGeom->end_TPC_id();
	          tpc_iter++)
	{
	    mf::LogVerbatim("PMAlgTracker")
	        << "Reconstruct tracks within Cryo:" << tpc_iter->Cryostat
	        << " / TPC:" << tpc_iter->TPC << ".";

	    if (fValidation != pma::PMAlgTracker::kHits) // initialize ADC images for all planes in this TPC (in "adc" and "calib")
	    {
	        mf::LogVerbatim("PMAlgTracker") << "Prepare validation ADC images...";
	        bool ok = true;
            for (size_t p = 0; p < nplanes; ++p) { ok &= fAdcImages[p].setWireDriftData(fWires, p, tpc_iter->TPC, tpc_iter->Cryostat); }
            if (ok) { mf::LogVerbatim("PMAlgTracker") << "  ...done."; }
            else { mf::LogVerbatim("PMAlgTracker") << "  ...failed."; continue; }
        }

        // find reasonably large parts
		fromMaxCluster_tpc(tracks[tpc_iter->TPC], fMinSeedSize1stPass, tpc_iter->TPC, tpc_iter->Cryostat);
		// loop again to find small things
		fromMaxCluster_tpc(tracks[tpc_iter->TPC], fMinSeedSize2ndPass, tpc_iter->TPC, tpc_iter->Cryostat);

        mf::LogVerbatim("PMAlgTracker") << "Found tracks: " << tracks[tpc_iter->TPC].size();
        if (tracks[tpc_iter->TPC].empty()) { continue; }

	    // add 3D ref.points for clean endpoints of wire-plane parallel track
		guideEndpoints(tracks[tpc_iter->TPC]);
		// try correcting single-view sections spuriously merged on 2D clusters level
		reassignSingleViewEnds_1(tracks[tpc_iter->TPC]);

	    if (fMergeWithinTPC)
	    {
			mf::LogVerbatim("PMAlgTracker") << "Merge co-linear tracks within TPC " << tpc_iter->TPC << ".";
			while (mergeCoLinear(tracks[tpc_iter->TPC]))
			{
				mf::LogVerbatim("PMAlgTracker") << "  found co-linear tracks";
			}
	    }
	}

	if (fStitchBetweenTPCs)
	{
		mf::LogVerbatim("PMAlgTracker") << "Stitch co-linear tracks between TPCs.";
		mergeCoLinear(tracks);
	}

	for (auto & tpc_entry : tracks) // put tracks in the single collection
		for (auto & trk : tpc_entry.second.tracks())
		{
			if (trk.Track()->HasTwoViews() && (trk.Track()->Nodes().size() > 1)) { fResult.push_back(trk); }
			else { trk.DeleteTrack(); }
		}

    if (fTagCosmicTracks)
    {
        mf::LogVerbatim("PMAlgTracker") << "Tag cosmic tracks activity.";
        fCosmicTagger.tag(fResult);
    }

	if (fRunVertexing)
	{
		mf::LogVerbatim("PMAlgTracker") << "Vertex finding / track-vertex reoptimization.";
		fPMAlgVertexing.run(fResult);

		//reassignSingleViewEnds(result); // final check for correct hit-track assignments
	}

	fResult.setTreeIds();

  if(fMatchT0inCPACrossing)
  {
		mf::LogVerbatim("PMAlgTracker") << "Find co-linear CPA-crossing tracks with any T0.";
    fStitcher.StitchTracksCPA(fResult);
  }

	if (fMatchT0inAPACrossing)
	{
		mf::LogVerbatim("PMAlgTracker") << "Find co-linear APA-crossing tracks with any T0.";
    fStitcher.StitchTracksAPA(fResult);
	}

  if (fTagCosmicTracks)
  {
    mf::LogVerbatim("PMAlgTracker") << "Second pass cosmic tagging for stitched tracks";
    fCosmicTagger.tag(fResult);
  }

	//double dQdxFlipThr = 0.0;
	//if (fFlipToBeam) dQdxFlipThr = 0.4;
	if (fFlipToBeam) fResult.flipTreesToCoordinate(2);        // flip the tracks / trees to the beam direction (Z)
	else if (fFlipDownward) fResult.flipTreesToCoordinate(1); // flip the tracks / trees to point downward (-Y)

	if (fAutoFlip_dQdx) fResult.flipTreesByDQdx();            // flip the tracks / trees to get best dQ/dx sequences

	fResult.setParentDaughterConnections();

	listUsedClusters();
	return fResult.size();
}
// ------------------------------------------------------
// ------------------------------------------------------

void pma::PMAlgTracker::fromMaxCluster_tpc(pma::TrkCandidateColl & result,
	size_t minBuildSize, unsigned int tpc, unsigned int cryo)
{
	fInitialClusters.clear();

	size_t minSizeCompl = minBuildSize / 8;  // smaller minimum required in complementary views
	if (minSizeCompl < 2) minSizeCompl = 2;  // but at least two hits!

	int max_first_idx = 0;
	while (max_first_idx >= 0) // loop over clusters, any view, starting from the largest
	{
		mf::LogVerbatim("PMAlgTracker") << "Find max cluster...";
		max_first_idx = maxCluster(minBuildSize, geo::kUnknown, tpc, cryo); // any view, but must be track-like
		if ((max_first_idx >= 0) && !fCluHits[max_first_idx].empty())
		{
			geo::View_t first_view = fCluHits[max_first_idx].front()->View();

			pma::TrkCandidate candidate = matchCluster(max_first_idx,
				minSizeCompl, tpc, cryo, first_view);

			if (candidate.IsGood()) result.push_back(candidate);
		}
		else mf::LogVerbatim("PMAlgTracker") << "small clusters only";
	}

	fInitialClusters.clear();
}
// ------------------------------------------------------

pma::TrkCandidate pma::PMAlgTracker::matchCluster(
	int first_clu_idx, const std::vector< art::Ptr<recob::Hit> > & first_hits,
	size_t minSizeCompl, unsigned int tpc, unsigned int cryo, geo::View_t first_view)
{
	pma::TrkCandidate result;

    for (auto av : fAvailableViews) { fTriedClusters[av].clear(); }

	if (first_clu_idx >= 0)
	{
		fTriedClusters[first_view].push_back((size_t)first_clu_idx);
		fInitialClusters.push_back((size_t)first_clu_idx);
	}

    unsigned int nFirstHits = first_hits.size(), first_plane_idx = first_hits.front()->WireID().Plane;
	mf::LogVerbatim("PMAlgTracker") << std::endl << "--- start new candidate ---";
	mf::LogVerbatim("PMAlgTracker") << "use view  *** " << first_view << " *** plane idx " << first_plane_idx << " ***  size: " << nFirstHits;

	float x, xmax = fDetProp->ConvertTicksToX(first_hits.front()->PeakTime(), first_plane_idx, tpc, cryo), xmin = xmax;
	//mf::LogVerbatim("PMAlgTracker") << "  *** x max0: " << xmax;
	for (size_t j = 1; j < first_hits.size(); ++j)
	{
		x = fDetProp->ConvertTicksToX(first_hits[j]->PeakTime(), first_plane_idx, tpc, cryo);
		if (x > xmax) { xmax = x; }
		if (x < xmin) { xmin = x; }
	}
	//mf::LogVerbatim("PMAlgTracker") << "  *** x max: " << xmax << " min:" << xmin;

	pma::TrkCandidateColl candidates; // possible solutions of the selected cluster and clusters in complementary views

	size_t imatch = 0;
	bool try_build = true;
	while (try_build) // loop over complementary views
	{
		pma::TrkCandidate candidate;
		if (first_clu_idx >= 0) candidate.Clusters().push_back((size_t)first_clu_idx);

        try_build = false;
        int idx = -1, av_idx = -1;
        unsigned int nMaxHits = 0, nHits = 0;
        unsigned int testView = geo::kUnknown, bestView = geo::kUnknown;
        for (auto av : fAvailableViews)
        {
            if (av == first_view) continue;

            av_idx = maxCluster(first_clu_idx, candidates, xmin, xmax, minSizeCompl, av, tpc, cryo);
            if (av_idx >= 0)
            {
                nHits = fCluHits[av_idx].size();
                if ((nHits > nMaxHits) && (nHits >= minSizeCompl))
                {
                    nMaxHits = nHits; idx = av_idx; bestView = av;
                    fTriedClusters[av].push_back(idx);
                    try_build = true;
                }
                else { testView = av; } // not selected, and not "first_view" -> use as a test view
            }
        }

		if (try_build)
		{
            mf::LogVerbatim("PMAlgTracker") << "--> " << imatch++ << " match with:";
		    mf::LogVerbatim("PMAlgTracker") << "    cluster in view  *** " << bestView << " ***  size: " << nMaxHits;

			if (!fGeom->TPC(tpc, cryo).HasPlane(testView)) { testView = geo::kUnknown; }
			else { mf::LogVerbatim("PMAlgTracker") << "    validation plane  *** " << testView << " ***"; }

			double m0 = 0.0, v0 = 0.0;
			double mseThr = 0.15, validThr = 0.7; // cuts for a good track candidate

			candidate.Clusters().push_back(idx);
			candidate.SetTrack(fProjectionMatchingAlg.buildTrack(first_hits, fCluHits[idx]));

			if (candidate.IsValid() && // no track if hits from 2 views do not alternate
			    fProjectionMatchingAlg.isContained(*(candidate.Track()), 2.0F)) // sticks out of TPC's?
			{
				m0 = candidate.Track()->GetMse();
				if (m0 < mseThr) // check validation only if MSE is OK - thanks for Tracy for noticing this
				{ v0 = validate(*(candidate.Track()), testView); }
			}

			if (candidate.Track() && (m0 < mseThr) && (v0 > validThr)) // good candidate, try to extend it
			{
				mf::LogVerbatim("PMAlgTracker") << "  good track candidate, MSE = " << m0 << ", v = " << v0;

				candidate.SetMse(m0);
				candidate.SetValidation(v0);
				candidate.SetGood(true);

				size_t minSize = 5;      // min size for clusters matching
				double fraction = 0.5;   // min fraction of close hits

				idx = 0;
				while (idx >= 0) // try to collect matching clusters, use **any** plane except validation
				{
					idx = matchCluster(candidate, minSize, fraction, geo::kUnknown, testView, tpc, cryo);
					if (idx >= 0)
					{
						// try building extended copy:
						//                src,        hits,      valid.plane, add nodes
						if (extendTrack(candidate, fCluHits[idx],  testView,    true))
						{
							candidate.Clusters().push_back(idx);
						}
						else idx = -1;
					}
				}

				mf::LogVerbatim("PMAlgTracker") << "merge clusters from the validation plane";
				fraction = 0.7; // only well matching the existing track

				idx = 0;
				bool extended = false;
				while ((idx >= 0) && (testView != geo::kUnknown))
				{	//                     match clusters from the plane used previously for the validation
					idx = matchCluster(candidate, minSize, fraction, testView, geo::kUnknown, tpc, cryo);
					if (idx >= 0)
					{
						// validation not checked here, no new nodes:
						if (extendTrack(candidate, fCluHits[idx], geo::kUnknown, false))
						{
							candidate.Clusters().push_back(idx);
							extended = true;
						}
						else idx = -1;
					}
				}
				// need to calculate again only if trk was extended w/o checking validation:
				if (extended) candidate.SetValidation(validate(*(candidate.Track()), testView));
			}
			else
			{
				mf::LogVerbatim("PMAlgTracker") << "track REJECTED, MSE = " << m0 << "; v = " << v0;
				candidate.SetGood(false); // save also bad matches to avoid trying again the same pair of clusters
			}
			candidates.push_back(candidate);
		}
		else
		{
			mf::LogVerbatim("PMAlgTracker") << "no matching clusters";
		}
	} // end loop over complementary views

	if (!candidates.empty()) // return best candidate, release other tracks and clusters
	{
		int best_trk = -1;
		double f, max_f = 0., min_mse = 10., max_v = 0.;
		for (size_t t = 0; t < candidates.size(); t++)
			if (candidates[t].IsGood() &&
			    (candidates[t].Track()->Nodes().size() > 1) &&
			    candidates[t].Track()->HasTwoViews())
		{
			f = fProjectionMatchingAlg.twoViewFraction(*(candidates[t].Track()));

			if ((f > max_f) || ((f == max_f) &&
				((candidates[t].Validation() > max_v) || (candidates[t].Mse() < min_mse))))
			{
				max_f = f;
				min_mse = candidates[t].Mse();
				max_v = candidates[t].Validation();
				best_trk = t;
			}
		}

		if ((best_trk > -1) && candidates[best_trk].IsGood() && (max_f > fMinTwoViewFraction))
		{
			candidates[best_trk].Track()->ShiftEndsToHits();

			for (auto c : candidates[best_trk].Clusters())
				fUsedClusters.push_back(c);

			result = candidates[best_trk];
		}

		for (size_t t = 0; t < candidates.size(); t++)
		{
			if (int(t) != best_trk) candidates[t].DeleteTrack();
		}
	}

	return result;
}
// ------------------------------------------------------

bool pma::PMAlgTracker::extendTrack(pma::TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes)
{
	double m_max = 2.0 * candidate.Mse(); // max acceptable MSE value
	if (m_max < 0.05) m_max = 0.05;     // this is still good, low MSE value

	double v_min1 = 0.98 * candidate.Validation();
	double v_min2 = 0.9 * candidate.Validation();

	pma::Track3D* copy = fProjectionMatchingAlg.extendTrack(*(candidate.Track()), hits, add_nodes);
	double m1 = copy->GetMse();
	double v1 = validate(*copy, testView);

	if (((m1 < candidate.Mse()) && (v1 >= v_min2)) ||
	    ((m1 < 0.5) && (m1 <= m_max) && (v1 >= v_min1)))
	{
		mf::LogVerbatim("PMAlgTracker")
			<< "  track EXTENDED, MSE = " << m1 << ", v = " << v1;
		candidate.SetTrack(copy);  // replace with the new track (deletes old one)
		copy->SortHits();          // sort hits in the new track

		candidate.SetMse(m1);      // save info
		candidate.SetValidation(v1);

		return true;
	}
	else
	{
		mf::LogVerbatim("PMAlgTracker")
			<< "  track NOT extended, MSE = " << m1 << ", v = " << v1;
		delete copy;
		return false;
	}
}
// ------------------------------------------------------

int pma::PMAlgTracker::matchCluster(const pma::TrkCandidate& trk,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo) const
{
	double f, fmax = 0.0;
	unsigned int n, max = 0;
	int idx = -1;
	for (size_t i = 0; i < fCluHits.size(); ++i)
	{
		if (fCluHits[i].empty()) continue;

		unsigned int view = fCluHits[i].front()->View();
		unsigned int nhits = fCluHits[i].size();

		if (has(fUsedClusters, i) ||                             // don't try already used clusters
			has(trk.Clusters(), i) ||                            // don't try clusters from this candidate
		    (view == testView) ||                                // don't use clusters from validation view
		    ((preferedView != geo::kUnknown)&&(view != preferedView)) || // only prefered view if specified
		    (nhits < minSize))                                   // skip small clusters
		    continue;

		n = fProjectionMatchingAlg.testHits(*(trk.Track()), fCluHits[i]);
		f = n / (double)nhits;
		if ((f > fraction) && (n > max))
		{
			max = n; fmax = f; idx = i;
		}
	}

	if (idx >= 0) mf::LogVerbatim("PMAlgTracker") << "max matching hits: " << max << " (" << fmax << ")";
	else mf::LogVerbatim("PMAlgTracker") << "no clusters to extend the track";

	return idx;
}
// ------------------------------------------------------

int pma::PMAlgTracker::maxCluster(int first_idx_tag,
	const pma::TrkCandidateColl & candidates,
	float xmin, float xmax, size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo) const
{
	int idx = -1;
	size_t s_max = 0, s;
	double fraction = 0.0;
	float x;

	size_t first_idx = 0;
	bool has_first = false;
	if (first_idx_tag >= 0)
	{
		first_idx = (size_t)first_idx_tag;
		has_first = true;
	}

	for (size_t i = 0; i < fCluHits.size(); ++i)
	{
		if ((fCluHits[i].size() <  min_clu_size) || (fCluHits[i].front()->View() != view) ||
		    has(fUsedClusters, i) || has(fInitialClusters, i) || has(fTriedClusters[view], i))
			continue;

		bool pair_checked = false;
		for (auto const & c : candidates.tracks())
			if (has_first && has(c.Clusters(), first_idx) && has(c.Clusters(), i))
			{
				pair_checked = true; break;
			}
		if (pair_checked) continue;
		    
		const auto & v = fCluHits[i];

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = 0;
			for (size_t j = 0; j < v.size(); ++j)
			{
				x = fDetProp->ConvertTicksToX(v[j]->PeakTime(), v[j]->WireID().Plane, tpc, cryo);
				if ((x >= xmin) && (x <= xmax)) s++;
			}

			if (s > s_max)
			{
				s_max = s; idx = i; fraction = s / (double)v.size();
			}
		}
	}
	if (fraction > 0.4) return idx;
	else return -1;
}
// ------------------------------------------------------

int pma::PMAlgTracker::maxCluster(size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo) const
{
	int idx = -1;
	size_t s_max = 0, s;

	for (size_t i = 0; i < fCluHits.size(); ++i)
	{
		const auto & v = fCluHits[i];

		if (v.empty() || (fCluWeights[i] < fTrackLikeThreshold) ||
		    has(fUsedClusters, i) || has(fInitialClusters, i) || has(fTriedClusters[view], i) ||
		   ((view != geo::kUnknown) && (v.front()->View() != view)))
		continue;

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = v.size();
			if ((s >= min_clu_size) && (s > s_max))
			{
				s_max = s; idx = i;
			}
		}
	}
	return idx;
}
// ------------------------------------------------------
// ------------------------------------------------------

void pma::PMAlgTracker::listUsedClusters(void) const
{
	mf::LogVerbatim("PMAlgTracker") << std::endl << "----------- matched clusters: -----------";
	for (size_t i = 0; i < fCluHits.size(); ++i)
		if (!fCluHits[i].empty() && has(fUsedClusters, i))
		{
			mf::LogVerbatim("PMAlgTracker")
				<< "    tpc: " << fCluHits[i].front()->WireID().TPC
				<< ";\tview: " << fCluHits[i].front()->View()
				<< ";\tsize: " << fCluHits[i].size();
		}
	mf::LogVerbatim("PMAlgTracker") << "--------- not matched clusters: ---------";
	for (size_t i = 0; i < fCluHits.size(); ++i)
		if (!fCluHits[i].empty() && !has(fUsedClusters, i))
		{
			mf::LogVerbatim("PMAlgTracker")
				<< "    tpc: " << fCluHits[i].front()->WireID().TPC
				<< ";\tview: " << fCluHits[i].front()->View()
				<< ";\tsize: " << fCluHits[i].size();
		}
	mf::LogVerbatim("PMAlgTracker") << "-----------------------------------------";
}
// ------------------------------------------------------
// ------------------------------------------------------

