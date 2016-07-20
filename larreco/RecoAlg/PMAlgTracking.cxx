////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTracking
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), June 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgTracking.h"

#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

recob::Track pma::convertFrom(const pma::Track3D& src, unsigned int tidx)
{
	std::vector< TVector3 > xyz, dircos;
	xyz.reserve(src.size()); dircos.reserve(src.size());

	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = src.FrontCryo();
	unsigned int tpc = src.FrontTPC();

	art::ServiceHandle<geo::Geometry> geom;
	std::map< unsigned int, pma::dedx_map > src_dQdx;
	if (geom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		src_dQdx[geo::kU] = pma::dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kU], geo::kU);
	}
	if (geom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		src_dQdx[geo::kV] = pma::dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kV], geo::kV);
	}
	if (geom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		src_dQdx[geo::kZ] = pma::dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kZ], geo::kZ);
	}

	TVector3 p3d;
	double xshift = src.GetXShift();
	bool has_shift = (xshift != 0.0);
	for (size_t i = 0; i < src.size(); i++)
		if (src[i]->IsEnabled())
	{
		p3d = src[i]->Point3D();
		if (has_shift) p3d.SetX(p3d.X() + xshift);
		xyz.push_back(p3d);

		if (i < src.size() - 1)
		{
			size_t j = i + 1;
			double mag = 0.0;
			TVector3 dc(0., 0., 0.);

			while ((mag == 0.0) && (j < src.size()))
			{
				dc = src[j]->Point3D();
				dc -= src[i]->Point3D();
				mag = dc.Mag();
				j++;
			}

			if (mag > 0.0) dc *= 1.0 / mag;
			else if (!dircos.empty()) dc = dircos.back();

			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());

		double dQ = 0., dx = 0.;
		dst_dQdx[geo::kU].push_back(0.);
		dst_dQdx[geo::kV].push_back(0.);
		dst_dQdx[geo::kZ].push_back(0.);

		double dQdx;
		for (auto const& m : src_dQdx)
		{
			auto it = m.second.find(i);
			if (it != m.second.end())
			{
				dQ = it->second[5];
				dx = it->second[6];
				if (dx > 0.) dQdx = dQ/dx;
				else dQdx = 0.;

				size_t backIdx = dst_dQdx[m.first].size() - 1;
				dst_dQdx[m.first][backIdx] = dQdx;

				break;
			}
		}
	}

	// 0 is track-like (long and/or very straight, well matching 2D hits);
	// 0x10000 is EM shower-like trajectory *** to be replaced with a better tag ***
	unsigned int pidTag = 0;
	if (src.GetTag() == pma::Track3D::kEmLike) pidTag = 0x10000;

	if (xyz.size() != dircos.size())
	{
		mf::LogError("PMAlgTrackMaker") << "pma::Track3D to recob::Track conversion problem.";
	}

	return recob::Track(xyz, dircos, dst_dQdx, std::vector< double >(2, util::kBogusD), tidx + pidTag);
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

void pma::PMAlgTrackingBase::guideEndpoints(void)
{
	for (auto const & t : fResult.tracks())
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
                              const pma::PMAlgVertexing::Config& pmvtxConfig) :
	PMAlgTrackingBase(allhitlist, pmalgConfig, pmvtxConfig)
{
	mf::LogVerbatim("PMAlgTrajFitter") << "Found " << allhitlist.size() << "hits in the event.";
	mf::LogVerbatim("PMAlgTrajFitter") << "Sort hits by clusters assigned to PFParticles...";

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

		if (vtxFromPfps.at(i).size())
		{
			double xyz[3];
			vtxFromPfps.at(i).front()->XYZ(xyz);
			fPfpVtx[i] = pma::Vector3D(xyz[0], xyz[1], xyz[2]);
		}
	}

	mf::LogVerbatim("PMAlgTrajFitter") << "...done, "
		<< fCluHits.size() << " clusters from "
		<< fPfpClusters.size() << " pfparticles for 3D tracking.";
}
// ------------------------------------------------------

int pma::PMAlgFitter::build(bool runVertexing,
	const std::vector<int> & trackingOnlyPdg,
	const std::vector<int> & trackingSkipPdg)
{
    if (!fPfpClusters.empty() && !fCluHits.empty())
    {
			// build pm tracks
			buildTracks(trackingOnlyPdg, trackingSkipPdg);

			// add 3D ref.points for clean endpoints of wire-plae parallel tracks
			guideEndpoints();

			if (runVertexing) fPMAlgVertexing.run(fResult);

			// build segment of shower
			buildShowers(trackingOnlyPdg, trackingSkipPdg);
    }
    else
    {
        mf::LogWarning("PMAlgTrajFitter") << "no clusters, no pfparticles";
        return -1;
    }
    
    return fResult.size();
}
// ------------------------------------------------------

void pma::PMAlgFitter::buildTracks(
	const std::vector<int> & trackingOnlyPdg,
	const std::vector<int> & trackingSkipPdg)
{
		bool skipPdg = true;
		if (!trackingSkipPdg.empty() && (trackingSkipPdg.front() == 0)) skipPdg = false;

		bool selectPdg = true;
		if (!trackingOnlyPdg.empty() && (trackingOnlyPdg.front() == 0)) selectPdg = false;

		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (pdg == 11) continue;
			if (skipPdg && has(trackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(trackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgTrajFitter") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

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
	   			fResult.push_back(candidate);
			}
			else
			{
				candidate.DeleteTrack();
			}
		}
}
// ------------------------------------------------------

void pma::PMAlgFitter::buildShowers(
	const std::vector<int> & trackingOnlyPdg,
	const std::vector<int> & trackingSkipPdg)
{
		bool skipPdg = true;
		if (!trackingSkipPdg.empty() && (trackingSkipPdg.front() == 0))
			skipPdg = false;

		bool selectPdg = true;
		if (!trackingOnlyPdg.empty() && (trackingOnlyPdg.front() == 0))
			selectPdg = false;

		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (pdg != 11) continue;
			if (skipPdg && has(trackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(trackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgTrajFitter") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

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

			mf::LogVerbatim("PMAlgTrajFitter") << "building..." << ", pdg:" << pdg;

			auto search = fPfpVtx.find(pfPartIdx);
			if (search != fPfpVtx.end())
			{
				candidate.SetTrack(fProjectionMatchingAlg.buildShowerSeg(allHits, search->second));
				if (candidate.IsValid()
						&& candidate.Track()->HasTwoViews() 
						&& (candidate.Track()->Nodes().size() > 1)) 
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
