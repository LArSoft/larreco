////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks and vertices using Projection Matching Algorithm,
// please see RecoAlg/ProjectionMatchingAlg.h for basics of the PMA algorithm and its settings.
//
// Progress:
//    May-June 2015:   track finding and validation, growing tracks by iterative merging of matching
//                     clusters, no attempts to build multi-track structures, however cosmic tracking
//                     works fine as they are sets of independent tracks
//    June-July 2015:  merging track parts within a single tpc and stitching tracks across tpc's
//    August 2015:     optimization of track-vertex structures (so 3D vertices are also produced)
//    November 2015:   use track-shower splitting at 2D level, then tag low-E EM cascades in 3D
//                     note: the splitter is not finished and not as good as we want it
//    January 2016:    output of track-vertex finding as a tree of PFParticles, refined vertexing
//                     code, put vertex at front of each track, flip whole track structures to
//                     selected, direction (beam/down/dQdx), use any pattern reco stored in
//                     PFParticles as a source of associated clusters
//    Mar-Apr 2016:    kinks finding, EM shower direction reconstructed for PFPaarticles tagged as
//                     electrons
//    July 2016:       redesign module: extract trajectory fitting-only to separate module, move
//                     tracking functionality to algorithm classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/TrackHitMeta.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/AnalysisBase/T0.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
//#include "lardata/Utilities/PtrMaker.h"

#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlgTracking.h"
#include "larreco/RecoAlg/PMAlgVertexing.h"

#include <memory>

namespace trkf {

class PMAlgTrackMaker : public art::EDProducer {
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<pma::ProjectionMatchingAlg::Config> ProjectionMatchingAlg {
			Name("ProjectionMatchingAlg")
		};

		fhicl::Table<pma::PMAlgTracker::Config> PMAlgTracking {
			Name("PMAlgTracking")
		};

		fhicl::Table<pma::PMAlgVertexing::Config> PMAlgVertexing {
			Name("PMAlgVertexing")
		};

		fhicl::Atom<bool> SaveOnlyBranchingVtx {
			Name("SaveOnlyBranchingVtx"),
			Comment("use true to save only vertices interconnecting many tracks, otherwise vertex is added to the front of each track")
		};

		fhicl::Atom<bool> SavePmaNodes {
			Name("SavePmaNodes"),
			Comment("save track nodes (only for algorithm development purposes)")
		};

		fhicl::Atom<art::InputTag> HitModuleLabel {
			Name("HitModuleLabel"),
			Comment("tag of unclustered hits, which were used to validate tracks")
		};

		fhicl::Atom<art::InputTag> ClusterModuleLabel {
			Name("ClusterModuleLabel"),
			Comment("tag of cluster collection, these clusters are used for track building")
		};

		fhicl::Atom<art::InputTag> EmClusterModuleLabel {
			Name("EmClusterModuleLabel"),
			Comment("EM-like clusters, will be excluded from tracking if provided")
		};
    };
    using Parameters = art::EDProducer::Table<Config>;

	explicit PMAlgTrackMaker(Parameters const& config);

	PMAlgTrackMaker(PMAlgTrackMaker const &) = delete;
	PMAlgTrackMaker(PMAlgTrackMaker &&) = delete;
	PMAlgTrackMaker & operator = (PMAlgTrackMaker const &) = delete;
	PMAlgTrackMaker & operator = (PMAlgTrackMaker &&) = delete;

	void produce(art::Event & e) override;

private:
	// ******************** fcl parameters **********************
	art::InputTag fHitModuleLabel; // tag for hits collection (used for trk validation)
	art::InputTag fCluModuleLabel; // tag for input cluster collection
	art::InputTag fEmModuleLabel;  // tag for em-like cluster collection

	pma::ProjectionMatchingAlg::Config fPmaConfig;
	pma::PMAlgTracker::Config fPmaTrackerConfig;
	pma::PMAlgVertexing::Config fPmaVtxConfig;

	bool fSaveOnlyBranchingVtx;  // for debugging, save only vertices which connect many tracks
	bool fSavePmaNodes;          // for debugging, save only track nodes

	// ********** instance names (collections, assns) ************
	static const std::string kKinksName;        // kinks on tracks
	static const std::string kNodesName;        // pma nodes

	// *********************** geometry **************************
	art::ServiceHandle< geo::Geometry > fGeom;
};
// -------------------------------------------------------------
const std::string PMAlgTrackMaker::kKinksName = "kink";
const std::string PMAlgTrackMaker::kNodesName = "node";
// -------------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(PMAlgTrackMaker::Parameters const& config) :
	fHitModuleLabel(config().HitModuleLabel()),
	fCluModuleLabel(config().ClusterModuleLabel()),
	fEmModuleLabel(config().EmClusterModuleLabel()),

	fPmaConfig(config().ProjectionMatchingAlg()),
	fPmaTrackerConfig(config().PMAlgTracking()),
	fPmaVtxConfig(config().PMAlgVertexing()),

	fSaveOnlyBranchingVtx(config().SaveOnlyBranchingVtx()),
	fSavePmaNodes(config().SavePmaNodes())
{
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< std::vector<recob::Vertex> >(); // no instance name for interaction vertices
	produces< std::vector<recob::Vertex> >(kKinksName); // collection of kinks on tracks
	produces< std::vector<recob::Vertex> >(kNodesName); // collection of pma nodes
	produces< std::vector<anab::T0> >();

	produces< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
	produces< art::Assns<recob::Vertex, recob::Track> >(); // no instance name for assns of tracks to interaction vertices
	produces< art::Assns<recob::Track, recob::Vertex> >(kKinksName);  // assns of kinks to tracks
	produces< art::Assns<recob::Track, anab::T0> >();

	produces< std::vector<recob::PFParticle> >();
	produces< art::Assns<recob::PFParticle, recob::Cluster> >();
	produces< art::Assns<recob::PFParticle, recob::Vertex> >();
	produces< art::Assns<recob::PFParticle, recob::Track> >();
}
// ------------------------------------------------------

void PMAlgTrackMaker::produce(art::Event& evt)
{
	// ---------------- Create data products ------------------
	auto tracks = std::make_unique< std::vector< recob::Track > >();
	auto allsp = std::make_unique< std::vector< recob::SpacePoint > >();
	auto vtxs = std::make_unique< std::vector< recob::Vertex > >();  // interaction vertices
	auto kinks = std::make_unique< std::vector< recob::Vertex > >(); // kinks on tracks (no new particles start in kinks)
	auto nodes = std::make_unique< std::vector< recob::Vertex > >(); // pma nodes
	auto t0s = std::make_unique< std::vector< anab::T0 > >();

	auto trk2hit_oldway = std::make_unique< art::Assns< recob::Track, recob::Hit > >(); // ****** REMEMBER to remove when FindMany improved ******
	auto trk2hit = std::make_unique< art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta > >();

	auto trk2sp = std::make_unique< art::Assns< recob::Track, recob::SpacePoint > >();
	auto trk2t0 = std::make_unique< art::Assns< recob::Track, anab::T0 > >();

	auto sp2hit = std::make_unique< art::Assns< recob::SpacePoint, recob::Hit > >();
	auto vtx2trk = std::make_unique< art::Assns< recob::Vertex, recob::Track > >();  // one or more tracks (particles) start in the vertex
	auto trk2kink = std::make_unique< art::Assns< recob::Track, recob::Vertex > >(); // one or more kinks on the track

	auto pfps = std::make_unique< std::vector< recob::PFParticle > >();

    auto pfp2clu = std::make_unique< art::Assns<recob::PFParticle, recob::Cluster> >();
    auto pfp2vtx = std::make_unique< art::Assns<recob::PFParticle, recob::Vertex> >();
	auto pfp2trk = std::make_unique< art::Assns< recob::PFParticle, recob::Track > >();


	// -------------- Collect hits and clusters ---------------
	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle, splitCluHandle;
	std::vector< art::Ptr<recob::Hit> > allhitlist;

	if (!(evt.getByLabel(fHitModuleLabel, allHitListHandle) && // all hits associated used to create clusters
    	  evt.getByLabel(fCluModuleLabel, cluListHandle)))     // clusters used to build 3D tracks
	{
		mf::LogError("PMAlgTrackMaker") << "Not all required data products found in the event.";
		return;
	}

	art::fill_ptr_vector(allhitlist, allHitListHandle);
	art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fCluModuleLabel);

	// -------------- PMA Tracker for this event --------------
	auto pmalgTracker = pma::PMAlgTracker(allhitlist,
		fPmaConfig, fPmaTrackerConfig, fPmaVtxConfig);


	if (fEmModuleLabel != "") // --- Exclude EM parts ---------
	{
		if (evt.getByLabel(fEmModuleLabel, splitCluHandle))
		{
			art::FindManyP< recob::Hit > hitsFromEmParts(splitCluHandle, evt, fEmModuleLabel);
			pmalgTracker.init(hitsFromClusters, hitsFromEmParts);
		}
		else
		{
			mf::LogError("PMAlgTrackMaker") << "EM-tagged clusters not found in the event.";
			return;
		}
	}
	else // ------------------------ Use ALL clusters ---------
	{
		pmalgTracker.init(hitsFromClusters);
	}

	// ------------------ Do the job here: --------------------
	int retCode = pmalgTracker.build();
	// --------------------------------------------------------
	switch (retCode)
	{
		case -2: mf::LogError("Summary") << "problem"; break;
		case -1: mf::LogWarning("Summary") << "no input"; break;
		case  0: mf::LogVerbatim("Summary") << "no tracks done"; break;
		default:
			if (retCode < 0) mf::LogVerbatim("Summary") << "unknown result";
			else if (retCode == 1) mf::LogVerbatim("Summary") << retCode << " track ready";
			else mf::LogVerbatim("Summary") << retCode << " tracks ready";
			break;
	}

	// ---------- Translate output to data products: ----------
	auto const & result = pmalgTracker.Result();
	if (!result.empty()) // ok, there is something to save
	{
		const detinfo::DetectorProperties* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6];
		for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

		//auto const make_trkptr = lar::PtrMaker<recob::Track>(evt, *this); // PtrMaker Step #1
		//auto const make_t0ptr = lar::PtrMaker<anab::T0>(evt, *this);

		tracks->reserve(result.size());
		for (size_t trkIndex = 0; trkIndex < result.size(); ++trkIndex)
		{
			pma::Track3D* trk = result[trkIndex].Track();

			trk->SelectHits();  // just in case, set all to enabled
			unsigned int itpc = trk->FrontTPC(), icryo = trk->FrontCryo();
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kU)) trk->CompleteMissingWires(geo::kU);
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kV)) trk->CompleteMissingWires(geo::kV);
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kZ)) trk->CompleteMissingWires(geo::kZ);

			tracks->push_back(pma::convertFrom(*trk, trkIndex));

			//auto const trkPtr = make_trkptr(tracks->size() - 1); // PtrMaker Step #2

			double xShift = trk->GetXShift();
			if (xShift > 0.0)
			{
				double tisk2time = 1.0; // what is the coefficient, offset?
				double t0time = tisk2time * xShift / detProp->GetXTicksCoefficient(trk->FrontTPC(), trk->FrontCryo());

				// TriggBits=3 means from 3d reco (0,1,2 mean something else)
				t0s->push_back(anab::T0(t0time, 0, 3, tracks->back().ID()));

				util::CreateAssn(*this, evt, *tracks, *t0s, *trk2t0, t0s->size() - 1, t0s->size());
				//auto const t0Ptr = make_t0ptr(t0s->size() - 1);  // PtrMaker Step #3
				//trk2t0->addSingle(trkPtr, t0Ptr);
			}

			size_t trkIdx = tracks->size() - 1; // stuff for assns:
			art::ProductID trkId = getProductID< std::vector<recob::Track> >(evt);
			art::Ptr<recob::Track> trkPtr(trkId, trkIdx, evt.productGetter(trkId));

			// which idx from start, except disabled, really....
			unsigned int hIdxs[trk->size()];
			for (size_t h = 0, cnt = 0; h < trk->size(); h++)
			{
				if ((*trk)[h]->IsEnabled()) hIdxs[h] = cnt++;
				else hIdxs[h] = 0;
			}

			art::PtrVector< recob::Hit > sp_hits;
			spStart = allsp->size();
			for (int h = trk->size() - 1; h >= 0; h--)
			{
				pma::Hit3D* h3d = (*trk)[h];
				if (!h3d->IsEnabled()) continue;

				recob::TrackHitMeta metadata(hIdxs[h], h3d->Dx());
				trk2hit->addSingle(trkPtr, h3d->Hit2DPtr(), metadata);
				trk2hit_oldway->addSingle(trkPtr, h3d->Hit2DPtr()); // ****** REMEMBER to remove when FindMany improved ******

				double hx = h3d->Point3D().X() + xShift;
				double hy = h3d->Point3D().Y();
				double hz = h3d->Point3D().Z();

				if ((h == 0) || (sp_pos[0] != hx) || (sp_pos[1] != hy) || (sp_pos[2] != hz))
				{
					if (sp_hits.size()) // hits assigned to the previous sp
					{
						util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
						sp_hits.clear();
					}
					sp_pos[0] = hx; sp_pos[1] = hy; sp_pos[2] = hz;
					allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
				}
				sp_hits.push_back(h3d->Hit2DPtr());
			}

			if (sp_hits.size()) // hits assigned to the last sp
			{
				util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
			}
			spEnd = allsp->size();

			if (spEnd > spStart) util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);
		}

		auto pfpid = getProductID< std::vector<recob::PFParticle> >(evt);
		auto vid = getProductID< std::vector<recob::Vertex> >(evt);
		auto kid = getProductID< std::vector<recob::Vertex> >(evt, kKinksName);
		auto const* kinkGetter = evt.productGetter(kid);

		auto tid = getProductID< std::vector<recob::Track> >(evt);
		auto const* trkGetter = evt.productGetter(tid);

		auto vsel = fPMAlgVertexing.getVertices(result, fSaveOnlyBranchingVtx); // vtx pos's with vector of connected track idxs
		auto ksel = fPMAlgVertexing.getKinks(result); // pairs of kink position - associated track idx 
		std::map< size_t, art::Ptr<recob::Vertex> > frontVtxs; // front vertex ptr for each track index

		if (fPmaTrackerConfig.RunVertexing()) // save vertices and vtx-trk assns
		{
			double xyz[3];
			for (auto const & v : vsel)
			{
				xyz[0] = v.first.X(); xyz[1] = v.first.Y(); xyz[2] = v.first.Z();
				mf::LogVerbatim("Summary")
					<< "  vtx:" << xyz[0] << ":" << xyz[1] << ":" << xyz[2]
					<< "  (" << v.second.size() << " tracks)";

				size_t vidx = vtxs->size();
				vtxs->push_back(recob::Vertex(xyz, vidx));

				art::Ptr<recob::Vertex> vptr(vid, vidx, evt.productGetter(vid));
				if (vptr.isNull()) mf::LogWarning("PMAlgTrackMaker") << "Vertex ptr is null.";
				if (!v.second.empty())
				{
					for (const auto & vEntry : v.second)
					{
						size_t tidx = vEntry.first;
						bool isFront = vEntry.second;

						if (isFront) frontVtxs[tidx] = vptr; // keep ptr of the front vtx

						art::Ptr<recob::Track> tptr(tid, tidx, trkGetter);
						vtx2trk->addSingle(vptr, tptr);
					}
				}
				else mf::LogWarning("PMAlgTrackMaker") << "No tracks found at this vertex.";
			}
			mf::LogVerbatim("Summary") << vtxs->size() << " vertices ready";

			for (auto const & k : ksel)
			{
				xyz[0] = k.first.X(); xyz[1] = k.first.Y(); xyz[2] = k.first.Z();
				mf::LogVerbatim("Summary") << "  kink:" << xyz[0] << ":" << xyz[1] << ":" << xyz[2];

				size_t kidx = kinks->size();
				size_t tidx = k.second; // track idx on which this kink was found

				kinks->push_back(recob::Vertex(xyz, tidx)); // save index of track (will have color of trk in evd)

				art::Ptr<recob::Track> tptr(tid, tidx, trkGetter);
				art::Ptr<recob::Vertex> kptr(kid, kidx, kinkGetter);
				trk2kink->addSingle(tptr, kptr);
			}
			mf::LogVerbatim("Summary") << ksel.size() << " kinks ready";
		}

		if (fSavePmaNodes)
		{
			double xyz[3];
			for (size_t t = 0; t < result.size(); ++t)
			{
				auto const & trk = *(result[t].Track());
				for (auto const * node : trk.Nodes())
				{
					xyz[0] = node->Point3D().X(); xyz[1] = node->Point3D().Y(); xyz[2] = node->Point3D().Z();
					nodes->push_back(recob::Vertex(xyz, t));
				}
			}
		}

		// first particle, to be replaced with nu reco when possible
		pfps->emplace_back(0, 0, 0, std::vector< size_t >());
		for (size_t t = 0; t < result.size(); ++t)
		{
			size_t parentIdx = 0;
			if (result[t].Parent() >= 0) parentIdx = (size_t)result[t].Parent() + 1;

			std::vector< size_t > daughterIdxs;
			for (size_t idx : result[t].Daughters()) daughterIdxs.push_back(idx + 1);

			size_t pfpidx = pfps->size();
			pfps->emplace_back(0, pfpidx, parentIdx, daughterIdxs);

			art::Ptr<recob::PFParticle> pfpptr(pfpid, pfpidx, evt.productGetter(pfpid));
			art::Ptr<recob::Track> tptr(tid, t, trkGetter);
			pfp2trk->addSingle(pfpptr, tptr);

			if (fRunVertexing) // vertexing was used, so add assns to front vertex of each particle
			{
				art::Ptr<recob::Vertex> vptr = frontVtxs[t];
				if (!vptr.isNull()) pfp2vtx->addSingle(pfpptr, vptr);
				else mf::LogWarning("PMAlgTrackMaker") << "Front vertex for PFParticle is missing.";
			}
		}
		mf::LogVerbatim("Summary") << pfps->size() << " PFParticles created";
	}


	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(vtxs));
	evt.put(std::move(kinks), kKinksName);
	evt.put(std::move(nodes), kNodesName);
	evt.put(std::move(t0s));

	evt.put(std::move(trk2hit_oldway)); // ****** REMEMBER to remove when FindMany improved ******
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(trk2t0));

	evt.put(std::move(sp2hit));
	evt.put(std::move(vtx2trk));
	evt.put(std::move(trk2kink), kKinksName);

	evt.put(std::move(pfps));
	evt.put(std::move(pfp2clu));
	evt.put(std::move(pfp2vtx));
	evt.put(std::move(pfp2trk));
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

DEFINE_ART_MODULE(PMAlgTrackMaker)

} // namespace trkf

