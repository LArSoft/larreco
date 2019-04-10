////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrajFitter
// Module Type: producer
// File:        PMAlgTrajFitter_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks using Projection Matching Algorithm, please see
// RecoAlg/ProjectionMatchingAlg.h for basics of the PMA algorithm and its settings.
//
// Progress:
//    July 2016:   fit-only module separated from PMAlgTrackMaker_module, common tracking
//                 functionality moved to algorithm class;
//                 note, that some part of vertexing can be still reasonable in this module:
//                  - use hierarchy of input PFP's to create vertices, join tracks, reoptimize
//                  - detect vertices/kinks on fitted trajectories
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "canvas/Utilities/InputTag.h"
//#include "art/Persistency/Common/PtrMaker.h"

#include "larreco/RecoAlg/PMAlgTracking.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include <memory>

namespace trkf {

class PMAlgTrajFitter : public art::EDProducer {
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<pma::ProjectionMatchingAlg::Config> ProjectionMatchingAlg {
			Name("ProjectionMatchingAlg")
		};

		fhicl::Table<pma::PMAlgFitter::Config> PMAlgFitting {
			Name("PMAlgFitting")
		};

		fhicl::Table<pma::PMAlgVertexing::Config> PMAlgVertexing {
			Name("PMAlgVertexing")
		};

		fhicl::Atom<art::InputTag> HitModuleLabel {
			Name("HitModuleLabel"),
			Comment("tag of unclustered hits, which were used to produce PFPs and clusters")
		};

		fhicl::Atom<art::InputTag> PfpModuleLabel {
			Name("PfpModuleLabel"),
			Comment("tag of the input PFParticles and associated clusters")
		};

		fhicl::Atom<bool> SaveOnlyBranchingVtx {
			Name("SaveOnlyBranchingVtx"),
			Comment("use true to save only vertices interconnecting many tracks, otherwise vertex is added to the front of each track")
		};

		fhicl::Atom<bool> SavePmaNodes {
			Name("SavePmaNodes"),
			Comment("save track nodes (only for algorithm development purposes)")
		};
    };
    using Parameters = art::EDProducer::Table<Config>;

	explicit PMAlgTrajFitter(Parameters const& config);

	PMAlgTrajFitter(PMAlgTrajFitter const &) = delete;
	PMAlgTrajFitter(PMAlgTrajFitter &&) = delete;
	PMAlgTrajFitter & operator = (PMAlgTrajFitter const &) = delete;
	PMAlgTrajFitter & operator = (PMAlgTrajFitter &&) = delete;

	void produce(art::Event & e) override;

private:
  // ******************** fcl parameters ***********************
  art::InputTag fHitModuleLabel; // tag for unclustered hit collection
  art::InputTag fPfpModuleLabel; // tag for PFParticle and cluster collections

  pma::ProjectionMatchingAlg::Config fPmaConfig;
  pma::PMAlgFitter::Config fPmaFitterConfig;
  pma::PMAlgVertexing::Config fPmaVtxConfig;

  bool fSaveOnlyBranchingVtx;  // for debugging, save only vertices which connect many tracks
  bool fSavePmaNodes;          // for debugging, save only track nodes

  // ********** instance names (collections, assns) ************
  static const std::string kKinksName;        // kinks on tracks
  static const std::string kNodesName;        // pma nodes

  // *********************** geometry **************************
  art::ServiceHandle<geo::Geometry const> fGeom;
};
// -------------------------------------------------------------
const std::string PMAlgTrajFitter::kKinksName = "kink";
const std::string PMAlgTrajFitter::kNodesName = "node";
// -------------------------------------------------------------

PMAlgTrajFitter::PMAlgTrajFitter(PMAlgTrajFitter::Parameters const& config) :
    EDProducer{config},
	fHitModuleLabel(config().HitModuleLabel()),
	fPfpModuleLabel(config().PfpModuleLabel()),

	fPmaConfig(config().ProjectionMatchingAlg()),
	fPmaFitterConfig(config().PMAlgFitting()),
	fPmaVtxConfig(config().PMAlgVertexing()),

	fSaveOnlyBranchingVtx(config().SaveOnlyBranchingVtx()),
	fSavePmaNodes(config().SavePmaNodes())
{
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< std::vector<recob::Vertex> >(); // no instance name for interaction vertices
	produces< std::vector<recob::Vertex> >(kKinksName); // collection of kinks on tracks
	produces< std::vector<recob::Vertex> >(kNodesName); // collection of pma nodes

	produces< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
	produces< art::Assns<recob::Vertex, recob::Track> >(); // no instance name for assns of tracks to interaction vertices
	produces< art::Assns<recob::Track, recob::Vertex> >(kKinksName);  // assns of kinks to tracks

	produces< art::Assns<recob::PFParticle, recob::Track> >();
}
// ------------------------------------------------------

void PMAlgTrajFitter::produce(art::Event& evt)
{
	// ---------------- Create data products ------------------
	auto tracks = std::make_unique< std::vector<recob::Track> >();
	auto allsp = std::make_unique< std::vector<recob::SpacePoint> >();
	auto vtxs = std::make_unique< std::vector<recob::Vertex> >();  // interaction vertices
	auto kinks = std::make_unique< std::vector<recob::Vertex> >(); // kinks on tracks (no new particles start in kinks)
	auto nodes = std::make_unique< std::vector<recob::Vertex> >(); // pma nodes

	auto trk2hit_oldway = std::make_unique< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	auto trk2hit = std::make_unique< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	auto trk2sp = std::make_unique< art::Assns<recob::Track, recob::SpacePoint> >();

	auto sp2hit = std::make_unique< art::Assns<recob::SpacePoint, recob::Hit> >();
	auto vtx2trk = std::make_unique< art::Assns<recob::Vertex, recob::Track> >(); // one or more tracks (particles) start in the vertex
	auto trk2kink = std::make_unique< art::Assns<recob::Track, recob::Vertex> >(); // one or more kinks on the track

	auto pfp2trk = std::make_unique< art::Assns< recob::PFParticle, recob::Track> >();

	// ------------------- Collect inputs ---------------------
	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
	std::vector< art::Ptr<recob::Hit> > allhitlist;
	if (!(evt.getByLabel(fHitModuleLabel, allHitListHandle) &&  // all hits used to make clusters and PFParticles
	      evt.getByLabel(fPfpModuleLabel, cluListHandle) &&     // clusters associated to PFParticles
	      evt.getByLabel(fPfpModuleLabel, pfparticleHandle)))   // and finally PFParticles
	{
		throw cet::exception("PMAlgTrajFitter") << "Not all required data products found in the event." << std::endl;
	}

	art::fill_ptr_vector(allhitlist, allHitListHandle);

	art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fPfpModuleLabel);
	art::FindManyP< recob::Cluster > clustersFromPfps(pfparticleHandle, evt, fPfpModuleLabel);
	art::FindManyP< recob::Vertex > vtxFromPfps(pfparticleHandle, evt, fPfpModuleLabel);

	// -------------- PMA Fitter for this event ---------------
	auto pmalgFitter = pma::PMAlgFitter(allhitlist,
		*cluListHandle, *pfparticleHandle,
		hitsFromClusters, clustersFromPfps, vtxFromPfps,
		fPmaConfig, fPmaFitterConfig, fPmaVtxConfig);

	// ------------------ Do the job here: --------------------
	int retCode = pmalgFitter.build();
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
	auto const & result = pmalgFitter.result();
	if (!result.empty()) // ok, there is something to save
	{
		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6];
		for (size_t i = 0; i < 3; i++) sp_pos[i] = 1.0;
		for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

		// use the following to create PFParticle <--> Track associations;
		std::map< size_t, std::vector< art::Ptr<recob::Track> > > pfPartToTrackVecMap;

		//auto const make_trkptr = art::PtrMaker<recob::Track>(evt, *this); // PtrMaker Step #1

		tracks->reserve(result.size());
		for (size_t trkIndex = 0; trkIndex < result.size(); ++trkIndex)
		{
			pma::Track3D* trk = result[trkIndex].Track();

			trk->SelectHits();  // just in case, set all to enabled
			unsigned int itpc = trk->FrontTPC(), icryo = trk->FrontCryo();
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kU)) trk->CompleteMissingWires(geo::kU);
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kV)) trk->CompleteMissingWires(geo::kV);
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kZ)) trk->CompleteMissingWires(geo::kZ);

			//gc: make sure no tracks are created with less than 2 points
			if (trk->size()<2) continue;

			tracks->push_back(pma::convertFrom(*trk, trkIndex));

			//auto const trkPtr = make_trkptr(tracks->size() - 1); // PtrMaker Step #2

			size_t trkIdx = tracks->size() - 1; // stuff for assns:
                        art::ProductID trkId = evt.getProductID< std::vector<recob::Track> >();
			art::Ptr<recob::Track> trkPtr(trkId, trkIdx, evt.productGetter(trkId));

			//gc: save associated hits in the same order as trajectory points
                        for (size_t h = 0, cnt = 0; h < trk->size(); h++)
			{
				pma::Hit3D* h3d = (*trk)[h];
				if (!h3d->IsEnabled()) continue;

				recob::TrackHitMeta metadata(cnt++, h3d->Dx());
				trk2hit->addSingle(trkPtr, h3d->Hit2DPtr(), metadata);
				trk2hit_oldway->addSingle(trkPtr, h3d->Hit2DPtr()); // ****** REMEMBER to remove when FindMany improved ******
			}

			art::PtrVector< recob::Hit > sp_hits;
			spStart = allsp->size();
			for (size_t h = 0; h < trk->size(); ++h)
			{
				pma::Hit3D* h3d = (*trk)[h];
				if (!h3d->IsEnabled()) continue;

				double hx = h3d->Point3D().X();
				double hy = h3d->Point3D().Y();
				double hz = h3d->Point3D().Z();

                                if ((h == 0) ||
                                    (std::fabs(sp_pos[0] - hx) > 1.0e-5) ||
                                    (std::fabs(sp_pos[1] - hy) > 1.0e-5) ||
                                    (std::fabs(sp_pos[2] - hz) > 1.0e-5))
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

			// if there is a PFParticle collection then recover PFParticle and add info to map
			if (result[trkIndex].Key() > -1)
			{
				size_t trackIdx = tracks->size() - 1;
                                art::ProductID trackId = evt.getProductID< std::vector<recob::Track> >();
				art::Ptr<recob::Track> trackPtr(trackId, trackIdx, evt.productGetter(trackId));
				pfPartToTrackVecMap[result[trkIndex].Key()].push_back(trackPtr);
			}
		}

                auto vid = evt.getProductID< std::vector<recob::Vertex> >();
                auto kid = evt.getProductID< std::vector<recob::Vertex> >(kKinksName);
		auto const* kinkGetter = evt.productGetter(kid);

                auto tid = evt.getProductID< std::vector<recob::Track> >();
		auto const* trkGetter = evt.productGetter(tid);

		auto vsel = pmalgFitter.getVertices(fSaveOnlyBranchingVtx); // vtx pos's with vector of connected track idxs
		auto ksel = pmalgFitter.getKinks(); // pairs of kink position - associated track idx 
		std::map< size_t, art::Ptr<recob::Vertex> > frontVtxs; // front vertex ptr for each track index

		if (fPmaFitterConfig.RunVertexing()) // save vertices and vtx-trk assns
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
				if (vptr.isNull()) mf::LogWarning("PMAlgTrajFitter") << "Vertex ptr is null.";
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
				else mf::LogWarning("PMAlgTrajFitter") << "No tracks found at this vertex.";
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

		for (const auto & pfParticleItr : pfPartToTrackVecMap)
		{
			art::Ptr<recob::PFParticle> pfParticle(pfparticleHandle, pfParticleItr.first);
			mf::LogVerbatim("PMAlgTrajFitter") << "PFParticle key: " << pfParticle.key()
				<< ", self: " << pfParticle->Self() << ", #tracks: " << pfParticleItr.second.size();

			if (!pfParticle.isNull()) util::CreateAssn(*this, evt, pfParticle, pfParticleItr.second, *pfp2trk);
			else mf::LogError("PMAlgTrajFitter") << "Error in PFParticle lookup, pfparticle index: "
				<< pfParticleItr.first << ", key: " << pfParticle.key();
		}
	}

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(vtxs));
	evt.put(std::move(kinks), kKinksName);
	evt.put(std::move(nodes), kNodesName);

	evt.put(std::move(trk2hit_oldway)); // ****** REMEMBER to remove when FindMany improved ******
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));

	evt.put(std::move(sp2hit));
	evt.put(std::move(vtx2trk));
	evt.put(std::move(trk2kink), kKinksName);

	evt.put(std::move(pfp2trk));
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

DEFINE_ART_MODULE(PMAlgTrajFitter)

} // namespace trkf
