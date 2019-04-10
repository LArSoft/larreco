////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Authors      D.Stefan (Dorota.Stefan@ncbj.gov.pl),         from DUNE, CERN/NCBJ, since May 2015
//              R.Sulej (Robert.Sulej@cern.ch),               from DUNE, FNAL/NCBJ, since May 2015
//              L.Whitehead (leigh.howard.whitehead@cern.ch), from DUNE, CERN,      since Jan 2017
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
//    ~Jan-May 2017:   track stitching on APA and CPA, cosmics tagging
//    July 2017:       validation on 2D ADC image
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"

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
#include "lardataobj/AnalysisBase/T0.h" 
#include "lardataobj/AnalysisBase/CosmicTag.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlgTracking.h"
#include "larreco/RecoAlg/PMAlgVertexing.h"
#include "larreco/RecoAlg/PMAlgStitching.h"

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

        fhicl::Table<pma::PMAlgCosmicTagger::Config> PMAlgCosmicTagging {
            Name("PMAlgCosmicTagging")
        };

        fhicl::Table<pma::PMAlgStitching::Config> PMAlgStitching {
            Name("PMAlgStitching")
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

		fhicl::Atom< art::InputTag > WireModuleLabel {
			Name("WireModuleLabel"),
			Comment("tag of recob::Wire producer.")
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
    // will try to get EM- and track-like values from various lenght MVA vectors
    template <size_t N> bool init(const art::Event & evt, pma::PMAlgTracker & pmalgTracker) const;

    // calculate EM/track value for hits in track, in its best 2D projection
    // (tracks are built starting from track-like cluster, some electrons
    // still may look track-like)
    template <size_t N> int getPdgFromCnnOnHits(const art::Event& evt, const pma::Track3D& trk) const;

    // convert to a vector of LArSoft's cosmic tags
    std::vector<anab::CosmicTagID_t> getCosmicTag(const pma::Track3D::ETag pmaTag) const;

	// ******************** fcl parameters **********************
	art::InputTag fHitModuleLabel;  // tag for hits collection (used for trk validation)
	art::InputTag fWireModuleLabel; // tag for recob::Wire collection (used for trk validation)
	art::InputTag fCluModuleLabel;  // tag for input cluster collection
	art::InputTag fEmModuleLabel;   // tag for em-like cluster collection

	pma::ProjectionMatchingAlg::Config fPmaConfig;
	pma::PMAlgTracker::Config fPmaTrackerConfig;
	pma::PMAlgCosmicTagger::Config fPmaTaggingConfig;
	pma::PMAlgVertexing::Config fPmaVtxConfig;
    pma::PMAlgStitching::Config fPmaStitchConfig;

	bool fSaveOnlyBranchingVtx;  // for debugging, save only vertices which connect many tracks
	bool fSavePmaNodes;          // for debugging, save only track nodes

	// ********** instance names (collections, assns) ************
	static const std::string kKinksName;        // kinks on tracks
	static const std::string kNodesName;        // pma nodes

	// *********************** geometry **************************
	geo::GeometryCore const* fGeom;

    // histograms created only for the calibration of the ADC-based track validation mode
    std::vector< TH1F* > fAdcInPassingPoints, fAdcInRejectedPoints;
};
// -------------------------------------------------------------
const std::string PMAlgTrackMaker::kKinksName = "kink";
const std::string PMAlgTrackMaker::kNodesName = "node";
// -------------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(PMAlgTrackMaker::Parameters const& config) :
    EDProducer{config},
	fHitModuleLabel(config().HitModuleLabel()),
	fWireModuleLabel(config().WireModuleLabel()),
	fCluModuleLabel(config().ClusterModuleLabel()),
	fEmModuleLabel(config().EmClusterModuleLabel()),

	fPmaConfig(config().ProjectionMatchingAlg()),
	fPmaTrackerConfig(config().PMAlgTracking()),
	fPmaTaggingConfig(config().PMAlgCosmicTagging()),
	fPmaVtxConfig(config().PMAlgVertexing()),
    fPmaStitchConfig(config().PMAlgStitching()),

	fSaveOnlyBranchingVtx(config().SaveOnlyBranchingVtx()),
	fSavePmaNodes(config().SavePmaNodes()),
	
	fGeom( &*(art::ServiceHandle<geo::Geometry const>()))
{
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< std::vector<recob::Vertex> >(); // no instance name for interaction vertices
	produces< std::vector<recob::Vertex> >(kKinksName); // collection of kinks on tracks
	produces< std::vector<recob::Vertex> >(kNodesName); // collection of pma nodes
	produces< std::vector<anab::T0> >();
	produces< std::vector<anab::CosmicTag> >(); // Cosmic ray tags

	produces< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
	produces< art::Assns<recob::Vertex, recob::Track> >(); // no instance name for assns of tracks to interaction vertices
	produces< art::Assns<recob::Track, recob::Vertex> >(kKinksName);  // assns of kinks to tracks
	produces< art::Assns<recob::Track, anab::T0> >();
	produces< art::Assns<recob::Track, anab::CosmicTag> >(); // Cosmic ray tags associated to tracks

	produces< std::vector<recob::PFParticle> >();
	produces< art::Assns<recob::PFParticle, recob::Cluster> >();
	produces< art::Assns<recob::PFParticle, recob::Vertex> >();
	produces< art::Assns<recob::PFParticle, recob::Track> >();

    if (fPmaTrackerConfig.Validation() == "calib") // create histograms only in the calibration mode
    {
    	art::ServiceHandle<art::TFileService const> tfs;
	    for (size_t p = 0; p < fGeom->MaxPlanes(); ++p)
	    {
	        std::ostringstream ss1; ss1 << "adc_plane_" << p ;
	        fAdcInPassingPoints.push_back( tfs->make<TH1F>((ss1.str() + "_passing").c_str(), "max adc around the point on track", 100., 0., 5.) );
	        fAdcInRejectedPoints.push_back( tfs->make<TH1F>((ss1.str() + "_rejected").c_str(), "max adc around spurious point ", 100., 0., 5.) );
	    }
	}
}
// ------------------------------------------------------

template <size_t N>
int PMAlgTrackMaker::getPdgFromCnnOnHits(const art::Event& evt, const pma::Track3D& trk) const
{
    int pdg = 0;
    if (fPmaTrackerConfig.TrackLikeThreshold() > 0)
    {
        auto hitResults = anab::MVAReader<recob::Hit, N>::create(evt, fCluModuleLabel);
        if (hitResults)
        {
            int trkLikeIdx = hitResults->getIndex("track");
            int emLikeIdx = hitResults->getIndex("em");
            if ((trkLikeIdx < 0) || (emLikeIdx < 0))
            {
                throw cet::exception("PMAlgTrackMaker") << "No em/track labeled columns in MVA data products." << std::endl;
            }

            size_t nh[3] = { 0, 0, 0 };
            for (size_t hidx = 0; hidx < trk.size(); ++hidx) { ++nh[trk[hidx]->View2D()]; }

            size_t best_view = 2; // collection
            if ((nh[0] >= nh[1]) && (nh[0] > 1.25 * nh[2])) best_view = 0; // ind1
            if ((nh[1] >= nh[0]) && (nh[1] > 1.25 * nh[2])) best_view = 1; // ind2

            std::vector< art::Ptr<recob::Hit> > trkHitPtrList;
            trkHitPtrList.reserve(nh[best_view]);
            for (size_t hidx = 0; hidx < trk.size(); ++hidx)
            {
                if (trk[hidx]->View2D() == best_view) { trkHitPtrList.emplace_back(trk[hidx]->Hit2DPtr()); }
            }
            auto vout = hitResults->getOutput(trkHitPtrList);
            double trk_like = -1, trk_or_em = vout[trkLikeIdx] + vout[emLikeIdx];
            if (trk_or_em > 0)
            {
                trk_like = vout[trkLikeIdx] / trk_or_em;
                if (trk_like < fPmaTrackerConfig.TrackLikeThreshold()) pdg = 11; // tag if EM-like
                // (don't set pdg for track-like, for the moment don't like the idea of using "13")
            }
            //std::cout << "trk:" << best_view << ":" << trk.size() << ":" << trkHitPtrList.size() << " p=" << trk_like << std::endl;
        }
    }
    return pdg;
}

template <size_t N>
bool PMAlgTrackMaker::init(const art::Event & evt, pma::PMAlgTracker & pmalgTracker) const
{
    auto cluResults = anab::MVAReader< recob::Cluster, N >::create(evt, fCluModuleLabel);
    if (!cluResults) { return false; }

    int trkLikeIdx = cluResults->getIndex("track");
    int emLikeIdx = cluResults->getIndex("em");
    if ((trkLikeIdx < 0) || (emLikeIdx < 0)) { return false; }

    const art::FindManyP< recob::Hit > hitsFromClusters(cluResults->dataHandle(), evt, cluResults->dataTag());
    const auto & cnnOuts = cluResults->outputs();
    std::vector< float > trackLike(cnnOuts.size());
    for (size_t i = 0; i < cnnOuts.size(); ++i)
    {
        double trkOrEm = cnnOuts[i][trkLikeIdx] + cnnOuts[i][emLikeIdx];
        if (trkOrEm > 0) { trackLike[i] = cnnOuts[i][trkLikeIdx] / trkOrEm; }
        else { trackLike[i] = 0; }
    }
    pmalgTracker.init(hitsFromClusters, trackLike);
    return true;
}

std::vector<anab::CosmicTagID_t> PMAlgTrackMaker::getCosmicTag(const pma::Track3D::ETag pmaTag) const
{
    std::vector<anab::CosmicTagID_t> anabTags;
    
    if (pmaTag & pma::Track3D::kOutsideDrift_Partial)  { anabTags.push_back(anab::CosmicTagID_t::kOutsideDrift_Partial);   }
    if (pmaTag & pma::Track3D::kOutsideDrift_Complete) { anabTags.push_back(anab::CosmicTagID_t::kOutsideDrift_Complete);  }
    if (pmaTag & pma::Track3D::kBeamIncompatible)      { anabTags.push_back(anab::CosmicTagID_t::kFlash_BeamIncompatible); }
    if (pmaTag & pma::Track3D::kGeometry_XX)           { anabTags.push_back(anab::CosmicTagID_t::kGeometry_XX);            }
    if (pmaTag & pma::Track3D::kGeometry_YY)           { anabTags.push_back(anab::CosmicTagID_t::kGeometry_YY);            }
    if (pmaTag & pma::Track3D::kGeometry_ZZ)           { anabTags.push_back(anab::CosmicTagID_t::kGeometry_ZZ);            }
    if (pmaTag & pma::Track3D::kGeometry_YZ)           { anabTags.push_back(anab::CosmicTagID_t::kGeometry_YZ);            }
    if (pmaTag & pma::Track3D::kGeometry_Y)            { anabTags.push_back(anab::CosmicTagID_t::kGeometry_Y);             }

    if (anabTags.empty())                              { anabTags.push_back(anab::CosmicTagID_t::kUnknown);                }
    
    return anabTags;
}

void PMAlgTrackMaker::produce(art::Event& evt)
{
	// ---------------- Create data products --------------------------
	auto tracks = std::make_unique< std::vector< recob::Track > >();
	auto allsp = std::make_unique< std::vector< recob::SpacePoint > >();
	auto vtxs = std::make_unique< std::vector< recob::Vertex > >();  // interaction vertices
	auto kinks = std::make_unique< std::vector< recob::Vertex > >(); // kinks on tracks (no new particles start in kinks)
	auto nodes = std::make_unique< std::vector< recob::Vertex > >(); // pma nodes
	auto t0s = std::make_unique< std::vector< anab::T0 > >();
	auto cosmicTags = std::make_unique< std::vector< anab::CosmicTag > >();

	auto trk2hit_oldway = std::make_unique< art::Assns< recob::Track, recob::Hit > >(); // ****** REMEMBER to remove when FindMany improved ******
	auto trk2hit = std::make_unique< art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta > >();

	auto trk2sp = std::make_unique< art::Assns< recob::Track, recob::SpacePoint > >();
	auto trk2t0 = std::make_unique< art::Assns< recob::Track, anab::T0 > >();
	auto trk2ct = std::make_unique< art::Assns< recob::Track, anab::CosmicTag > >();

	auto sp2hit = std::make_unique< art::Assns< recob::SpacePoint, recob::Hit > >();
	auto vtx2trk = std::make_unique< art::Assns< recob::Vertex, recob::Track > >();  // one or more tracks (particles) start in the vertex
	auto trk2kink = std::make_unique< art::Assns< recob::Track, recob::Vertex > >(); // one or more kinks on the track

	auto pfps = std::make_unique< std::vector< recob::PFParticle > >();

	auto pfp2clu = std::make_unique< art::Assns<recob::PFParticle, recob::Cluster> >();
	auto pfp2vtx = std::make_unique< art::Assns<recob::PFParticle, recob::Vertex> >();
	auto pfp2trk = std::make_unique< art::Assns< recob::PFParticle, recob::Track > >();


	// --------------------- Wires & Hits -----------------------------
	auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireModuleLabel);
	auto allHitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	std::vector< art::Ptr<recob::Hit> > allhitlist;
	art::fill_ptr_vector(allhitlist, allHitListHandle);


	// -------------- PMA Tracker for this event ----------------------
	auto pmalgTracker = pma::PMAlgTracker(allhitlist, *wireHandle,
		fPmaConfig, fPmaTrackerConfig, fPmaVtxConfig, fPmaStitchConfig, fPmaTaggingConfig, fAdcInPassingPoints, fAdcInRejectedPoints);

    size_t mvaLength = 0;
	if (fEmModuleLabel != "") // ----------- Exclude EM parts ---------
	{
	    auto cluListHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fCluModuleLabel);
	    auto splitCluHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fEmModuleLabel);

	    art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fCluModuleLabel);
		art::FindManyP< recob::Hit > hitsFromEmParts(splitCluHandle, evt, fEmModuleLabel);
		pmalgTracker.init(hitsFromClusters, hitsFromEmParts);
	}
	else if (fPmaTrackerConfig.TrackLikeThreshold() > 0) // --- CNN EM/trk separation ----
	{
	    // try to dig out 4- or 3-output MVA data product
	    if (init<4>(evt, pmalgTracker) )      { mvaLength = 4; } // e.g.: EM / track / Michel / none
	    else if (init<3>(evt, pmalgTracker))  { mvaLength = 3; } // e.g.: EM / track / none
	    else if (init<2>(evt, pmalgTracker))  { mvaLength = 2; } // just EM / track (LArIAT starts with this style)
	    else
	    {
	        throw cet::exception("PMAlgTrackMaker") << "No EM/track MVA data products." << std::endl;
	    }
	}
	else // ------------------------ Use ALL clusters -----------------
	{
	    auto cluListHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fCluModuleLabel);
	    art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fCluModuleLabel);
		pmalgTracker.init(hitsFromClusters);
	}


	// ------------------ Do the job here: ----------------------------
	int retCode = pmalgTracker.build();
	// ----------------------------------------------------------------
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

	// ---------- Translate output to data products: ------------------
	auto const & result = pmalgTracker.result();
	if (!result.empty()) // ok, there is something to save
	{
		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6];
		for (size_t i = 0; i < 3; ++i) sp_pos[i] = 0.0;
		for (size_t i = 0; i < 6; ++i) sp_err[i] = 1.0;

        auto const make_pfpptr = art::PtrMaker<recob::PFParticle>(evt);
                auto const make_trkptr = art::PtrMaker<recob::Track>(evt); // PtrMaker Step #1
                auto const make_vtxptr = art::PtrMaker<recob::Vertex>(evt);
                auto const make_kinkptr = art::PtrMaker<recob::Vertex>(evt, kKinksName);
                auto const make_t0ptr = art::PtrMaker<anab::T0>(evt);
                auto const make_ctptr = art::PtrMaker<anab::CosmicTag>(evt);

		tracks->reserve(result.size());
		for (size_t trkIndex = 0; trkIndex < result.size(); ++trkIndex)
		{
			pma::Track3D* trk = result[trkIndex].Track();

			trk->SelectHits();  // just in case, set all to enabled
			unsigned int itpc = trk->FrontTPC(), icryo = trk->FrontCryo();
			pma::dedx_map dedx_tmp;
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kU)) { trk->CompleteMissingWires(geo::kU); trk->GetRawdEdxSequence(dedx_tmp, geo::kU, 1); }
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kV)) { trk->CompleteMissingWires(geo::kV); trk->GetRawdEdxSequence(dedx_tmp, geo::kV, 1); }
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kX)) { trk->CompleteMissingWires(geo::kX); trk->GetRawdEdxSequence(dedx_tmp, geo::kX, 1); }
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kY)) { trk->CompleteMissingWires(geo::kY); trk->GetRawdEdxSequence(dedx_tmp, geo::kY, 1); }
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kZ)) { trk->CompleteMissingWires(geo::kZ); trk->GetRawdEdxSequence(dedx_tmp, geo::kZ, 1); }

			//gc: make sure no tracks are created with less than 2 points
			if (trk->size()<2) continue;

		    int pdg = 0;
		    if (mvaLength == 4) pdg = getPdgFromCnnOnHits<4>(evt, *(result[trkIndex].Track()));
		    else if (mvaLength == 3) pdg = getPdgFromCnnOnHits<3>(evt, *(result[trkIndex].Track()));
		    else if (mvaLength == 2) pdg = getPdgFromCnnOnHits<2>(evt, *(result[trkIndex].Track()));
		    //else mf::LogInfo("PMAlgTrackMaker") << "Not using PID from CNN.";

			tracks->push_back(pma::convertFrom(*trk, trkIndex, pdg));

			auto const trkPtr = make_trkptr(tracks->size() - 1); // PtrMaker Step #2

			if (trk->HasT0())
			{
				// TriggBits=3 means from 3d reco (0,1,2 mean something else)
				t0s->push_back(anab::T0(trk->GetT0(), 0, 3, tracks->back().ID()));

				auto const t0Ptr = make_t0ptr(t0s->size() - 1);  // PtrMaker Step #3
				trk2t0->addSingle(trkPtr, t0Ptr);
			}

			// Check if this is a cosmic ray and create an association if it is.
			if(trk->HasTagFlag(pma::Track3D::kCosmic)){
				// Get the track end points
				std::vector<float> trkEnd0;
				std::vector<float> trkEnd1;
				// Get the drift direction, but don't care about the sign
				// Also need to subtract 1 due to the definition.
				int driftDir = abs(fGeom->TPC(trk->FrontTPC(), trk->FrontCryo()).DetectDriftDirection()) - 1;
				
				for(int i = 0; i < 3; ++i){
					// Get the drift direction and apply the opposite of the drift shift in order to
					// give the CosmicTag the drift coordinate assuming T0 = T_beam as it requests. 
					double shift = 0.0;
					if(i == driftDir){
						shift = trk->Nodes()[0]->GetDriftShift();
					}
					trkEnd0.push_back(trk->Nodes()[0]->Point3D()[i] - shift);
					trkEnd1.push_back(trk->Nodes()[trk->Nodes().size()-1]->Point3D()[i] - shift);
				}
				// Make the tag object. For now, let's say this is very likely a cosmic (3rd argument = 1).
				// Set the type of cosmic to the value saved in pma::Track.
				auto tags = getCosmicTag(trk->GetTag());
				for (const auto t : tags)
				{
					cosmicTags->emplace_back(trkEnd0, trkEnd1, 1, t);
					auto const cosmicPtr = make_ctptr(cosmicTags->size()-1);
					trk2ct->addSingle(trkPtr,cosmicPtr);
				}
			}

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
			//rs: so make also associations to space points in the order of trajectory points
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
		}

		auto vsel = pmalgTracker.getVertices(fSaveOnlyBranchingVtx); // vtx pos's with vector of connected track idxs
		auto ksel = pmalgTracker.getKinks(); // pairs of kink position - associated track idx 
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

                auto const vptr = make_vtxptr(vidx);
				if (!v.second.empty())
				{
					for (const auto & vEntry : v.second)
					{
						size_t tidx = vEntry.first;
						bool isFront = vEntry.second;

						if (isFront) frontVtxs[tidx] = vptr; // keep ptr of the front vtx

						auto const tptr = make_trkptr(tidx);
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

				auto const tptr = make_trkptr(tidx);
				auto const kptr = make_kinkptr(kidx);
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

		for (size_t t = 0; t < result.size(); ++t)
		{
			size_t parentIdx = recob::PFParticle::kPFParticlePrimary;
			if (result[t].Parent() >= 0) parentIdx = (size_t)result[t].Parent();

			std::vector< size_t > daughterIdxs;
			for (size_t idx : result[t].Daughters()) { daughterIdxs.push_back(idx); }

			size_t pfpidx = pfps->size();
			pfps->emplace_back((*tracks)[t].ParticleId(), pfpidx, parentIdx, daughterIdxs);

			auto const pfpptr = make_pfpptr(pfpidx);
			auto const tptr = make_trkptr(t);
			pfp2trk->addSingle(pfpptr, tptr);

			// add assns to FRONT vertex of each particle
			if (fPmaTrackerConfig.RunVertexing())
			{
				art::Ptr<recob::Vertex> vptr = frontVtxs[t];
				if (!vptr.isNull()) pfp2vtx->addSingle(pfpptr, vptr);
				else mf::LogWarning("PMAlgTrackMaker") << "Front vertex for PFParticle is missing.";
			}
		}
		mf::LogVerbatim("Summary") << pfps->size() << " PFParticles created for reconstructed tracks.";
        mf::LogVerbatim("Summary") << "Adding " << result.parents().size() << " primary PFParticles.";
		for (size_t t = 0; t < result.parents().size(); ++t)
		{
			std::vector< size_t > daughterIdxs;
			for (size_t idx : result.parents()[t].Daughters()) { daughterIdxs.push_back(idx); }

			size_t pfpidx = pfps->size();
			size_t parentIdx = recob::PFParticle::kPFParticlePrimary;
			pfps->emplace_back(0, pfpidx, parentIdx, daughterIdxs);

			// add assns to END vertex of primary
			if (fPmaTrackerConfig.RunVertexing() && !daughterIdxs.empty())
			{
			    auto const pfpptr = make_pfpptr(pfpidx);
				art::Ptr<recob::Vertex> vptr = frontVtxs[daughterIdxs.front()]; // same vertex for all daughters
				if (!vptr.isNull()) pfp2vtx->addSingle(pfpptr, vptr);
				else mf::LogWarning("PMAlgTrackMaker") << "Front vertex for PFParticle is missing.";
			}
		}
		mf::LogVerbatim("Summary") << pfps->size() << " PFParticles created in total.";
	}

    // for (const auto & ct : *cosmicTags) { std::cout << "Cosmic tag: " << ct << std::endl; }

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(vtxs));
	evt.put(std::move(kinks), kKinksName);
	evt.put(std::move(nodes), kNodesName);
	evt.put(std::move(t0s));
	evt.put(std::move(cosmicTags));

	evt.put(std::move(trk2hit_oldway)); // ****** REMEMBER to remove when FindMany improved ******
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(trk2t0));
	evt.put(std::move(trk2ct));

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
