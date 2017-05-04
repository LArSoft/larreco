///////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTracking
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), June 2016
//
// Single track reconstruction toolkit based on Projection Matching Algorithm. Uses cluster collections
// to find single tracks (modularized version of our original code).
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgTracking_h
#define PMAlgTracking_h


// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlgCosmicTagger.h"
#include "larreco/RecoAlg/PMAlgVertexing.h"
#include "larreco/RecoAlg/PMAlgStitching.h"

// ROOT & C++
#include <memory>

namespace pma
{
	typedef std::map< size_t, pma::TrkCandidateColl > tpc_track_map;

	recob::Track convertFrom(const pma::Track3D& src, unsigned int tidx, int pdg = 0);

	class PMAlgTrackingBase;
	class PMAlgFitter;
	class PMAlgTracker;
}

class pma::PMAlgTrackingBase
{
public:

	const pma::TrkCandidateColl & result(void) { return fResult; }

	std::vector< std::pair< TVector3, std::vector< std::pair< size_t, bool > > > >
	getVertices(bool onlyBranching = false) const
	{ return fPMAlgVertexing.getVertices(fResult, onlyBranching); }

	std::vector< std::pair< TVector3, size_t > > getKinks(void) const
	{ return fPMAlgVertexing.getKinks(fResult); }

protected:

	PMAlgTrackingBase(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
		const pma::ProjectionMatchingAlg::Config& pmalgConfig,
		const pma::PMAlgVertexing::Config& pmvtxConfig);
	~PMAlgTrackingBase(void);

	void guideEndpoints(pma::TrkCandidateColl & tracks);

	pma::cryo_tpc_view_hitmap fHitMap;

	pma::ProjectionMatchingAlg fProjectionMatchingAlg;
	pma::PMAlgVertexing fPMAlgVertexing;

	pma::TrkCandidateColl fResult;
};


class pma::PMAlgFitter : public pma::PMAlgTrackingBase
{
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Sequence<int> TrackingOnlyPdg {
			Name("TrackingOnlyPdg"),
			Comment("PDG list to select which PFParticles should be reconstructed; all PFP's are used if the list is empty or starts with 0")
		};

		fhicl::Sequence<int> TrackingSkipPdg {
			Name("TrackingSkipPdg"),
			Comment("PDG list to select which PFParticles should NOT be reconstructed, e.g. skip EM-like if contains 11; no skipping if the list is empty or starts with 0")
		};

		fhicl::Atom<bool> RunVertexing {
			Name("RunVertexing"),
			Comment("find vertices from PFP hierarchy, join with tracks, reoptimize track-vertex structure")
		};
    };

	PMAlgFitter(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
		const std::vector< recob::Cluster > & clusters,
		const std::vector< recob::PFParticle > & pfparticles,
		const art::FindManyP< recob::Hit > & hitsFromClusters,
		const art::FindManyP< recob::Cluster > & clusFromPfps,
		const art::FindManyP< recob::Vertex > & vtxFromPfps,
		const pma::ProjectionMatchingAlg::Config& pmalgConfig,
		const pma::PMAlgFitter::Config& pmalgFitterConfig,
		const pma::PMAlgVertexing::Config& pmvtxConfig);

	int build(void);

private:

	void buildTracks(void);
	void buildShowers(void);

	bool has(const std::vector<int> & v, int i) const
	{
		for (auto c : v) { if (c == i) return true; }
		return false;
	}

	std::vector< std::vector< art::Ptr<recob::Hit> > > fCluHits;
	std::map< int, std::vector< art::Ptr<recob::Cluster> > > fPfpClusters;
	std::map< int, pma::Vector3D > fPfpVtx;
	std::map< int, int > fPfpPdgCodes;

	// ******************** fcl parameters ***********************
	std::vector<int> fTrackingOnlyPdg; // make tracks only for this pdg's when using input from PFParticles
	std::vector<int> fTrackingSkipPdg; // skip tracks with this pdg's when using input from PFParticles
	bool fRunVertexing;                // run vertex finding
};

class pma::PMAlgTracker : public pma::PMAlgTrackingBase
{
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Atom<size_t> MinSeedSize1stPass {
			Name("MinSeedSize1stPass"),
			Comment("min. cluster size used to start building a track in the 1st pass")
		};

		fhicl::Atom<size_t> MinSeedSize2ndPass {
			Name("MinSeedSize2ndPass"),
			Comment("min. cluster size used to start building a track in the 2nd pass")
		};

		fhicl::Atom<float> TrackLikeThreshold {
			Name("TrackLikeThreshold"),
			Comment("Threshold for track-like recognition")
		};

		fhicl::Atom<bool> RunVertexing {
			Name("RunVertexing"),
			Comment("find vertices from PFP hierarchy, join with tracks, reoptimize track-vertex structure")
		};

		fhicl::Atom<bool> FlipToBeam {
			Name("FlipToBeam"),
			Comment("set the track direction to increasing Z values")
		};

		fhicl::Atom<bool> FlipDownward {
			Name("FlipDownward"),
			Comment("set the track direction to decreasing Y values (like cosmic rays)")
		};

		fhicl::Atom<bool> AutoFlip_dQdx {
			Name("AutoFlip_dQdx"),
			Comment("set the track direction to increasing dQ/dx (overrides FlipToBeam and FlipDownward if significant rise of dQ/dx at the track end)")
		};

		fhicl::Atom<bool> MergeWithinTPC {
			Name("MergeWithinTPC"),
			Comment("merge witnin single TPC; finds tracks best matching by angle and displacement")
		};

		fhicl::Atom<double> MergeTransverseShift {
			Name("MergeTransverseShift"),
			Comment("max. transverse displacement [cm] between tracks")
		};

		fhicl::Atom<double> MergeAngle {
			Name("MergeAngle"),
			Comment("max. angle [degree] between tracks (nearest segments)")
		};

		fhicl::Atom<bool> StitchBetweenTPCs {
			Name("StitchBetweenTPCs"),
			Comment("stitch between TPCs; finds tracks best matching by angle and displacement")
		};

		fhicl::Atom<double> StitchDistToWall {
			Name("StitchDistToWall"),
			Comment("max. track endpoint distance [cm] to TPC boundary")
		};

		fhicl::Atom<double> StitchTransverseShift {
			Name("StitchTransverseShift"),
			Comment("max. transverse displacement [cm] between tracks")
		};

		fhicl::Atom<double> StitchAngle {
			Name("StitchAngle"),
			Comment("max. angle [degree] between tracks (nearest segments)")
		};

		fhicl::Atom<bool> MatchT0inAPACrossing {
			Name("MatchT0inAPACrossing"),
			Comment("match T0 of APA-crossing tracks using PMAlgStitcher")
		};

		fhicl::Atom<bool> MatchT0inCPACrossing {
			Name("MatchT0inCPACrossing"),
			Comment("match T0 of CPA-crossing tracks using PMAlgStitcher")
		};

  };

	PMAlgTracker(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
		const pma::ProjectionMatchingAlg::Config& pmalgConfig,
		const pma::PMAlgTracker::Config& pmalgTrackerConfig,
		const pma::PMAlgVertexing::Config& pmvtxConfig,
		const pma::PMAlgStitching::Config& pmstitchConfig,
		const pma::PMAlgCosmicTagger::Config& pmtaggerConfig) :

		PMAlgTrackingBase(allhitlist, pmalgConfig, pmvtxConfig),

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

        fGeom(&*(art::ServiceHandle<geo::Geometry>())),
		fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
	{}

	void init(const art::FindManyP< recob::Hit > & hitsFromClusters);

    void init(const art::FindManyP< recob::Hit > & hitsFromClusters,
        const std::vector< float > & trackLike);

	void init(const art::FindManyP< recob::Hit > & hitsFromClusters,
		const art::FindManyP< recob::Hit > & hitsFromEmParts);

	int build(void);

private:

	double collectSingleViewEnd(pma::Track3D & trk, std::vector< art::Ptr<recob::Hit> > & hits);
	double collectSingleViewFront(pma::Track3D & trk, std::vector< art::Ptr<recob::Hit> > & hits);

	bool reassignHits_1(const std::vector< art::Ptr<recob::Hit> > & hits,
		pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2);
	bool reassignSingleViewEnds_1(pma::TrkCandidateColl & tracks); // use clusters

	bool reassignHits_2(const std::vector< art::Ptr<recob::Hit> > & hits,
		pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2);
	bool reassignSingleViewEnds_2(pma::TrkCandidateColl & tracks);

	bool areCoLinear(pma::Track3D* trk1, pma::Track3D* trk2,
		double& dist, double& cos3d, bool& reverseOrder,
		double distThr, double distThrMin,
		double distProjThr,
		double cosThr);

	void freezeBranchingNodes(pma::TrkCandidateColl & tracks);
	void releaseAllNodes(pma::TrkCandidateColl & tracks);

	bool mergeCoLinear(pma::TrkCandidateColl & tracks);
	void mergeCoLinear(pma::tpc_track_map& tracks);

	double validate(pma::Track3D& trk, unsigned int testView);

	void fromMaxCluster_tpc(pma::TrkCandidateColl & result,
		size_t minBuildSize, unsigned int tpc, unsigned int cryo);

	pma::TrkCandidate matchCluster(int first_clu_idx, const std::vector< art::Ptr<recob::Hit> > & first_hits,
		size_t minSizeCompl, unsigned int tpc, unsigned int cryo, geo::View_t first_view);

	pma::TrkCandidate matchCluster(int first_clu_idx, size_t minSizeCompl,
		unsigned int tpc, unsigned int cryo, geo::View_t first_view)
	{
		return matchCluster(first_clu_idx, fCluHits[first_clu_idx], minSizeCompl, tpc, cryo, first_view);
	}

	int matchCluster(const pma::TrkCandidate& trk,
		size_t minSize, double fraction,
		unsigned int preferedView, unsigned int testView,
		unsigned int tpc, unsigned int cryo) const;

	bool extendTrack(pma::TrkCandidate& candidate,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		unsigned int testView, bool add_nodes);

	int maxCluster(int first_idx_tag,
		const pma::TrkCandidateColl & candidates,
		float xmin, float xmax, size_t min_clu_size,
		geo::View_t view, unsigned int tpc, unsigned int cryo) const;

	int maxCluster(size_t min_clu_size,
		geo::View_t view, unsigned int tpc, unsigned int cryo) const;

	void listUsedClusters(void) const;

	bool has(const std::vector<size_t>& v, size_t idx) const
	{
		for (auto c : v) if (c == idx) return true;
		return false;
	}

	std::vector< std::vector< art::Ptr<recob::Hit> > > fCluHits;
	std::vector< float > fCluWeights;

	/// --------------------------------------------------------------
	std::vector< size_t > fUsedClusters, fInitialClusters;
	mutable std::map< unsigned int, std::vector<size_t> > fTriedClusters;
	/// --------------------------------------------------------------

	// ******************** fcl parameters **********************
	size_t fMinSeedSize1stPass;  // min. cluster size used to start building a track in the 1st pass
	size_t fMinSeedSize2ndPass;  // min. cluster size used to start building a track in the 2nd pass
	float  fTrackLikeThreshold;  // trk-like threshold on cnn output
	double fMinTwoViewFraction;

	bool fFlipToBeam;            // set the track direction to increasing Z values
	bool fFlipDownward;          // set the track direction to decreasing Y values
	bool fAutoFlip_dQdx;         // set the track direction to increasing dQ/dx

	bool fMergeWithinTPC;          // merge witnin single TPC; finds tracks best matching by angle, with limits:
	double fMergeTransverseShift;  //   - max. transverse displacement [cm] between tracks
	double fMergeAngle;            //   - max. angle [degree] between tracks (nearest segments)

    pma::PMAlgCosmicTagger fCosmicTagger; // cosmic tagger alg
    bool fTagCosmicTracks;         // do any tagging of cosmic rays (simple or of tagger flags)

	bool fStitchBetweenTPCs;       // stitch between TPCs; finds tracks best matching by angle, with limits:
	double fStitchDistToWall;      //   - max. track endpoint distance [cm] to TPC boundary
	double fStitchTransverseShift; //   - max. transverse displacement [cm] between tracks
	double fStitchAngle;           //   - max. angle [degree] between tracks (nearest segments)

	bool fMatchT0inAPACrossing;    // match T0 of APA-crossing tracks using PMAlgStitcher
	bool fMatchT0inCPACrossing;    // match T0 of CPA-crossing tracks using PMAlgStitcher

    pma::PMAlgStitching fStitcher;

	bool fRunVertexing;          // run vertex finding

	// *********************** services *************************
	geo::GeometryCore const* fGeom;
	const detinfo::DetectorProperties* fDetProp;
};

#endif
