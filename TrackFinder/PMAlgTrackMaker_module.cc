////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks using Projection Matching Algorithm, see RecoAlg/ProjectionMatchingAlg.h
// for basics of the PMA algorithm and its possible settings.
//
// Progress:
//    May-June 2015:   track finding and validation, growing tracks by iterative merging of matching
//                     clusters, no attempts to build multi-track structures, however cosmic tracking
//                     works fine as they are sets of independent tracks
//    June-July 2015:  merging track parts within a single tpc and stitching tracks across tpc's
//    August 2015:     optimization of track-vertex structures (so 3D vertices are also produced)
//    November 2015:   use track-shower splitting at 2D level, then tag low-E EM cascades in 3D
//                     note: the splitter is not finished and not as good as we want it
//    January 2016:    output results of track-vertex finding as a tree of PFParticles, also refined
//                     vertexing code
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Track.h"
#include "RecoBase/TrackHitMeta.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/SpacePoint.h"
#include "AnalysisBase/T0.h" 
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "MCCheater/BackTracker.h"
#include "Simulation/ParticleList.h"
#include "SimulationBase/MCParticle.h"

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlgVertexing.h"
#include "RecoAlg/PMAlg/Utilities.h"
#include "RecoAlg/PMAlg/PmaTrkCandidate.h"

#include "TTree.h"
#include "TMath.h"

#include <memory>

namespace trkf {

typedef std::map< size_t, std::vector<double> > dedx_map;
typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

typedef std::map< size_t, pma::trk_candidates > tpc_track_map;

class PMAlgTrackMaker : public art::EDProducer {
public:
  explicit PMAlgTrackMaker(fhicl::ParameterSet const & p);

  PMAlgTrackMaker(PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker(PMAlgTrackMaker &&) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker &&) = delete;

  void beginJob() override;

  void reconfigure(fhicl::ParameterSet const& p) override;

  void produce(art::Event & e) override;


private:
  // *** methods to create tracks from clusters ***

  // loop over all clusters and build as much as possible to find
  // the logic implemented here for sure is not exhaustive, it was
  // checked on long cosmic tracks and low energy stopping tracks
  // and seems to be a good example how to use the algorithm
  int fromMaxCluster(const art::Event& evt, pma::trk_candidates& result);

  void fromMaxCluster_tpc(
    pma::trk_candidates& result,
    const std::vector< art::Ptr<recob::Cluster> >& clusters,
    size_t minBuildSize, unsigned int tpc, unsigned int cryo,
	int pfParticleIdx = -1);

  bool extendTrack(
	pma::TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes);

  int matchCluster(
    const pma::TrkCandidate& trk,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo);

  // display what was used and what is left
  void listUsedClusters(const std::vector< art::Ptr<recob::Cluster> >& clusters) const;
  // ------------------------------------------------------

  // build tracks from clusters associated to PFParticles (use internal pattern recognition
  // on the subset of clusters selected with PFParticle)
  int fromPfpClusterSubset(const art::Event& evt, pma::trk_candidates& result);
  // ------------------------------------------------------

  // build tracks from straight from clusters associated to PFParticle (no pattern recognition)
  int fromPfpDirect(const art::Event& evt, pma::trk_candidates& result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  void reset(const art::Event& evt);

  cryo_tpc_view_hitmap fHitMap;
  std::vector< std::vector< art::Ptr<recob::Hit> > > fCluHits;
  std::map< int, std::vector< art::Ptr<recob::Cluster> > > fPfpClusters;
  bool sortHits(const art::Event& evt);
  bool sortHitsPfp(const art::Event& evt);

  std::vector< size_t > used_clusters, initial_clusters;
  std::map< unsigned int, std::vector<size_t> > tried_clusters;
  bool has(const std::vector<size_t>& v, size_t idx) const
  {
  	for (auto c : v) if (c == idx) return true;
  	return false;
  }

  // temporary set of possible solutions of the selected cluster and clusters in complementary views
  pma::trk_candidates fCandidates;

  int maxCluster(
    const std::vector< art::Ptr<recob::Cluster> >& clusters,
    size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);
  int maxCluster(
	size_t first_idx,
    const std::vector< art::Ptr<recob::Cluster> >& clusters,
    float tmin, float tmax, size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);

  void freezeBranchingNodes(pma::trk_candidates& tracks);
  void releaseAllNodes(pma::trk_candidates& tracks);

  bool areCoLinear(
	pma::Track3D* trk1, pma::Track3D* trk2,
	double& dist, double& cos, bool& reverseOrder,
	double distThr, double distThrMin,
	double distProjThr,
	double cosThr);
  void mergeCoLinear(pma::trk_candidates& tracks);
  void mergeCoLinear(tpc_track_map& tracks);

  bool areCoLinear(
		double& cos3d,
		TVector3 f0, TVector3 b0, TVector3 f1, TVector3 b1,
		double distProjThr);
  void matchCoLinearAnyT0(pma::trk_candidates& tracks);

  bool reassignHits(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	pma::trk_candidates& tracks, size_t trk_idx);
  bool reassignSingleViewEnds(pma::trk_candidates& tracks);
  void guideEndpoints(pma::trk_candidates& tracks);

  double validate(pma::Track3D& trk, unsigned int testView);
  recob::Track convertFrom(const pma::Track3D& src);

  bool fIsRealData;
  bool isMcStopping(void) const; // to be moved to the testing module
  int getMcPdg(const pma::Track3D& trk) const; // to be moved to the testing module
  // ------------------------------------------------------

  art::ServiceHandle< geo::Geometry > fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

  // ******************* tree output **********************
  int fEvNumber;        // event number
  int fTrkIndex;        // track index in the event, same for all dQ/dx points of the track
  int fPlaneIndex;      // wire plane index of the dQ/dx data point
  int fPlaneNPoints;    // number of data points in this plane
  int fIsStopping;      // tag tracks of stopping particles
  int fMcPdg;           // MC truth PDG matched for a track
  int fPidTag;          // Tag: 0=trk-like, 1=cascade-like
  double fdQdx;         // dQ/dx data point stored for each plane
  double fRange;        // residual range at dQ/dx data point, tracks are auto-flipped
  double fLength;       // track length
  double fHitsMse;      // MSE of hits: mean dist^2 of hit to 2D track projection
  double fSegAngMean;   // Mean segment-segment 3D angle.
  TTree* fTree_dQdx;    // dQ/dx info (saved optionally)
  TTree* fTree_trk;     // overall info

  // ******************** parameters **********************
  std::string fHitModuleLabel; // label for hits collection (used for trk validation)
  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association

  bool fMakePFPs;              // output track-vertex net as a tree of PFParticles

  size_t fMinSeedSize1stPass;  // min. cluster size used to start building a track in the 1st pass
  size_t fMinSeedSize2ndPass;  // min. cluster size used to start building a track in the 2nd pass

  bool fFlipToBeam;            // set the track direction to increasing Z values
  bool fFlipDownward;          // set the track direction to decreasing Y values
  bool fAutoFlip_dQdx;         // set the track direction to increasing dQ/dx
  bool fSave_dQdx;             // for debugging purposes, off by default

  bool fMergeWithinTPC;          // merge witnin single TPC; finds tracks best matching by angle, with limits:
  double fMergeTransverseShift;  //   - max. transverse displacement [cm] between tracks
  double fMergeAngle;            //   - max. angle [degree] between tracks (nearest segments)

  bool fStitchBetweenTPCs;       // stitch between TPCs; finds tracks best matching by angle, with limits:
  double fStitchDistToWall;      //   - max. track endpoint distance [cm] to TPC boundary
  double fStitchTransverseShift; //   - max. transverse displacement [cm] between tracks
  double fStitchAngle;           //   - max. angle [degree] between tracks (nearest segments)

  bool fMatchT0inAPACrossing;    // match T0 of APA-crossing tracks, TPC stitching limits are used

  pma::ProjectionMatchingAlg fProjectionMatchingAlg;
  double fMinTwoViewFraction;    // ProjectionMatchingAlg parameter used also in the module

  pma::PMAlgVertexing fPMAlgVertexing;
  bool fRunVertexing;          // run vertex finding
};
// ------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(fhicl::ParameterSet const & p) :
	fProjectionMatchingAlg(p.get< fhicl::ParameterSet >("ProjectionMatchingAlg")),
	fPMAlgVertexing(p.get< fhicl::ParameterSet >("PMAlgVertexing"))
{
	this->reconfigure(p);

	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< std::vector<recob::Vertex> >();
	produces< std::vector<anab::T0> >();

	produces< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
	produces< art::Assns<recob::Vertex, recob::Track> >();
	produces< art::Assns<recob::Track, anab::T0> >();

	produces< std::vector<recob::PFParticle> >();
	produces< art::Assns<recob::PFParticle, recob::Cluster> >();
	produces< art::Assns<recob::PFParticle, recob::Vertex> >();
	produces< art::Assns<recob::PFParticle, recob::Track> >();
}
// ------------------------------------------------------

void PMAlgTrackMaker::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree_dQdx = tfs->make<TTree>("PMAlgTrackMaker_dQdx", "tracks dQ/dx info");
	fTree_dQdx->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fTree_dQdx->Branch("fTrkIndex", &fTrkIndex, "fTrkIndex/I");
	fTree_dQdx->Branch("fPlaneIndex", &fPlaneIndex, "fPlaneIndex/I");
	fTree_dQdx->Branch("fPlaneNPoints", &fPlaneNPoints, "fPlaneNPoints/I");
	fTree_dQdx->Branch("fIsStopping", &fIsStopping, "fIsStopping/I");
	fTree_dQdx->Branch("fMcPdg", &fMcPdg, "fMcPdg/I");
	fTree_dQdx->Branch("fdQdx", &fdQdx, "fdQdx/D");
	fTree_dQdx->Branch("fRange", &fRange, "fRange/D");
	fTree_dQdx->Branch("fLength", &fLength, "fLength/D");

	fTree_trk = tfs->make<TTree>("PMAlgTrackMaker_trk", "tracks overall info");
	fTree_trk->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fTree_trk->Branch("fTrkIndex", &fTrkIndex, "fTrkIndex/I");
	fTree_trk->Branch("fMcPdg", &fMcPdg, "fMcPdg/I");
	fTree_trk->Branch("fLength", &fLength, "fLength/D");
	fTree_trk->Branch("fHitsMse", &fHitsMse, "fHitsMse/D");
	fTree_trk->Branch("fSegAngMean", &fSegAngMean, "fSegAngMean/D");
	fTree_trk->Branch("fPidTag", &fPidTag, "fPidTag/I");
}

void PMAlgTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
	fCluMatchingAlg = pset.get< int >("CluMatchingAlg");

	fMakePFPs = pset.get< bool >("MakePFPs");

	fMinSeedSize1stPass = pset.get< size_t >("MinSeedSize1stPass");
	fMinSeedSize2ndPass = pset.get< size_t >("MinSeedSize2ndPass");

	fFlipToBeam = pset.get< bool >("FlipToBeam");
	fFlipDownward = pset.get< bool >("FlipDownward");
	fAutoFlip_dQdx = pset.get< bool >("AutoFlip_dQdx");
	fSave_dQdx = pset.get< bool >("Save_dQdx");

	fMergeWithinTPC = pset.get< bool >("MergeWithinTPC");
	fMergeTransverseShift = pset.get< double >("MergeTransverseShift");
	fMergeAngle = pset.get< double >("MergeAngle");

	fStitchBetweenTPCs = pset.get< bool >("StitchBetweenTPCs");
	fStitchDistToWall = pset.get< double >("StitchDistToWall");
	fStitchTransverseShift = pset.get< double >("StitchTransverseShift");
	fStitchAngle = pset.get< double >("StitchAngle");

	fMatchT0inAPACrossing = pset.get< bool >("MatchT0inAPACrossing");

	fProjectionMatchingAlg.reconfigure(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"));
	fMinTwoViewFraction = pset.get< double >("ProjectionMatchingAlg.MinTwoViewFraction");

	fPMAlgVertexing.reconfigure(pset.get< fhicl::ParameterSet >("PMAlgVertexing"));
	fRunVertexing = pset.get< bool >("RunVertexing");
}

void PMAlgTrackMaker::reset(const art::Event& evt)
{
	fHitMap.clear();
	fCluHits.clear();
	fPfpClusters.clear();
	fEvNumber = evt.id().event();
	fIsRealData = evt.isRealData();
	fTrkIndex = 0;
	fPlaneIndex = 0;
	fPlaneNPoints = 0;
	fIsStopping = 0;
	fMcPdg = 0;
	fPidTag = 0;
	fdQdx = 0.0;
	fRange = 0.0;
	fLength = 0.0;
	fHitsMse = 0.0;
	fSegAngMean = 0.0;

	fPMAlgVertexing.reset();
}
// ------------------------------------------------------

recob::Track PMAlgTrackMaker::convertFrom(const pma::Track3D& src)
{
	std::vector< TVector3 > xyz, dircos;
	xyz.reserve(src.size()); dircos.reserve(src.size());

	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = src.FrontCryo();
	unsigned int tpc = src.FrontTPC();

	std::map< unsigned int, dedx_map > src_dQdx;
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		src_dQdx[geo::kU] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kU], geo::kU);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		src_dQdx[geo::kV] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kV], geo::kV);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		src_dQdx[geo::kZ] = dedx_map();
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
			TVector3 dc(src[i + 1]->Point3D());
			dc -= src[i]->Point3D();
			dc *= 1.0 / dc.Mag();
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

	fLength = src.Length();
	fHitsMse = src.GetMse();
	fSegAngMean = src.GetMeanAng();

	// 0 is track-like (long and/or very straight, well matching 2D hits);
	// 0x10000 is EM shower-like trajectory
	if (src.GetTag() == pma::Track3D::kEmLike) fPidTag = 0x10000;
	else fPidTag = 0;

	fMcPdg = getMcPdg(src);

	if (fSave_dQdx)
	{
		fIsStopping = (int)isMcStopping();
		for (unsigned int view = 0; view < fGeom->Nviews(); view++)
			if (fGeom->TPC(tpc, cryo).HasPlane(view))
			{
				fPlaneIndex = view;
				fPlaneNPoints = src_dQdx[view].size();
				for (auto const& data : src_dQdx[view])
				{
					double dQ = data.second[5];
					double dx = data.second[6];
					fRange = data.second[7];
					if (dx > 0.)
					{
						fdQdx = dQ/dx;
						fTree_dQdx->Fill();
					}
				}
			}
	}
	fTree_trk->Fill();

	if (xyz.size() != dircos.size())
	{
		mf::LogError("PMAlgTrackMaker") << "pma::Track3D to recob::Track conversion problem.";
	}

	return recob::Track(xyz, dircos, dst_dQdx, std::vector< double >(2, util::kBogusD), fTrkIndex + fPidTag);
}
// ------------------------------------------------------

bool PMAlgTrackMaker::isMcStopping(void) const
{
	if (fIsRealData)
	{
		// maybe possible to apply check here
		return false;
	}
	else
	{
		art::ServiceHandle< cheat::BackTracker > bt;
		const sim::ParticleList& plist = bt->ParticleList();
		const simb::MCParticle* particle = plist.Primary(0);

		if (particle)
		{
		//	std::cout << "...:SIM:... " << particle->EndProcess()
		//		<< " n:" << particle->NumberDaughters()
		//		<< " m:" << particle->Mass()
		//		<< " E:" << particle->EndE() << std::endl;
			return (particle->EndE() - particle->Mass() < 0.001);
		//	return (particle->NumberDaughters() == 0);
		}
		else return false;
	}
}

int PMAlgTrackMaker::getMcPdg(const pma::Track3D& trk) const
{
	if (fIsRealData) return 0;

	art::ServiceHandle< cheat::BackTracker > bt;
	std::map< int, size_t > pdg_counts;
	for (size_t i = 0; i < trk.size(); i++)
		if (trk[i]->IsEnabled() && (trk[i]->View2D() == geo::kZ))
	{
		art::Ptr< recob::Hit > hit = trk[i]->Hit2DPtr();
		std::vector< sim::TrackIDE > tids = bt->HitToTrackID(hit);
		for (auto const& tid : tids)
		{
			const simb::MCParticle* pi = bt->TrackIDToParticle(tid.trackID);
			int pdg = pi->PdgCode();
			auto it = pdg_counts.find(pdg);
			if (it != pdg_counts.end()) pdg_counts[pdg]++;
			else pdg_counts[pdg] = 1;
        }
	}

	int best_pdg = 0;
	size_t count, max_count = 0;
	for (auto const& p : pdg_counts)
	{
		count = p.second;
		if (count > max_count)
		{
			max_count = count; best_pdg = p.first;
		}
	}
	return best_pdg;
}
// ------------------------------------------------------

double PMAlgTrackMaker::validate(pma::Track3D& trk, unsigned int testView)
{
	if ((trk.FirstElement()->GetDistToWall() < -3.0) ||
	    (trk.LastElement()->GetDistToWall() < -3.0))
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "first or last node too far out of its initial TPC";
		return 0.0;
	}

	if (testView != geo::kUnknown)
		mf::LogVerbatim("PMAlgTrackMaker") << "validation in plane: " << testView;
	else return 1.0;

	std::vector< art::Ptr<recob::Hit> >& hits = fHitMap[trk.FrontCryo()][trk.FrontTPC()][testView];

	// always validate (needed for disambiguation postponed to 3D step):
	return fProjectionMatchingAlg.validate(trk, hits, testView);

	// in case of usual disambig on hit level one may validate only if there are at least a few hits:
/*	if (hits.size() > 10) return fProjectionMatchingAlg.validate(trk, hits, testView);
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "   too few hits (" << hits.size() << ")";
		return 1.0;
	}
*/
}
// ------------------------------------------------------

bool PMAlgTrackMaker::extendTrack(pma::TrkCandidate& candidate,
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
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track EXTENDED, MSE = " << m1 << ", v = " << v1;
		candidate.SetTrack(copy);  // replace with the new track (deletes old one)
		copy->SortHits();          // sort hits in the new track

		candidate.SetMse(m1);      // save info
		candidate.SetValidation(v1);

		return true;
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track NOT extended, MSE = " << m1 << ", v = " << v1;
		delete copy;
		return false;
	}
}
// ------------------------------------------------------

void PMAlgTrackMaker::freezeBranchingNodes(pma::trk_candidates& tracks)
{
	for (auto const & trk : tracks)
		for (auto node : trk.Track()->Nodes())
			if (node->IsBranching()) node->SetFrozen(true);
}
// ------------------------------------------------------

void PMAlgTrackMaker::releaseAllNodes(pma::trk_candidates& tracks)
{
	for (auto const & trk : tracks)
		for (auto node : trk.Track()->Nodes())
			node->SetFrozen(false);
}
// ------------------------------------------------------

bool PMAlgTrackMaker::areCoLinear(pma::Track3D* trk1, pma::Track3D* trk2,
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
			default: mf::LogError("PMAlgTrackMaker") << "Should never happen.";
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

		TVector3 dir1 = trk1->Segments().back()->GetDirection3D();
		TVector3 dir2 = trk2->Segments().front()->GetDirection3D();

		cos3d = dir1 * dir2;

		if ((cos3d > cosThr) && (distProj1 < distProjThr) && (distProj2 < distProjThr))
			return true;
		else // check if parallel to wires & colinear in 2D
		{
			const double maxCosXZ = 0.996195; // 5 deg

			TVector3 dir1_xz(dir1.X(), 0., dir1.Z());
			dir1_xz *= 1.0 / dir1_xz.Mag();

			TVector3 dir2_xz(dir2.X(), 0., dir2.Z());
			dir2_xz *= 1.0 / dir2_xz.Mag();

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
				double cosXZ = dir1_xz * dir2_xz;
				if ((cosXZ > cosThrXZ) && (distProj1 < distProjThrXZ) && (distProj2 < distProjThrXZ))
					return true;
			}
		}
	}
	return false;
}
// ------------------------------------------------------

void PMAlgTrackMaker::mergeCoLinear(pma::trk_candidates& tracks)
{
	double distThr = 0.05;    // max gap as a fraction of the longer track length
	double distThrMin = 0.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fMergeTransverseShift;
	double cosThr = cos(TMath::Pi() * fMergeAngle / 180.0);

	bool r;
	double d, c;
	size_t t = 0, u = 0;
	while (t < tracks.size())
	{
		pma::Track3D* trk1 = tracks[t].Track();

		pma::Track3D* trk2 = 0;
		for (u = t + 1; u < tracks.size(); u++)
		{
			trk2 = tracks[u].Track();

			if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr)) break;

			trk2 = 0;
		}

		if (trk2)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Merge track ("
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
			tracks.erase(tracks.begin() + u);
		}
		else t++;
	}
}
// ------------------------------------------------------

void PMAlgTrackMaker::mergeCoLinear(tpc_track_map& tracks)
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
		pma::trk_candidates& tracks1 = tpc_entry1.second;

		size_t t = 0;
		while (t < tracks1.size())
		{
			bool r, reverse = false;
			double d, c, cmax = 0.0;
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

					pma::trk_candidates& tracks2 = tpc_entry2.second;

					for (size_t u = 0; u < tracks2.size(); u++)
					{
						pma::Track3D* trk2 = tracks2[u].Track();
						dfront2 = trk2->Nodes().front()->GetDistToWall();
						dback2 = trk2->Nodes().back()->GetDistToWall();
						if ((dfront2 < wallDistThr) || (dback2 < wallDistThr))
						{
							if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr))
							{
								if (c > cmax)
								{
									cmax = c;
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
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge track ("
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
				tracks[best_tpc].erase(tracks[best_tpc].begin() + best_idx);
			}
			else t++;
		}
	}

	//for (auto & tpc_entry : tracks) releaseAllNodes(tpc_entry.second);
}
// ------------------------------------------------------

bool PMAlgTrackMaker::areCoLinear(double& cos3d,
		TVector3 f0, TVector3 b0, TVector3 f1, TVector3 b1,
		double distProjThr)
{
	TVector3 s0 = b0 - f0, s1 = b1 - f1;
	cos3d = s0 * s1 / (s0.Mag() * s1.Mag());

	TVector3 proj0 = pma::GetProjectionToSegment(b0, f1, b1);
	double distProj0 = sqrt( pma::Dist2(b0, proj0) );

	TVector3 proj1 = pma::GetProjectionToSegment(b1, f0, b0);
	double distProj1 = sqrt( pma::Dist2(b1, proj1) );

	double d = sqrt( pma::Dist2(b0, f1) );
	double dThr = (1 + 0.02 * d) * distProjThr;

	mf::LogVerbatim("PMAlgTrackMaker")
		<< "   dThr:" << dThr << " d0:" << distProj0 << " d1:" << distProj1 << " c:" << cos3d;

	if ((distProj0 < dThr) && (distProj1 < dThr))
		return true;
	else return false;
}
// ------------------------------------------------------

void PMAlgTrackMaker::matchCoLinearAnyT0(pma::trk_candidates& tracks)
{
	double distProjThr = fStitchTransverseShift;
	double cosThr = cos(TMath::Pi() * fStitchAngle / 180.0);
	double xApaDistDiffThr = 1.0;

	for (size_t u = 0; u < tracks.size(); u++)
	{
		pma::Track3D* trk1 = tracks[u].Track();

		unsigned int tpcFront1 = trk1->FrontTPC(), cryoFront1 = trk1->FrontCryo();
		unsigned int tpcBack1 = trk1->BackTPC(), cryoBack1 = trk1->BackCryo();

		unsigned int firstPlane = 0;
		while ((firstPlane < 3) && !fGeom->TPC(tpcFront1, cryoFront1).HasPlane(firstPlane)) firstPlane++;

		double dxFront1 = trk1->front()->Point3D().X();
		dxFront1 -= fGeom->TPC(tpcFront1, cryoFront1).PlaneLocation(firstPlane)[0];

		firstPlane = 0;
		while ((firstPlane < 3) && !fGeom->TPC(tpcBack1, cryoBack1).HasPlane(firstPlane)) firstPlane++;

		double dxBack1 = trk1->back()->Point3D().X();
		dxBack1 -= fGeom->TPC(tpcBack1, cryoBack1).PlaneLocation(firstPlane)[0];

		//size_t best_idx = 0;
		pma::Track3D* best_trk2 = 0;
		double dx1 = 0.0, dx2 = 0.0, c, cmax = cosThr;
		bool reverse = false, flip1 = false, flip2 = false;
		TVector3 f0, b0, f1, b1;
		for (size_t t = u + 1; t < tracks.size(); t++)
		{
			pma::Track3D* trk2 = tracks[t].Track();
			if (trk2->GetXShift() != 0.0) continue;

			unsigned int tpcFront2 = trk2->FrontTPC(), cryoFront2 = trk2->FrontCryo();
			unsigned int tpcBack2 = trk2->BackTPC(), cryoBack2 = trk2->BackCryo();

			firstPlane = 0;
			while ((firstPlane < 3) && !fGeom->TPC(tpcFront2, cryoFront2).HasPlane(firstPlane)) firstPlane++;

			double dxFront2 = trk2->front()->Point3D().X();
			dxFront2 -= fGeom->TPC(tpcFront2, cryoFront2).PlaneLocation(firstPlane)[0];

			firstPlane = 0;
			while ((firstPlane < 3) && !fGeom->TPC(tpcBack2, cryoBack2).HasPlane(firstPlane)) firstPlane++;

			double dxBack2 = trk2->back()->Point3D().X();
			if (dxBack2 > 0.0) dxBack2 -= fGeom->TPC(tpcBack2, cryoBack2).PlaneLocation(firstPlane)[0];

			mf::LogVerbatim("PMAlgTrackMaker")
				<< "   xf1:" << dxFront1 << " xb1:" << dxBack1
				<< "   xf2:" << dxFront2 << " xb2:" << dxBack2;

			if ((cryoFront1 == cryoFront2) && (dxFront1 * dxFront2 < 0.0) &&
			    (fabs(dxFront1 + dxFront2) < xApaDistDiffThr))
			{
				if (fabs(dxFront1) < fabs(dxFront2)) dxFront2 = -dxFront1;
				else dxFront1 = -dxFront2;

				f0 = trk1->Nodes()[1]->Point3D(); f0.SetX(f0.X() - dxFront1);
				b0 = trk1->Nodes()[0]->Point3D(); b0.SetX(b0.X() - dxFront1);
				f1 = trk2->Nodes()[0]->Point3D(); f1.SetX(f1.X() - dxFront2);
				b1 = trk2->Nodes()[1]->Point3D(); b1.SetX(b1.X() - dxFront2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = false; flip1 = true; flip2 = false;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxFront1; dx2 = dxFront2;
				}
			}
			else if ((cryoFront1 == cryoBack2) && (dxFront1 * dxBack2 < 0.0) &&
			    (fabs(dxFront1 + dxBack2) < xApaDistDiffThr))
			{
				if (fabs(dxFront1) < fabs(dxBack2)) dxBack2 = -dxFront1;
				else dxFront1 = -dxBack2;

				f0 = trk1->Nodes()[1]->Point3D(); f0.SetX(f0.X() - dxFront1);
				b0 = trk1->Nodes()[0]->Point3D(); b0.SetX(b0.X() - dxFront1);
				f1 = trk2->Nodes()[trk2->Nodes().size() - 1]->Point3D(); f1.SetX(f1.X() - dxBack2);
				b1 = trk2->Nodes()[trk2->Nodes().size() - 2]->Point3D(); b1.SetX(b1.X() - dxBack2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = true; flip1 = false; flip2 = false;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxFront1; dx2 = dxBack2;
				}
			}
			else if ((cryoBack1 == cryoFront2) && (dxBack1 * dxFront2 < 0.0) &&
			    (fabs(dxBack1 + dxFront2) < xApaDistDiffThr))
			{
				if (fabs(dxBack1) < fabs(dxFront2)) dxFront2 = -dxBack1;
				else dxBack1 = -dxFront2;

				f0 = trk1->Nodes()[trk1->Nodes().size() - 2]->Point3D(); f0.SetX(f0.X() - dxBack1);
				b0 = trk1->Nodes()[trk1->Nodes().size() - 1]->Point3D(); b0.SetX(b0.X() - dxBack1);
				f1 = trk2->Nodes()[0]->Point3D(); f1.SetX(f1.X() - dxFront2);
				b1 = trk2->Nodes()[1]->Point3D(); b1.SetX(b1.X() - dxFront2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = false; flip1 = false; flip2 = false;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxBack1; dx2 = dxFront2;
				}
			}
			else if ((cryoBack1 == cryoBack2) && (dxBack1 * dxBack2 < 0.0) &&
			    (fabs(dxBack1 + dxBack2) < xApaDistDiffThr))
			{
				if (fabs(dxBack1) < fabs(dxBack2)) dxBack2 = -dxBack1;
				else dxBack1 = -dxBack2;

				f0 = trk1->Nodes()[trk1->Nodes().size() - 2]->Point3D(); f0.SetX(f0.X() - dxBack1);
				b0 = trk1->Nodes()[trk1->Nodes().size() - 1]->Point3D(); b0.SetX(b0.X() - dxBack1);
				f1 = trk2->Nodes()[trk2->Nodes().size() - 1]->Point3D(); f1.SetX(f1.X() - dxBack2);
				b1 = trk2->Nodes()[trk2->Nodes().size() - 2]->Point3D(); b1.SetX(b1.X() - dxBack2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = false; flip1 = false; flip2 = true;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxBack1; dx2 = dxBack2;
				}
			}
		}
		if (best_trk2)
		{
/*			if (reverse)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Match cross-APA tracks: ("
					<< best_idx << ":" << best_trk2->size() << ") and ("
					<< u << ":" << trk1->size() << ")";
			}
			else
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Match cross-APA tracks: ("
					<< u << ":" << trk1->size() << ") and ("
					<< best_idx << ":" << best_trk2->size() << ")";
			}
*/
			if (flip1) trk1->Flip();
			if (flip2) best_trk2->Flip();

			trk1->GetRoot()->ApplyXShiftInTree(-dx1);
			best_trk2->GetRoot()->ApplyXShiftInTree(-dx2);

			if (reverse)
			{
				best_trk2->SetSubsequentTrack(trk1);
				trk1->SetPrecedingTrack(best_trk2);
			}
			else
			{
				trk1->SetSubsequentTrack(best_trk2);
				best_trk2->SetPrecedingTrack(trk1);
			}
		}
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::reassignHits(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	pma::trk_candidates& tracks, size_t trk_idx)
{
	pma::Track3D* trk1 = tracks[trk_idx].Track();
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
		mf::LogVerbatim("PMAlgTrackMaker") << "  Reassign " << n_max << " hits." << std::endl;

		trk1->RemoveHits(hits);
		trk1->CleanupTails();
		trk1->ShiftEndsToHits();

		pma::Track3D* ext = fProjectionMatchingAlg.extendTrack(*best_trk, hits,	false);
		ext->SortHits(); ext->ShiftEndsToHits();
		if (fProjectionMatchingAlg.isContained(*ext))
		{
			tracks[best_u].SetTrack(ext); // and this deletes best_trk stored at best_u
			return true;
		}
		else delete ext;
	}
	return false;
}

bool PMAlgTrackMaker::reassignSingleViewEnds(pma::trk_candidates& tracks)
{
	bool result = false;
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D* trk = tracks[t].Track();
		if (trk->size() < 5) continue;

		trk->DisableSingleViewEnds();

		std::vector< art::Ptr<recob::Hit> > hits;

		size_t idx = 0;
		while ((idx < trk->size() - 1) && !((*trk)[idx]->IsEnabled()))
		{
			hits.push_back((*trk)[idx++]->Hit2DPtr());
		}

		double d2;
		if (idx > 0)
		{
			if ((idx < trk->size() - 1) &&
			    ((*trk)[idx]->View2D() == (*trk)[idx - 1]->View2D()))
			{
				double dprev = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx - 1]->Point3D());
				double dnext = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx + 1]->Point3D());
				if (dprev < dnext)
				{
					hits.push_back((*trk)[idx++]->Hit2DPtr());
				}
			}
			d2 = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx - 1]->Point3D());
		}
		else d2 = 0.0;

		if ((hits.size() > 1) || (d2 > 1.0)) // min. 2 hits or single hit separated from the rest
		{
			result |= reassignHits(hits, tracks, t);
		}

		hits.clear();
		idx = trk->size() - 1;
		while ((idx > 0) && !((*trk)[idx]->IsEnabled()))
		{
			hits.push_back((*trk)[idx--]->Hit2DPtr());
		}

		if (idx < trk->size() - 1)
		{
			if ((idx > 0) &&
			    ((*trk)[idx]->View2D() == (*trk)[idx + 1]->View2D()))
			{
				double dprev = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx + 1]->Point3D());
				double dnext = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx - 1]->Point3D());
				if (dprev < dnext)
				{
					hits.push_back((*trk)[idx--]->Hit2DPtr());
				}
			}
			d2 = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx + 1]->Point3D());
		}
		else d2 = 0.0;

		if ((hits.size() > 1) || (d2 > 1.0))  // min. 2 hits or single hit separated from the rest
		{
			result |= reassignHits(hits, tracks, t);
		}

		trk->SelectHits();
	}
	return result;
}
// ------------------------------------------------------

void PMAlgTrackMaker::guideEndpoints(pma::trk_candidates& tracks)
{
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D& trk = *(tracks[t].Track());
		fProjectionMatchingAlg.guideEndpoints(trk, fHitMap[trk.FrontCryo()][trk.FrontTPC()]);
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::sortHits(const art::Event& evt)
{
	fHitMap.clear(); fCluHits.clear();

	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle, splitCluHandle;
	std::vector< art::Ptr<recob::Hit> > allhitlist;
	if (evt.getByLabel(fHitModuleLabel, splitCluHandle) &&   // clusters that tag em-like hits
	    evt.getByLabel(fCluModuleLabel, allHitListHandle) && // all hits associated to both cluster sets
	    evt.getByLabel(fCluModuleLabel, cluListHandle))      // clusters used to build 3D tracks
	{
		art::fill_ptr_vector(allhitlist, allHitListHandle);

		mf::LogVerbatim("PMAlgTrackMaker") << "Sort all hits for validation...";
		unsigned int cryo, tpc, view;
		for (auto const& h : allhitlist) // all hits used for validation
		{
			cryo = h->WireID().Cryostat;
			tpc = h->WireID().TPC;
			view = h->WireID().Plane;

			fHitMap[cryo][tpc][view].push_back(h);
		}
		mf::LogVerbatim("PMAlgTrackMaker") << "...done.";

		mf::LogVerbatim("PMAlgTrackMaker") << "Filter track-like clusters...";
		fCluHits.reserve(cluListHandle->size());
		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);
		art::FindManyP< recob::Hit > fem(splitCluHandle, evt, fHitModuleLabel);
		for (size_t i = 0; i < cluListHandle->size(); ++i)
		{
			auto v = fbp.at(i);

			fCluHits.emplace_back(std::vector< art::Ptr<recob::Hit> >());

			for (auto const& h : v)
			{
				bool trkLike = true;
				if (fCluModuleLabel != fHitModuleLabel)
				{
					for (size_t j = 0; j < splitCluHandle->size(); ++j)
					{
						auto u = fem.at(j);
						for (auto const& g : u) // is hit clustered in one of em-like?
						{
							if (g.key() == h.key())
							{
								trkLike = false; break;
							}
						}
					}
				}
				if (trkLike) fCluHits.back().push_back(h);
			}
		}
		if (fCluHits.size() != cluListHandle->size())
		{
			mf::LogError("PMAlgTrackMaker") << "Hit-cluster map incorrect, better skip this event.";
			return false;
		}

		mf::LogVerbatim("PMAlgTrackMaker") << "...done, " << fCluHits.size() << " clusters for 3D tracking.";
		return true;
	}
	else return false;
}
// ------------------------------------------------------

bool PMAlgTrackMaker::sortHitsPfp(const art::Event& evt)
{
	fHitMap.clear(); fCluHits.clear(); fPfpClusters.clear();

	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
	std::vector< art::Ptr<recob::Hit> > allhitlist;
	if (evt.getByLabel(fHitModuleLabel, allHitListHandle) && // all hits used to make clusters and PFParticles
	    evt.getByLabel(fCluModuleLabel, cluListHandle) &&    // clusters associated to PFParticles
	    evt.getByLabel(fCluModuleLabel, pfparticleHandle))   // and finally PFParticles
	{
		art::fill_ptr_vector(allhitlist, allHitListHandle);

		mf::LogVerbatim("PMAlgTrackMaker") << "Sort all hits for validation...";
		unsigned int cryo, tpc, view;
		for (auto const& h : allhitlist) // all hits used for validation
		{
			cryo = h->WireID().Cryostat;
			tpc = h->WireID().TPC;
			view = h->WireID().Plane;

			fHitMap[cryo][tpc][view].push_back(h);
		}
		mf::LogVerbatim("PMAlgTrackMaker") << "...done.";

		mf::LogVerbatim("PMAlgTrackMaker") << "Sort hits by clusters assigned to PFParticles...";
		fCluHits.reserve(cluListHandle->size());
		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);
		art::FindManyP< recob::Cluster > fpf(pfparticleHandle, evt, fCluModuleLabel);
		for (size_t i = 0; i < cluListHandle->size(); ++i)
		{
			fCluHits.emplace_back(std::vector< art::Ptr<recob::Hit> >());
		}
		for (size_t i = 0; i < pfparticleHandle->size(); ++i)
		{
			auto cv = fpf.at(i);
			for (const auto & c : cv)
			{
				fPfpClusters[i].push_back(c);

				if (fCluHits[c.key()].empty())
				{
					auto hv = fbp.at(c.key());
					fCluHits[c.key()].reserve(hv.size());
					for (auto const & h : hv) fCluHits[c.key()].push_back(h);
				}
			}
		}

		mf::LogVerbatim("PMAlgTrackMaker") << "...done, "
			<< fCluHits.size() << " clusters from "
			<< fPfpClusters.size() << " pfparticles for 3D tracking.";
		return true;
	}
	else return false;
}
// ------------------------------------------------------

void PMAlgTrackMaker::produce(art::Event& evt)
{
	reset(evt); // set default values, clear containers at the beginning of each event

	pma::trk_candidates result;

	std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
	std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);
	std::unique_ptr< std::vector< recob::Vertex > > vtxs(new std::vector< recob::Vertex >);
	std::unique_ptr< std::vector< anab::T0 > > t0s(new std::vector< anab::T0 >);

	std::unique_ptr< art::Assns< recob::Track, recob::Hit > > trk2hit_oldway(new art::Assns< recob::Track, recob::Hit >); // ****** REMEMBER to remove when FindMany improved ******
	std::unique_ptr< art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta > > trk2hit(new art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta >);

	std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
	std::unique_ptr< art::Assns< recob::Track, anab::T0 > > trk2t0(new art::Assns< recob::Track, anab::T0 >);

	std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);
	std::unique_ptr< art::Assns< recob::Vertex, recob::Track > > vtx2trk(new art::Assns< recob::Vertex, recob::Track >);


	std::unique_ptr< std::vector< recob::PFParticle > > pfps(new std::vector< recob::PFParticle >);

    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> > pfp2clu( new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > pfp2vtx( new art::Assns<recob::PFParticle, recob::Vertex> );
	std::unique_ptr< art::Assns< recob::PFParticle, recob::Track > > pfp2trk(new art::Assns< recob::PFParticle, recob::Track >);

	bool sortHitsClustersOK = false;
	switch (fCluMatchingAlg)
	{
		default: // try to match from all clusters in the event
		case 1: sortHitsClustersOK = sortHits(evt); break;

		case 2: // take clusters-hit assns from PFP, keep all hits for validation
		case 3: sortHitsClustersOK = sortHitsPfp(evt); break;
	}

	if (sortHitsClustersOK)
	{
		int retCode = 0;
		switch (fCluMatchingAlg)
		{
			default:
			case 1: retCode = fromMaxCluster(evt, result); break; // try to match from all clusters in the event
			case 2: retCode = fromPfpClusterSubset(evt, result); break; // each trk matched only from clusters assigned to PFP
			case 3: retCode = fromPfpDirect(evt, result); break; // no pattern recognition, just take clusters assigned to PFP
		}
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

		if (!result.empty()) // ok, there is something to save
		{
			size_t spStart = 0, spEnd = 0;
			double sp_pos[3], sp_err[6];
			for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

			double dQdxFlipThr = 0.0;
			if (fFlipToBeam) dQdxFlipThr = 0.4;

            // use the following to create PFParticle <--> Track associations;
			// note: these are assns to existing PFParticles, that are used for CluMatchingAlg = 2 or 3.
            std::map< size_t, std::vector< art::Ptr<recob::Track> > > pfPartToTrackVecMap;

			tracks->reserve(result.size());
			for (fTrkIndex = 0; fTrkIndex < (int)result.size(); ++fTrkIndex)
			{
				pma::Track3D* trk = result[fTrkIndex].Track();
				if (!(trk->HasTwoViews() && (trk->Nodes().size() > 1)))
				{
					mf::LogWarning("PMAlgTrackMaker") << "Skip degenerated track, code needs to be corrected.";
					continue;
				}

				if (trk->CanFlip())
				{
					if (fFlipToBeam)    // flip the track to the beam direction
					{
						double z0 = trk->front()->Point3D().Z();
						double z1 = trk->back()->Point3D().Z();
						if (z0 > z1) trk->Flip();
					}
					if (fFlipDownward)  // flip the track to point downward
					{
						double y0 = trk->front()->Point3D().Y();
						double y1 = trk->back()->Point3D().Y();
						if (y0 < y1) trk->Flip();
					}
					if (fAutoFlip_dQdx) // flip the track by dQ/dx
						fProjectionMatchingAlg.autoFlip(*trk, pma::Track3D::kForward, dQdxFlipThr);
						/* test code: fProjectionMatchingAlg.autoFlip(*trk, pma::Track3D::kBackward, dQdxFlipThr); */
				}

				trk->SelectHits();  // just in case, set all to enabled
				unsigned int itpc = trk->FrontTPC(), icryo = trk->FrontCryo();
				if (fGeom->TPC(itpc, icryo).HasPlane(geo::kU)) trk->CompleteMissingWires(geo::kU);
				if (fGeom->TPC(itpc, icryo).HasPlane(geo::kV)) trk->CompleteMissingWires(geo::kV);
				if (fGeom->TPC(itpc, icryo).HasPlane(geo::kZ)) trk->CompleteMissingWires(geo::kZ);

				tracks->emplace_back(convertFrom(*trk));

				double xShift = trk->GetXShift();
				if (xShift > 0.0)
				{
					double tisk2time = 1.0; // what is the coefficient, offset?
					double t0time = tisk2time * xShift / fDetProp->GetXTicksCoefficient(trk->FrontTPC(), trk->FrontCryo());

					// TriggBits=3 means from 3d reco (0,1,2 mean something else)
					t0s->push_back(anab::T0(t0time, 0, 3, tracks->back().ID()));
					util::CreateAssn(*this, evt, *tracks, *t0s, *trk2t0, t0s->size() - 1, t0s->size());
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

                // if there is a PFParticle collection then recover PFParticle and add info to map
                if (!fMakePFPs && (result[fTrkIndex].Key() > -1))
                {
                    size_t trackIdx = tracks->size() - 1;
                    art::ProductID trackId = getProductID< std::vector<recob::Track> >(evt);
                    art::Ptr<recob::Track> trackPtr(trackId, trackIdx, evt.productGetter(trackId));
                    pfPartToTrackVecMap[result[fTrkIndex].Key()].push_back(trackPtr);
                }
			}

			auto pfpid = getProductID< std::vector<recob::PFParticle> >(evt);
			auto vid = getProductID< std::vector<recob::Vertex> >(evt);
			auto tid = getProductID< std::vector<recob::Track> >(evt);
			auto const* trkGetter = evt.productGetter(tid);

			auto vsel = fPMAlgVertexing.getVertices(result); // vertex positions with vector of connected tracks idxs
			std::map< size_t, art::Ptr<recob::Vertex> > frontVtxs; // front vertex ptr for each track index

			if (fRunVertexing) // save vertices and vtx-trk assns
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
					if (!v.second.empty())
					{
						frontVtxs[v.second.front()] = vptr; // keep ptr of the front vtx
						for (size_t tidx : v.second)
						{
							art::Ptr<recob::Track> tptr(tid, tidx, trkGetter);
							vtx2trk->addSingle(vptr, tptr);
						}
					}
				}
				mf::LogVerbatim("Summary") << vtxs->size() << " vertices ready";
			}

			if (fMakePFPs)
			{
				// first particle, to be replaced with nu reco when possible
				pfps->emplace_back(recob::PFParticle(0, 0, 0, std::vector< size_t >()));

				pma::setParentDaughterConnections(result);
				for (size_t t = 0; t < result.size(); ++t)
				{
					size_t parentIdx = 0;
					if (result[t].Parent() >= 0) parentIdx = (size_t)result[t].Parent() + 1;

					std::vector< size_t > daughterIdxs;
					for (size_t idx : result[t].Daughters()) daughterIdxs.push_back(idx + 1);

					size_t pfpidx = pfps->size();
					pfps->emplace_back(recob::PFParticle(0, pfpidx, parentIdx, daughterIdxs));

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
			else
			{
	            // if we have used existing PFParticles then do the associations here
    	        if (!pfPartToTrackVecMap.empty())
    	        {
					art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
					evt.getByLabel(fCluModuleLabel, pfParticleHandle);
    	            for (const auto & pfParticleItr : pfPartToTrackVecMap)
    	            {
    	                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfParticleItr.first);
    	                mf::LogVerbatim("PMAlgTrackMaker") << "PFParticle key: " << pfParticle.key()
							<< ", self: " << pfParticle->Self() << ", #tracks: " << pfParticleItr.second.size();

    	                if (!pfParticle.isNull()) util::CreateAssn(*this, evt, pfParticle, pfParticleItr.second, *pfp2trk);
						else mf::LogError("PMAlgTrackMaker") << "Error in PFParticle lookup, pfparticle index: "
							<< pfParticleItr.first << ", key: " << pfParticle.key();
    	            }
    	        }
			}

			// data prods done, delete all pma::Track3D's
			for (auto t : result) t.DeleteTrack();
		}
	}
	else mf::LogError("PMAlgTrackMaker") << "Hits not found in the event.";

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(vtxs));
	evt.put(std::move(t0s));

	evt.put(std::move(trk2hit_oldway)); // ****** REMEMBER to remove when FindMany improved ******
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(trk2t0));

	evt.put(std::move(sp2hit));
	evt.put(std::move(vtx2trk));

	evt.put(std::move(pfps));
	evt.put(std::move(pfp2clu));
	evt.put(std::move(pfp2vtx));
	evt.put(std::move(pfp2trk));
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromMaxCluster(const art::Event& evt, pma::trk_candidates& result)
{
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		initial_clusters.clear();
		tried_clusters.clear();
		used_clusters.clear();

        std::vector< art::Ptr<recob::Cluster> > clusterVec;
        art::fill_ptr_vector(clusterVec, cluListHandle); // use all clusters

		tpc_track_map tracks; // track parts in tpc's

		// find reasonably large parts
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(tracks[tpc_iter->TPC], clusterVec, fMinSeedSize1stPass, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// loop again to find small things
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(tracks[tpc_iter->TPC], clusterVec, fMinSeedSize2ndPass, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// try correcting track ends:
		//   - 3D ref.points for clean endpoints of wire-plae parallel tracks
		//   - single-view sections spuriously merged on 2D clusters level
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			guideEndpoints(tracks[tpc_iter->TPC]);
			reassignSingleViewEnds(tracks[tpc_iter->TPC]);
		}

		if (fMergeWithinTPC)
		{
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge co-linear tracks within TPC " << tpc_iter->TPC << ".";
				mergeCoLinear(tracks[tpc_iter->TPC]);
			}
		}

		if (fStitchBetweenTPCs)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Stitch co-linear tracks between TPCs.";
			mergeCoLinear(tracks);
		}

		for (auto const & tpc_entry : tracks) // put tracks in the single collection
			for (auto & trk : tpc_entry.second)
				if (trk.Track()->HasTwoViews() && (trk.Track()->Nodes().size() > 1))
		{
			fProjectionMatchingAlg.setTrackTag(*(trk.Track())); // tag EM-like tracks
			result.push_back(trk);
		}

		if (fRunVertexing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Vertex finding / track-vertex reoptimization.";
			fPMAlgVertexing.run(result);
		}

		if (fMatchT0inAPACrossing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Find co-linear APA-crossing tracks with any T0.";
			matchCoLinearAnyT0(result);
		}

		listUsedClusters(clusterVec);
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "no clusters";
		return -1;
	}

	return result.size();
}

int PMAlgTrackMaker::matchCluster(const pma::TrkCandidate& trk,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo)
{
	double f, fmax = 0.0;
	unsigned int n, max = 0;
	int idx = -1;
	for (size_t i = 0; i < clusters.size(); ++i)
	{
		unsigned int view = clusters[i]->View();
		unsigned int nhits = fCluHits.at(clusters[i].key()).size();

		if (has(used_clusters, i) ||                             // don't try already used clusters
			has(trk.Clusters(), i) ||                            // don't try clusters from this candidate
		    (view == testView) ||                                // don't use clusters from validation view
		    ((preferedView != geo::kUnknown)&&(view != preferedView)) || // only prefered view if specified
		    (nhits < minSize))                                   // skip small clusters
		    continue;

		n = fProjectionMatchingAlg.testHits(*(trk.Track()), fCluHits.at(clusters[i].key()));
		f = n / (double)nhits;
		if ((f > fraction) && (n > max))
		{
			max = n; fmax = f; idx = i;
		}
	}

	if (idx >= 0) mf::LogVerbatim("PMAlgTrackMaker") << "max matching hits: " << max << " (" << fmax << ")";
	else mf::LogVerbatim("PMAlgTrackMaker") << "no clusters to extend the track";

	return idx;
}


void PMAlgTrackMaker::fromMaxCluster_tpc(
	pma::trk_candidates& result,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minBuildSize, unsigned int tpc, unsigned int cryo,
	int pfParticleIdx)
{
	initial_clusters.clear();

	size_t minSizeCompl = minBuildSize / 8;  // smaller minimum required in complementary views
	if (minSizeCompl < 2) minSizeCompl = 2;  // but at least two hits!

	int max_first_idx = 0;
	while (max_first_idx >= 0) // loop over clusters, any view, starting from the largest
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "Find max cluster...";
		max_first_idx = maxCluster(clusters, minBuildSize, geo::kUnknown, tpc, cryo); // any view
		if (max_first_idx >= 0)
		{
			geo::View_t first_view = clusters[max_first_idx]->View();
			geo::View_t sec_view_a, sec_view_b;
			switch (first_view)
			{
				case geo::kU: sec_view_a = geo::kZ; sec_view_b = geo::kV; break;
				case geo::kV: sec_view_a = geo::kZ; sec_view_b = geo::kU; break;
				case geo::kZ: sec_view_a = geo::kV; sec_view_b = geo::kU; break;
				default:
					mf::LogError("PMAlgTrackMaker") << "Not a 2D view.";
					return;
			}

			tried_clusters[geo::kU].clear();
			tried_clusters[geo::kV].clear();
			tried_clusters[geo::kZ].clear();

			tried_clusters[first_view].push_back(max_first_idx);
			initial_clusters.push_back(max_first_idx);

			const auto & v_first = fCluHits.at(clusters[max_first_idx].key());
			unsigned int nFirstHits = v_first.size();
			mf::LogVerbatim("PMAlgTrackMaker") << "--- start new candidate ---";
			mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << first_view << " ***  size: " << nFirstHits;

			float tmax = fDetProp->ConvertTicksToX(v_first.front()->PeakTime(), first_view, tpc, cryo);
			float t, tmin = tmax;
			for (size_t j = 1; j < v_first.size(); ++j)
			{
				t = fDetProp->ConvertTicksToX(v_first[j]->PeakTime(), first_view, tpc, cryo);
				if (t > tmax) { tmax = t; }
				if (t < tmin) { tmin = t; }
			}

			fCandidates.clear(); // temporary set of possible solutions of the selected cluster and clusters in complementary views

			size_t imatch = 0;
			bool try_build = true;
			while (try_build) // loop over complementary views
			{
				pma::TrkCandidate candidate;
				candidate.Clusters().push_back(max_first_idx);
				candidate.SetKey(pfParticleIdx);

				int idx, max_sec_a_idx, max_sec_b_idx;
				max_sec_a_idx = maxCluster(max_first_idx, clusters, tmin, tmax, minSizeCompl, sec_view_a, tpc, cryo);
				max_sec_b_idx = maxCluster(max_first_idx, clusters, tmin, tmax, minSizeCompl, sec_view_b, tpc, cryo);

				unsigned int nSecHitsA = 0, nSecHitsB = 0;
				if (max_sec_a_idx >= 0) nSecHitsA = fCluHits.at(clusters[max_sec_a_idx].key()).size();
				if (max_sec_b_idx >= 0) nSecHitsB = fCluHits.at(clusters[max_sec_b_idx].key()).size();

				unsigned int testView = geo::kUnknown;
				if ((nSecHitsA > nSecHitsB) && (nSecHitsA >= minSizeCompl))
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "--> " << imatch++ << " match with:";
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_a << " ***  size: " << nSecHitsA;
					tried_clusters[sec_view_a].push_back(max_sec_a_idx);
					idx = max_sec_a_idx; testView = sec_view_b;
				}
				else if (nSecHitsB >= minSizeCompl)
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "--> " << imatch++ << " match with:";
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_b << " ***  size: " << nSecHitsB;
					tried_clusters[sec_view_b].push_back(max_sec_b_idx);
					idx = max_sec_b_idx; testView = sec_view_a;
				}
				else try_build = false;

				if (try_build)
				{
					if (!fGeom->TPC(tpc, cryo).HasPlane(testView)) testView = geo::kUnknown;

					double m0 = 0.0, v0 = 0.0;
					double mseThr = 0.15, validThr = 0.7; // cuts for a good track candidate

					candidate.Clusters().push_back(idx);
					candidate.SetTrack(fProjectionMatchingAlg.buildTrack(v_first, fCluHits.at(clusters[idx].key())));

					if (candidate.IsValid() && // no track if hits from 2 views do not alternate
					    fProjectionMatchingAlg.isContained(*(candidate.Track()))) // sticks out of TPC's?
					{
						m0 = candidate.Track()->GetMse();
						if (m0 < mseThr) // check validation only if MSE is passing - thanks for Tracy for noticing this
							v0 = validate(*(candidate.Track()), testView);
					}
					if (candidate.Track() && (m0 < mseThr) && (v0 > validThr)) // good candidate, try to extend it
					{
						mf::LogVerbatim("PMAlgTrackMaker")
							<< "  good track candidate, MSE = " << m0 << ", v = " << v0;

						candidate.SetMse(m0);
						candidate.SetValidation(v0);
						candidate.SetGood(true);

						size_t minSize = 5;      // min size for clusters matching
						double fraction = 0.5;   // min fraction of close hits

						idx = 0;
						while (idx >= 0) // try to collect matching clusters, use **any** plane except validation
						{
							idx = matchCluster(candidate, clusters, minSize, fraction, geo::kUnknown, testView, tpc, cryo);
							if (idx >= 0)
							{
								// try building extended copy:
								//                src,        hits,    valid.plane, add nodes
								if (extendTrack(candidate, fCluHits.at(clusters[idx].key()),  testView,    true))
								{
									candidate.Clusters().push_back(idx);
								}
								else idx = -1;
							}
						}

						mf::LogVerbatim("PMAlgTrackMaker") << "merge clusters from the validation plane";
						fraction = 0.7; // only well matching the existing track

						idx = 0;
						bool extended = false;
						while ((idx >= 0) && (testView != geo::kUnknown))
						{	//                     match clusters from the plane used previously for the validation
							idx = matchCluster(candidate, clusters, minSize, fraction, testView, geo::kUnknown, tpc, cryo);
							if (idx >= 0)
							{
								// validation not checked here, no new nodes:
								if (extendTrack(candidate, fCluHits.at(clusters[idx].key()), geo::kUnknown, false))
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
						mf::LogVerbatim("PMAlgTrackMaker") << "track REJECTED, MSE = " << m0 << "; v = " << v0;
						candidate.SetGood(false); // save also bad matches to avoid trying again the same pair of clusters
					}
					fCandidates.emplace_back(candidate);
				}
				else
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "no matching clusters";
				}
			} // end loop over complementary views

			if (fCandidates.size()) // save best candidate, release other tracks and clusters
			{
				int best_trk = -1;
				double f, max_f = 0., min_mse = 10., max_v = 0.;
				for (size_t t = 0; t < fCandidates.size(); t++)
					if (fCandidates[t].IsGood() &&
					    (fCandidates[t].Track()->Nodes().size() > 1) &&
					    fCandidates[t].Track()->HasTwoViews())
				{
					f = fProjectionMatchingAlg.twoViewFraction(*(fCandidates[t].Track()));

					if ((f > max_f) || ((f == max_f) &&
						((fCandidates[t].Validation() > max_v) || (fCandidates[t].Mse() < min_mse))))
					{
						max_f = f;
						min_mse = fCandidates[t].Mse();
						max_v = fCandidates[t].Validation();
						best_trk = t;
					}
				}

				if ((best_trk > -1) && fCandidates[best_trk].IsGood() && (max_f > fMinTwoViewFraction))
				{
					fCandidates[best_trk].Track()->ShiftEndsToHits();

					result.push_back(fCandidates[best_trk]);

					for (auto c : fCandidates[best_trk].Clusters())
						used_clusters.push_back(c);
                }

				for (size_t t = 0; t < fCandidates.size(); t++)
				{
					if (int(t) != best_trk) fCandidates[t].DeleteTrack();
				}
				fCandidates.clear();
			}
		}
		else
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "small clusters only";
		}
	} // end loop over clusters, any view, from the largest
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;

	for (size_t i = 0; i < clusters.size(); ++i)
	{
		const auto & v = fCluHits.at(clusters[i].key());

		if (!v.size() ||
		    has(used_clusters, i) ||
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		   ((view != geo::kUnknown) && (clusters[i]->View() != view)))
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

int PMAlgTrackMaker::maxCluster(size_t first_idx,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	float tmin, float tmax, size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;
	double fraction = 0.0;
	float x;

	for (size_t i = 0; i < clusters.size(); ++i)
	{
		if (has(used_clusters, i) ||
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		    (fCluHits.at(clusters[i].key()).size() <  min_clu_size) ||
		    (clusters[i]->View() != view)) continue;

		bool pair_checked = false;
		for (auto const & c : fCandidates)
			if (has(c.Clusters(), first_idx) && has(c.Clusters(), i))
			{
				pair_checked = true; break;
			}
		if (pair_checked) continue;
		    
		const auto & v = fCluHits.at(clusters[i].key());

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = 0;
			for (size_t j = 0; j < v.size(); ++j)
			{
				x = fDetProp->ConvertTicksToX(v[j]->PeakTime(), view, tpc, cryo);
				if ((x >= tmin) && (x <= tmax)) s++;
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
// ------------------------------------------------------

int PMAlgTrackMaker::fromPfpClusterSubset(const art::Event& evt, pma::trk_candidates& result)
{
	// Code from Tracy merged with recent additions to PMA. Still to be changed in order to
	// skip not reasonalbe parts in this configuration.
    if (!fPfpClusters.empty() && !fCluHits.empty())
    {
		tpc_track_map tracks; // track parts in tpc's

		// Armed with all of this information we can begin looping through the PFParticles
		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			const auto & clusterVec = pfpCluEntry.second;

			initial_clusters.clear();
			tried_clusters.clear();
			used_clusters.clear();

			size_t minBuildSize = 2;
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				fromMaxCluster_tpc(tracks[tpc_iter->TPC], clusterVec, minBuildSize, tpc_iter->TPC, tpc_iter->Cryostat, pfPartIdx);
			}
   
			// used for development
			listUsedClusters(clusterVec);
		}

		// try correcting track ends:
		//   - single-view sections spuriously merged on 2D clusters level
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			reassignSingleViewEnds(tracks[tpc_iter->TPC]);
		}

		// try correcting track ends:
		//   - 3D ref.points for clean endpoints of wire-plae parallel tracks
		//   - single-view sections spuriously merged on 2D clusters level
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			guideEndpoints(tracks[tpc_iter->TPC]);
			reassignSingleViewEnds(tracks[tpc_iter->TPC]);
		}

		// merge co-linear parts inside each tpc
		if (fMergeWithinTPC)
		{
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge co-linear tracks within TPC " << tpc_iter->TPC << ".";
				mergeCoLinear(tracks[tpc_iter->TPC]);
			}
		}

		// merge co-linear parts between tpc's
		if (fStitchBetweenTPCs)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Stitch co-linear tracks between TPCs.";
			mergeCoLinear(tracks);
		}

		for (auto const & tpc_entry : tracks)
			for (auto & trk : tpc_entry.second)
				if (trk.Track()->HasTwoViews() && (trk.Track()->Nodes().size() > 1))
		{
			fProjectionMatchingAlg.setTrackTag(*(trk.Track()));
			result.push_back(trk);
		}

		if (fRunVertexing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Vertex finding / track-vertex reoptimization.";
			fPMAlgVertexing.run(result);
		}

		if (fMatchT0inAPACrossing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Find co-linear APA-crossing tracks with any T0.";
			matchCoLinearAnyT0(result);
		}
    }
    else
    {
        mf::LogWarning("PMAlgTrackMaker") << "no clusters, no pfparticles";
        return -1;
    }
    
    return result.size();
}
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromPfpDirect(const art::Event& evt, pma::trk_candidates& result)
{
    if (!fPfpClusters.empty() && !fCluHits.empty())
    {
		for (const auto & pfpCluEntry : fPfpClusters)
		{
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

			candidate.SetTrack(fProjectionMatchingAlg.buildTrack(allHits));

			if (candidate.IsValid() &&
			    candidate.Track()->HasTwoViews() &&
			    (candidate.Track()->Nodes().size() > 1))
			{
	   			result.push_back(candidate);
			}
			else
			{
				candidate.DeleteTrack();
			}
		}

		guideEndpoints(result); // add 3D ref.points for clean endpoints of wire-plae parallel tracks

		if (fRunVertexing) fPMAlgVertexing.run(result);
    }
    else
    {
        mf::LogWarning("PMAlgTrackMaker") << "no clusters, no pfparticles";
        return -1;
    }
    
    return result.size();
}
// ------------------------------------------------------
// ------------------------------------------------------

void PMAlgTrackMaker::listUsedClusters(const std::vector< art::Ptr<recob::Cluster> >& clusters) const
{
	mf::LogVerbatim("PMAlgTrackMaker") << std::endl << "----------- matched clusters: -----------";
	for (size_t i = 0; i < clusters.size(); ++i)
		if (has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    tpc: " << clusters[i]->Plane().TPC
				<< ";\tview: " << clusters[i]->View()
				<< ";\tsize: " << fCluHits.at(clusters[i].key()).size();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "--------- not matched clusters: ---------";
	for (size_t i = 0; i < clusters.size(); ++i)
		if (!has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    tpc: " << clusters[i]->Plane().TPC
				<< ";\tview: " << clusters[i]->View()
				<< ";\tsize: " << fCluHits.at(clusters[i].key()).size();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "-----------------------------------------";
}
// ------------------------------------------------------
// ------------------------------------------------------

DEFINE_ART_MODULE(PMAlgTrackMaker)

} // namespace trkf

