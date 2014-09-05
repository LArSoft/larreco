////////////////////////////////////////////////////////////////////////
//
// Class:       Track3DKalmanHit
// Module Type: producer
// File:        Track3DKalmanHit_module.cc
//
// This module produces RecoBase/Track objects using KalmanFilterAlgorithm.
//
// Configuration parameters:
//
// Hist               - Histogram flag (generate histograms if true).
// UseClusterHits     - Use clustered hits if true.
// UsePFParticleHits  - Use PFParticle hits if true.
// HitModuleLabel     - Module label for unclustered Hits.
// ClusterModuleLabel - Module label for Clusters.
// PFParticleModuleLabel - Module label for PFParticles.
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// DoDedx             - Global dE/dx enable flag.
// MinSeedHits        - Minimum number of hits per track seed.
// MinSeedChopHits:   - Potentially chop seeds that exceed this length.
// MaxChopHits        - Maximum number of hits to chop from each end of seed.
// MaxSeedChiDF       - Maximum seed track chisquare/dof.
// MinSeedSlope       - Minimum seed slope (dx/dz).
// InitialMomentum    - Initial momentum guess.
// KalmanFilterAlg    - Parameter set for KalmanFilterAlg.
// SeedFinderAlg      - Parameter set for seed finder algorithm object.
// SpacePointAlg      - Parmaeter set for space points.
//
// Usage Notes.
//
// 1.  It is an error if UseClusterHits and UsePFParticleHits are both
//     true (throw exception in that case).  However, both can be false,
//     in which case all hits are used as input.
//
// 2.  If clustered hits are used as input, then tracks can span clusters
//     (cluster structure of hits is ignored).
//
// 3.  If PFParticle hits are used as input, then tracks do not span
//     PFParticles.  However, it is possible for one the hits from
//     one PFParticle to give rise to multiple Tracks.
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Seed.h"
#include "RecoAlg/KalmanFilterAlg.h"
#include "RecoAlg/SeedFinderAlgorithm.h"
#include "RecoObjects/KHitContainerWireX.h"
#include "RecoObjects/SurfXYZPlane.h"
#include "RecoObjects/PropXYZPlane.h"
#include "RecoObjects/KHit.h"
#include "Utilities/AssociationUtil.h"

#include "TH1F.h"

// Local functions.

namespace {

  //----------------------------------------------------------------------------
  // Fill a collection of hits used by space points.
  //
  // Arguments:
  //
  // spts - Space point collection.
  // hits - Hit collection to fill.
  //
  void SpacePointsToHits(const std::vector<recob::SpacePoint>& spts,
			 art::PtrVector<recob::Hit>& hits, 
			 trkf::SpacePointAlg const& spalg)
  {
    hits.clear();
    std::set<const recob::Hit*> used_hits;

    // Loop over space points.

    for(std::vector<recob::SpacePoint>::const_iterator ispt = spts.begin();
	ispt != spts.end(); ++ispt) {
      const recob::SpacePoint& spt = *ispt;

      // Loop over hits in space point.

      const art::PtrVector<recob::Hit> spthits = spalg.getAssociatedHits(spt);
      for(art::PtrVector<recob::Hit>::const_iterator ihit = spthits.begin();
	  ihit != spthits.end(); ++ihit) {
	const art::Ptr<recob::Hit> phit = *ihit;
	const recob::Hit& hit = *phit;

	// Check if hit should be added to collection.

	if(used_hits.count(&hit) == 0) {
	  used_hits.insert(&hit);
	  hits.push_back(phit);
	}
      }
    }
  }

  //----------------------------------------------------------------------------
  // Filter a collection of hits (set difference).
  //
  // Arguments:
  //
  // hits      - Hit collection from which hits should be removed.
  // used_hits - Hits to remove.
  //
  void FilterHits(art::PtrVector<recob::Hit>& hits,
		  art::PtrVector<recob::Hit>& used_hits)
  {
    if(used_hits.size() > 0) {

      // Make sure both hit collections are sorted.

      std::stable_sort(hits.begin(), hits.end());
      std::stable_sort(used_hits.begin(), used_hits.end());

      // Do set difference operation.

      art::PtrVector<recob::Hit>::iterator it =
	std::set_difference(hits.begin(), hits.end(),
			    used_hits.begin(), used_hits.end(),
			    hits.begin());

      // Truncate hit collection.

      hits.erase(it, hits.end());
    }
  }
}

namespace trkf {

  class Propagator;

  class Track3DKalmanHit : public art::EDProducer {
  public:

    // Copnstructors, destructor.

    explicit Track3DKalmanHit(fhicl::ParameterSet const & pset);
    virtual ~Track3DKalmanHit();

    // Overrides.

    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

  private:

    // Fcl parameters.

    bool fHist;                         ///< Make histograms.
    bool fUseClusterHits;               ///< Use clustered hits as input.
    bool fUsePFParticleHits;            ///< Use PFParticle hits as input.
    std::string fHitModuleLabel;        ///< Unclustered Hits.
    std::string fClusterModuleLabel;    ///< Clustered Hits.
    std::string fPFParticleModuleLabel; ///< PFParticle label.
    double fMaxTcut;                    ///< Maximum delta ray energy in MeV for restricted dE/dx.
    bool fDoDedx;                       ///< Global dE/dx enable flag.
    int fMinSeedHits;                   ///< Minimum number of hits per track seed.
    int fMinSeedChopHits;               ///< Potentially chop seeds that exceed this length.
    int fMaxChopHits;                   ///< Maximum number of hits to chop from each end of seed.
    double fMaxSeedChiDF;               ///< Maximum seed track chisquare/dof.
    double fMinSeedSlope;               ///< Minimum seed slope (dx/dz).
    double fInitialMomentum;            ///< Initial (or constant) momentum.

    // Algorithm objects.

    KalmanFilterAlg fKFAlg;             ///< Kalman filter algorithm.
    SeedFinderAlgorithm fSeedFinderAlg; ///< Seed finder.
    SpacePointAlg fSpacePointAlg;       ///< Space point algorithm.

    /// Propagator.
    const Propagator* fProp;

    // Histograms.

    TH1F* fHIncChisq;   ///< Incremental chisquare.
    TH1F* fHPull;       ///< Hit pull.

    // Statistics.

    int fNumEvent;    ///< Number of events seen.
    int fNumTrack;    ///< Number of tracks produced.

  };

  DEFINE_ART_MODULE(Track3DKalmanHit)

} // namespace trkf

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// p - Fcl parameters.
///
trkf::Track3DKalmanHit::Track3DKalmanHit(fhicl::ParameterSet const & pset) :
  fHist(false),
  fUseClusterHits(false),
  fUsePFParticleHits(false),
  fMaxTcut(0.),
  fDoDedx(false),
  fMinSeedHits(0),
  fMinSeedChopHits(0),
  fMaxChopHits(0),
  fMaxSeedChiDF(0.),
  fMinSeedSlope(0.),
  fInitialMomentum(0.),
  fKFAlg(pset.get<fhicl::ParameterSet>("KalmanFilterAlg")),
  fSeedFinderAlg(pset.get<fhicl::ParameterSet>("SeedFinderAlg")),
  fSpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
  fProp(0),
  fHIncChisq(0),
  fHPull(0),
  fNumEvent(0),
  fNumTrack(0)
{
  reconfigure(pset);
  produces<std::vector<recob::Track> >();
  produces<std::vector<recob::SpacePoint>            >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  produces<art::Assns<recob::Track, recob::SpacePoint> >();
  produces<art::Assns<recob::SpacePoint, recob::Hit> >();
  produces<art::Assns<recob::PFParticle, recob::Track> >();

  // Report.

  mf::LogInfo("Track3DKalmanHit") 
    << "Track3DKalmanHit configured with the following parameters:\n"
    << "  UseClusterHits = " << fUseClusterHits << "\n"
    << "  HitModuleLabel = " << fHitModuleLabel << "\n"
    << "  ClusterModuleLabel = " << fClusterModuleLabel;
}

//----------------------------------------------------------------------------
/// Destructor.
trkf::Track3DKalmanHit::~Track3DKalmanHit()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// p - Fcl parameter set.
///
void trkf::Track3DKalmanHit::reconfigure(fhicl::ParameterSet const & pset)
{
  fHist = pset.get<bool>("Hist");
  fKFAlg.reconfigure(pset.get<fhicl::ParameterSet>("KalmanFilterAlg"));
  fSeedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("SeedFinderAlg"));
  fSpacePointAlg.reconfigure(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
  fUseClusterHits = pset.get<bool>("UseClusterHits");
  fUsePFParticleHits = pset.get<bool>("UsePFParticleHits");
  fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
  fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
  fMaxTcut = pset.get<double>("MaxTcut");
  fDoDedx = pset.get<bool>("DoDedx");
  fMinSeedHits = pset.get<int>("MinSeedHits");
  fMinSeedChopHits = pset.get<int>("MinSeedChopHits");
  fMaxChopHits = pset.get<int>("MaxChopHits");
  fMaxSeedChiDF = pset.get<double>("MaxSeedChiDF");
  fMinSeedSlope = pset.get<double>("MinSeedSlope");
  fInitialMomentum = pset.get<double>("InitialMomentum");
  if(fProp != 0)
    delete fProp;
  fProp = new PropXYZPlane(fMaxTcut, fDoDedx);
  if(fUseClusterHits && fUsePFParticleHits) {
    throw cet::exception("Track3DKalmanHit")
      << "Using input from both clustered and PFParticle hits.\n";
  }
}

//----------------------------------------------------------------------------
/// Begin job method.
void trkf::Track3DKalmanHit::beginJob()
{
  if(fHist) {

    // Book histograms.

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("hitkalman", "Track3DKalmanHit histograms");

    fHIncChisq = dir.make<TH1F>("IncChisq", "Incremental Chisquare", 100, 0., 20.);
    fHPull = dir.make<TH1F>("Pull", "Hit Pull", 100, -10., 10.);
  }
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// e - Art event.
///
/// This method extracts Hits from the event and produces and adds
/// Track objects.
///
void trkf::Track3DKalmanHit::produce(art::Event & evt)
{
  ++fNumEvent;

  // Make a collection of tracks, plus associations, that will
  // eventually be inserted into the event.

  std::unique_ptr<std::vector<recob::Track> > tracks(new std::vector<recob::Track>);
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > th_assn(new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > tsp_assn(new art::Assns<recob::Track, recob::SpacePoint>);
  std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> >  pfPartTrack_assns(  new art::Assns<recob::PFParticle, recob::Track>);

  // Make a collection of space points, plus associations, that will
  // be inserted into the event.

  std::unique_ptr<std::vector<recob::SpacePoint> > spts(new std::vector<recob::SpacePoint>);
  std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> > sph_assn(new art::Assns<recob::SpacePoint, recob::Hit>);

  // Make associations between PFParticles and the tracks they create
  // To facilitate this we'll recover the handle to the PFParticle collection - whether it exists or not
  art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
  evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle);

  // Make a map of a collection of KGTracks where we will save our results.
  // The first element of the map is the index of the corresponding hit colleciton.

  std::map<size_t, std::deque<KGTrack> > kalman_tracks_map;

  // Get Services.

  art::ServiceHandle<geo::Geometry> geom;

  // Reset space point algorithm.

  fSpacePointAlg.clearHitMap();

  // Get Hits.
  // There are three modes of operation:
  // 1.  Clustered hits (produces one hit collection).
  // 2.  PFParticle hits (products one hit collection for each PFParticle).
  // 3.  All hits (produces one hit collection).

  std::vector<art::PtrVector<recob::Hit> > hit_collections;

  // In the event we are making PFParticle/Track associations, we'll need the following map.
  // Hit collection index -> pfParticle index.
  std::map<size_t, size_t> hitColToPfPartMap;

  if(fUseClusterHits) {

    // Make one empty hit colleciton.

    hit_collections.emplace_back();
    art::PtrVector<recob::Hit>& hits = hit_collections.back();

    // Get clusters.

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    // Get hits from all clusters.
    art::FindManyP<recob::Hit> fm(clusterh, evt, fClusterModuleLabel);

    if(clusterh.isValid()) {
      int nclus = clusterh->size();

      for(int i = 0; i < nclus; ++i) {
	art::Ptr<recob::Cluster> pclus(clusterh, i);
	std::vector< art::Ptr<recob::Hit> > clushits = fm.at(i);
	int nhits = clushits.size();
	hits.reserve(hits.size() + nhits);

	for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = clushits.begin();
	    ihit != clushits.end(); ++ihit) {
	  hits.push_back(*ihit);
	}
      }
    }
  }
  else if(fUsePFParticleHits) {

    // Our program is to drive the track creation/fitting off the PFParticles in the data store
    // We'll use the hits associated to the PFParticles for each track - and only those hits.
    // Without a valid collection of PFParticles there is nothing to do here
    if (!pfParticleHandle.isValid()) return;
    
    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    evt.getByLabel(fPFParticleModuleLabel, clusterHandle);
    
    // If there are no clusters then something is really wrong
    if (!clusterHandle.isValid()) return;
    
    // Recover the collection of associations between PFParticles and clusters, this will
    // be the mechanism by which we actually deal with clusters
    art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleModuleLabel);
    
    // Likewise, recover the collection of associations to hits
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fPFParticleModuleLabel);
    
    // While PFParticle describes a hierarchal structure, for now we simply loop over the collection
    for(size_t partIdx = 0; partIdx < pfParticleHandle->size(); partIdx++) {
      art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, partIdx);

      // Add a new empty hit collection.
      hit_collections.emplace_back();
      art::PtrVector<recob::Hit>& hits = hit_collections.back();
        
      // Fill this vector by looping over associated clusters and finding the hits associated to them
      std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(pfParticle->Self());
        
      // Keep track of PFParticle to cluster collection
      hitColToPfPartMap[partIdx] = pfParticle->Self();
        
      for(const auto& cluster : clusterVec) {
	std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster->ID());
	hits.insert(hits.end(), hitVec.begin(), hitVec.end());
      }
    }
  }
  else {

    // Make one empty hit colleciton.

    hit_collections.emplace_back();
    art::PtrVector<recob::Hit>& hits = hit_collections.back();

    // Get unclustered hits.

    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(fHitModuleLabel, hith);
    if(hith.isValid()) {
      int nhits = hith->size();
      hits.reserve(nhits);

      for(int i = 0; i < nhits; ++i)
	hits.push_back(art::Ptr<recob::Hit>(hith, i));
    }
  }

  // Loop over hit collections.
  size_t hit_collection_cntr(0);

  for(auto& hits: hit_collections) {

    // Recover the kalman tracks double ended queue
    std::deque<KGTrack>& kalman_tracks = kalman_tracks_map[hit_collection_cntr++];

    // The hit collection "hits" (just filled), initially containing all
    // hits, represents hits available for making tracks.  Now we will
    // fill a second hit collection called "seederhits", also initially
    // containing all hits, which will represent hits available for
    // making track seeds.  These collections are not necessarily the
    // same, since hits that are not suitable for seeds may still be
    // suitable for tracks.

    art::PtrVector<recob::Hit> seederhits = hits;

    // Start of loop.

    bool done = false;
    while(!done) {

      // Use remaining seederhits to make seeds.

      std::vector<art::PtrVector<recob::Hit> > hitsperseed;
      std::vector<recob::Seed> seeds;
      if(seederhits.size()>0)
	seeds = fSeedFinderAlg.GetSeedsFromUnSortedHits(seederhits, hitsperseed);
      assert(seeds.size() == hitsperseed.size());
    
      if(seeds.size() == 0) {

	// Quit loop if we didn't find any new seeds.

	done = true;
	break;
      }
      else {

	// Loop over seeds.

	std::vector<recob::Seed>::const_iterator sit = seeds.begin();
	std::vector<art::PtrVector<recob::Hit> >::const_iterator hpsit = hitsperseed.begin();
	for(;sit != seeds.end() && hpsit != hitsperseed.end(); ++sit, ++hpsit) {

	  const recob::Seed& seed = *sit;

	  // Chop a couple of hits off each end of the seed.
	  // Seems like seeds often end at delta rays, Michel electrons,
	  // or other pathologies.

	  int nchopmax = std::max(0, int((hpsit->size() - fMinSeedChopHits)/2));
	  int nchop = std::min(nchopmax, fMaxChopHits);
	  art::PtrVector<recob::Hit>::const_iterator itb = hpsit->begin();
	  art::PtrVector<recob::Hit>::const_iterator ite = hpsit->end();
	  itb += nchop;
	  ite -= nchop;
	  art::PtrVector<recob::Hit> seedhits;
	  seedhits.reserve(hpsit->size());
	  for(art::PtrVector<recob::Hit>::const_iterator it = itb; it != ite; ++it)
	    seedhits.push_back(*it);

	  // Filter hits used by (chopped) seed from hits available to make future seeds.
	  // No matter what, we will never use these hits for another seed.
	  // This eliminates the possibility of an infinite loop.

	  size_t initial_seederhits = seederhits.size();
	  FilterHits(seederhits, seedhits);

	  // Require that this seed be fully disjoint from existing tracks.

	  if(seedhits.size() + seederhits.size() == initial_seederhits) {

	    // use mf::LogDebug instead of LOG_DEBUG because we reuse it in many lines;
	    // insertions are protected by mf::isDebugEnabled()
	    mf::LogDebug log("Track3DKalmanHit");

	    // Convert seed into initial KTracks on surface located at seed point, 
	    // and normal to seed direction.

	    double xyz[3];
	    double dir[3];
	    double err[3];   // Dummy.
	    seed.GetPoint(xyz, err);
	    seed.GetDirection(dir, err);

	    std::shared_ptr<const Surface> psurf(new SurfXYZPlane(xyz[0], xyz[1], xyz[2],
								  dir[0], dir[1], dir[2]));
	    TrackVector vec(5);
	    vec(0) = 0.;
	    vec(1) = 0.;
	    vec(2) = 0.;
	    vec(3) = 0.;
	    vec(4) = (fInitialMomentum != 0. ? 1./fInitialMomentum : 2.);
	  
	    if (mf::isDebugEnabled()) {
	      log << "Seed found with " << seedhits.size() <<" hits.\n"
		  << "(x,y,z) = " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "\n"
		  << "(dx,dy,dz) = " << dir[0] << ", " << dir[1] << ", " << dir[2] << "\n";
	    } // if debug

	    // Cut on the seed slope dx/dz.

	    if(std::abs(dir[0]) >= fMinSeedSlope * std::abs(dir[2])) {

	      // Make one or two initial KTracks for forward and backward directions.
	      // Assume muon (pdgid = 13).

	      int pdg = 13;
	      std::vector<KTrack> initial_tracks;

	      // The build_all flag specifies whether we should attempt to make
	      // tracks from all initial tracks, or alternatively, whether we 
	      // should declare victory and quit after getting a successful
	      // track from one initial track.

	      bool build_all = fDoDedx;
	      int ninit = 2;
	      initial_tracks.reserve(ninit);
	      initial_tracks.push_back(KTrack(psurf, vec, Surface::FORWARD, pdg));
	      if(ninit > 1)
		initial_tracks.push_back(KTrack(psurf, vec, Surface::BACKWARD, pdg));

	      // Loop over initial tracks.

	      int ntracks = kalman_tracks.size();   // Remember original track count.

	      for(std::vector<KTrack>::const_iterator itrk = initial_tracks.begin();
		  itrk != initial_tracks.end(); ++itrk) {
		const KTrack& trk = *itrk;

		// Fill hit container with current seed hits.

		KHitContainerWireX seedcont;
		seedcont.fill(seedhits, -1);

		// Set the preferred plane to be the one with the most hits.

		unsigned int prefplane = seedcont.getPreferredPlane();
		fKFAlg.setPlane(prefplane);
		if (mf::isDebugEnabled())
		  log << "Preferred plane = " << prefplane << "\n";

		// Build and smooth seed track.

		KGTrack trg0;
		bool ok = fKFAlg.buildTrack(trk, trg0, fProp, Propagator::FORWARD, seedcont);
		if(ok) {
		  KGTrack trg1;
		  ok = fKFAlg.smoothTrack(trg0, &trg1, fProp);
		  if(ok) {

		    // Now we have the seed track in the form of a KGTrack.
		    // Do additional quality cuts.

		    size_t n = trg1.numHits();
		    ok = (int(n) >= fMinSeedHits &&
			  trg0.startTrack().getChisq() <= n * fMaxSeedChiDF &&
			  trg0.endTrack().getChisq() <= n * fMaxSeedChiDF &&
			  trg1.startTrack().getChisq() <= n * fMaxSeedChiDF &&
			  trg1.endTrack().getChisq() <= n * fMaxSeedChiDF);
		    double mom0[3];
		    double mom1[3];
		    trg0.startTrack().getMomentum(mom0);
		    trg0.endTrack().getMomentum(mom1);
		    double dxdz0 = mom0[0] / mom0[2];
		    double dxdz1 = mom1[0] / mom1[2];
		    ok = ok && (std::abs(dxdz0) > fMinSeedSlope &&
				std::abs(dxdz1) > fMinSeedSlope);
		    if(ok) {

		      // Make a copy of the original hit collection of all
		      // available track hits.

		      art::PtrVector<recob::Hit> trackhits = hits;

		      // Do an extend + smooth loop here.
		      // Exit after two consecutive failures to extend (i.e. from each end),
		      // or if the iteration count reaches the maximum.

		      int niter = 6;
		      int nfail = 0;  // Number of consecutive extend fails.

		      for(int ix = 0; ok && ix < niter && nfail < 2; ++ix) {

			// Fill a collection of hits from the last good track
			// (initially the seed track).

			art::PtrVector<recob::Hit> goodhits;
			trg1.fillHits(goodhits);

			// Filter hits already on the track out of the available hits.

			FilterHits(trackhits, goodhits);

			// Fill hit container using filtered hits.

			KHitContainerWireX trackcont;
			trackcont.fill(trackhits, -1);

			// Extend the track.  It is not an error for the
			// extend operation to fail, meaning that no new hits
			// were added.
		      
			bool extendok = fKFAlg.extendTrack(trg1, fProp, trackcont);
			if(extendok)
			  nfail = 0;
			else
			  ++nfail;

			// Smooth the extended track, and make a new
			// unidirectionally fit track in the opposite
			// direction.

			KGTrack trg2;
			ok = fKFAlg.smoothTrack(trg1, &trg2, fProp);
			if(ok) {

			  // Skip momentum estimate for constant-momentum tracks.

			  if(fDoDedx) {
			    KETrack tremom;
			    bool pok = fKFAlg.fitMomentum(trg1, fProp, tremom);
			    if(pok)
			      fKFAlg.updateMomentum(tremom, fProp, trg2);
			  }
			  trg1 = trg2;
			}
		      }

		      // Do a final smooth.

		      if(ok) {
			ok = fKFAlg.smoothTrack(trg1, 0, fProp);
			if(ok) {

			  // Skip momentum estimate for constant-momentum tracks.

			  if(fDoDedx) {
			    KETrack tremom;
			    bool pok = fKFAlg.fitMomentum(trg1, fProp, tremom);
			    if(pok)
			      fKFAlg.updateMomentum(tremom, fProp, trg1);
			  }

			  // Save this track.

			  ++fNumTrack;
			  kalman_tracks.push_back(trg1);
			}
		      }
		    }
		  }
		}
		if (mf::isDebugEnabled())
		  log << (ok? "Find track succeeded.": "Find track failed.") << "\n";
		if(ok && !build_all)
		  break;	    
	      } // for initial track
	    
	      // Loop over newly added tracks and remove hits contained on
	      // these tracks from hits available for making additional
	      // tracks or track seeds.

	      for(unsigned int itrk = ntracks; itrk < kalman_tracks.size(); ++itrk) {
		const KGTrack& trg = kalman_tracks[itrk];
		art::PtrVector<recob::Hit> track_used_hits;
		trg.fillHits(track_used_hits);
		FilterHits(hits, track_used_hits);
		FilterHits(seederhits, track_used_hits);
	      }
	    }
	  }
	}
      }
    }
  }

  // Fill histograms.

  if(fHist) {

    // First loop over tracks.

    for(auto& kalman_tracks_pair : kalman_tracks_map) {
      std::deque<KGTrack>& kalman_tracks = kalman_tracks_pair.second;

      for(std::deque<KGTrack>::const_iterator k = kalman_tracks.begin();
	  k != kalman_tracks.end(); ++k) {
	const KGTrack& trg = *k;

	// Loop over measurements in this track.

	const std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();
	for(std::multimap<double, KHitTrack>::const_iterator ih = trackmap.begin();
	    ih != trackmap.end(); ++ih) {
	  const KHitTrack& trh = (*ih).second;
	  const std::shared_ptr<const KHitBase>& hit = trh.getHit();
	  double chisq = hit->getChisq();
	  fHIncChisq->Fill(chisq);
	  const KHit<1>* ph1 = dynamic_cast<const KHit<1>*>(&*hit);
	  if(ph1 != 0) {
	    double pull = ph1->getResVector()(0) / std::sqrt(ph1->getResError()(0, 0));
	    fHPull->Fill(pull);
	  }
	}
      }
    }
  }

  // Process Kalman filter tracks into persistent objects.

  size_t tracksSize(0);
  for(auto kalman_tracks_pair : kalman_tracks_map) {
    tracksSize += kalman_tracks_pair.second.size();
  }
  tracks->reserve(tracksSize);

  for(auto kalman_tracks_pair : kalman_tracks_map) {

    // Recover the kalman tracks double ended queue
    std::deque<KGTrack>& kalman_tracks = kalman_tracks_pair.second;
      
    size_t trackColStartIdx(tracks->size());
      
    for(std::deque<KGTrack>::const_iterator k = kalman_tracks.begin();
	k != kalman_tracks.end(); ++k) {
      const KGTrack& kalman_track = *k;

      // Add Track object to collection.

      tracks->push_back(recob::Track());
      kalman_track.fillTrack(tracks->back(), tracks->size() - 1);

      // Make Track to Hit associations.  

      art::PtrVector<recob::Hit> trhits;
      kalman_track.fillHits(trhits);
      util::CreateAssn(*this, evt, *tracks, trhits, *th_assn, tracks->size()-1);

      // Make space points from this track.

      int nspt = spts->size();
      fSpacePointAlg.fillSpacePoints(*spts, kalman_track.TrackMap());

      // Associate newly created space points with hits.
      // Also associate track with newly created space points.

      art::PtrVector<recob::SpacePoint> sptvec;

      // Loop over newly created space points.
      
      for(unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
	const recob::SpacePoint& spt = (*spts)[ispt];
	art::ProductID sptid = getProductID<std::vector<recob::SpacePoint> >(evt);
	art::Ptr<recob::SpacePoint> sptptr(sptid, ispt, evt.productGetter(sptid));
	sptvec.push_back(sptptr);

	// Make space point to hit associations.

	const art::PtrVector<recob::Hit>& sphits = 
	  fSpacePointAlg.getAssociatedHits(spt);
	util::CreateAssn(*this, evt, *spts, sphits, *sph_assn, ispt);
      }

      // Make track to space point associations.

      util::CreateAssn(*this, evt, *tracks, sptvec, *tsp_assn, tracks->size()-1);
    } // end of loop over a given collection

    if (fUsePFParticleHits) {
      size_t hitColIdx(kalman_tracks_pair.first);
      size_t pfPartIdx(hitColToPfPartMap[hitColIdx]);
      size_t trackColEndIdx(tracks->size());
          
      art::Ptr<recob::PFParticle> pfPartPtr(pfParticleHandle, pfPartIdx);
          
      // Create a vector of art ptrs to the tracks...
      std::vector<art::Ptr<recob::Track> > trackVec;
          
      for(size_t idx = trackColStartIdx; idx != trackColEndIdx; idx++) {
	art::ProductID trackId = getProductID<std::vector<recob::Track> >(evt);
	art::Ptr<recob::Track> trackPtr(trackId, idx, evt.productGetter(trackId));
	trackVec.push_back(trackPtr);
      }
          
      util::CreateAssn(*this, evt, pfPartPtr, trackVec, *pfPartTrack_assns);
    }
  }

  // Add tracks and associations to event.

  evt.put(std::move(tracks));
  evt.put(std::move(spts));
  evt.put(std::move(th_assn));
  evt.put(std::move(tsp_assn));
  evt.put(std::move(sph_assn));
  evt.put(std::move(pfPartTrack_assns));
}

//----------------------------------------------------------------------------
/// End job method.
void trkf::Track3DKalmanHit::endJob()
{
  mf::LogInfo("Track3DKalmanHit") 
    << "Track3DKalmanHit statistics:\n"
    << "  Number of events = " << fNumEvent << "\n"
    << "  Number of tracks created = " << fNumTrack;
}
