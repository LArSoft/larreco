////////////////////////////////////////////////////////////////////////
// Class:       TrackKalmanCheater
// Module Type: producer
// File:        TrackKalmanCheater.h
//
// This class produces RecoBase/Track objects using KalmanFilterService.
// MC truth information is used to associate Hits used as input to
// the Kalman filter.
//
// Configuration parameters:
//
// Hist               - Histogram flag (generate histograms if true).
// UseClusterHits     - Use clustered hits if true, use all hits if false.
// HitModuleLabel     - Module label for unclustered Hits.
// ClusterModuleLabel - Module label for Clusters.
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// KalmanFilterAlg    - Parameter set for KalmanFilterAlg.
// SpacePointAlg      - Parmaeter set for space points.
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Hit.h"
#include "MCCheater/BackTracker.h"
#include "RecoAlg/KalmanFilterAlg.h"
#include "RecoObjects/KHitContainerWireX.h"
#include "RecoObjects/SurfYZPlane.h"
#include "RecoObjects/PropYZPlane.h"
#include "RecoObjects/KHit.h"
#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCParticle.h"

#include "TH1F.h"

namespace trkf {

  class Propagator;

  class TrackKalmanCheater : public art::EDProducer {
  public:

    // Copnstructors, destructor.

    explicit TrackKalmanCheater(fhicl::ParameterSet const & pset);
    virtual ~TrackKalmanCheater();

    // Overrides.

    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

  private:

    // Fcl parameters.

    bool fHist;                        ///< Make histograms.
    KalmanFilterAlg fKFAlg;            ///< Kalman filter algorithm.
    SpacePointAlg fSpacePointAlg;      ///< Space point algorithm.
    bool fUseClusterHits;              ///< Use cluster hits or all hits?
    std::string fHitModuleLabel;       ///< Unclustered Hits.
    std::string fClusterModuleLabel;   ///< Clustered Hits.
    double fMaxTcut;                   ///< Maximum delta ray energy in MeV for restricted dE/dx.

    /// Propagator.
    const Propagator* fProp;

    // Histograms.

    TH1F* fHIncChisq;   ///< Incremental chisquare.
    TH1F* fHPull;       ///< Hit pull.

    // Statistics.

    int fNumEvent;    ///< Number of events seen.
    int fNumTrack;    ///< Number of tracks produced.

  };

  DEFINE_ART_MODULE(TrackKalmanCheater)

} // namespace trkf

//------------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// p - Fcl parameters.
///
trkf::TrackKalmanCheater::TrackKalmanCheater(fhicl::ParameterSet const & pset)
  : fHist(false)
  , fKFAlg(pset.get<fhicl::ParameterSet>("KalmanFilterAlg"))
  , fSpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"))
  , fUseClusterHits(false)
  , fMaxTcut(0.)
  , fProp(0)
  , fHIncChisq(0)
  , fHPull(0)
  , fNumEvent(0)
  , fNumTrack(0)
{
  reconfigure(pset);
  produces<std::vector<recob::Track>                   >();
  produces<std::vector<recob::SpacePoint>              >();
  produces<art::Assns<recob::Track, recob::Hit>        >();
  produces<art::Assns<recob::Track, recob::SpacePoint> >();
  produces<art::Assns<recob::SpacePoint, recob::Hit>   >();

  // Report.

  mf::LogInfo("TrackKalmanCheater") 
    << "TrackKalmanCheater configured with the following parameters:\n"
    << "  UseClusterHits = " << fUseClusterHits << "\n"
    << "  HitModuleLabel = " << fHitModuleLabel << "\n"
    << "  ClusterModuleLabel = " << fClusterModuleLabel;
}

//------------------------------------------------------------------------------
/// Destructor.
trkf::TrackKalmanCheater::~TrackKalmanCheater()
{}

//------------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// p - Fcl parameter set.
///
void trkf::TrackKalmanCheater::reconfigure(fhicl::ParameterSet const & pset)
{
  fHist = pset.get<bool>("Hist");
  fKFAlg.reconfigure(pset.get<fhicl::ParameterSet>("KalmanFilterAlg"));
  fSpacePointAlg.reconfigure(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
  fUseClusterHits = pset.get<bool>("UseClusterHits");
  fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
  fMaxTcut = pset.get<double>("MaxTcut");
  if(fProp != 0)
    delete fProp;
  fProp = new PropYZPlane(fMaxTcut);
}

//------------------------------------------------------------------------------
/// Begin job method.
void trkf::TrackKalmanCheater::beginJob()
{
  if(fHist) {

    // Book histograms.

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("hitkalman", "TrackKalmanCheater histograms");

    fHIncChisq = dir.make<TH1F>("IncChisq", "Incremental Chisquare", 100, 0., 20.);
    fHPull = dir.make<TH1F>("Pull", "Hit Pull", 100, -10., 10.);
  }
}

//------------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// e - Art event.
///
/// This method extracts Hit from the event and produces and adds
/// Track objects.
///
void trkf::TrackKalmanCheater::produce(art::Event & evt)
{
  ++fNumEvent;

  // Make a collection of tracks, plus associations, that will
  // eventually be inserted into the event.

  std::unique_ptr<std::vector<recob::Track> > tracks(new std::vector<recob::Track>);
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > th_assn(new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > tsp_assn(new art::Assns<recob::Track, recob::SpacePoint>);

  // Make a collection of space points, plus associations, that will
  // be inserted into the event.

  std::unique_ptr<std::vector<recob::SpacePoint> > spts(new std::vector<recob::SpacePoint>);
  std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> > sph_assn(new art::Assns<recob::SpacePoint, recob::Hit>);

  // Make a collection of KGTracks where we will save our results.

  std::deque<KGTrack> kalman_tracks;

  // Get Services.

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTracker> bt;

  // Reset space point algorithm.

  fSpacePointAlg.clearHitMap();


  // Get Hits.

  art::PtrVector<recob::Hit> hits;

  if(fUseClusterHits) {

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
  else {

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

  // Sort hits into separate PtrVectors based on track id.

  std::map<int, art::PtrVector<recob::Hit> > hitmap;

  // Loop over hits.

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {
    //const recob::Hit& hit = **ihit;

    // Get track ids for this hit.

    std::vector<cheat::TrackIDE> tids = bt->HitToTrackID(*ihit);

    // Loop over track ids.

    for(std::vector<cheat::TrackIDE>::const_iterator itid = tids.begin();
	itid != tids.end(); ++itid) {
      int trackID = itid->trackID;

      // Add hit to PtrVector corresponding to this track id.

      hitmap[trackID].push_back(*ihit);
    }
  }

  // Extract geant mc particles.

  sim::ParticleList plist = bt->ParticleList();

  // Loop over geant particles.

  {
    mf::LogDebug log("TrackKalmanCheater");
    for(sim::ParticleList::const_iterator ipart = plist.begin();
	ipart != plist.end(); ++ipart) {
      const simb::MCParticle* part = (*ipart).second;
      assert(part != 0);
      int pdg = part->PdgCode();

      // Ignore everything except stable charged nonshowering particles.

      int apdg = std::abs(pdg);
      if(apdg == 13 ||     // Muon
	 apdg == 211 ||    // Charged pion
	 apdg == 321 ||    // Charged kaon
	 apdg == 2212) {   // (Anti)proton

	int trackid = part->TrackId();
	int stat = part->StatusCode();
	int nhit = 0;
	if(hitmap.count(trackid) != 0)
	  nhit = hitmap[trackid].size();
	log << "Trackid=" << trackid 
	    << ", pdgid=" << pdg 
	    << ", status=" << stat
	    << ", NHit=" << nhit << "\n";
	const TLorentzVector& pos = part->Position();
	double x = pos.X();
	double y = pos.Y();
	double z = pos.Z();
	const TLorentzVector& mom = part->Momentum();
	double px = mom.Px();
	double py = mom.Py();
	double pz = mom.Pz();
	double p = std::sqrt(px*px + py*py + pz*pz);
	log << "  x = " << pos.X()
	    << ", y = " << pos.Y()
	    << ", z = " << pos.Z() << "\n"
	    << "  px= " << mom.Px()
	    << ", py= " << mom.Py()
	    << ", pz= " << mom.Pz() << "\n";

	if(nhit > 0 && pz != 0.) {
	  const art::PtrVector<recob::Hit>& trackhits = hitmap[trackid];
	  assert(trackhits.size() > 0);

	  // Make a seed track (KTrack).

	  std::shared_ptr<const Surface> psurf(new SurfYZPlane(y, z, 0.));
	  TrackVector vec(5);
	  vec(0) = x;
	  vec(1) = 0.;
	  vec(2) = px / pz;
	  vec(3) = py / pz;
	  vec(4) = 1. / p;
	  Surface::TrackDirection dir = (pz > 0. ? Surface::FORWARD : Surface::BACKWARD);
	  KTrack trk(psurf, vec, dir, pdg);

	  // Fill KHitContainer with hits.

	  KHitContainerWireX cont;
	  cont.fill(trackhits, -1);

	  // Count hits in each plane.  Set the preferred plane to be
	  // the one with the most hits.

	  std::vector<unsigned int> planehits(3, 0);
	  for(art::PtrVector<recob::Hit>::const_iterator ih = trackhits.begin();
	      ih != trackhits.end(); ++ih) {
	    const recob::Hit& hit = **ih;
	    unsigned int plane = hit.WireID().Plane;

	    assert(plane >= 0 && plane < planehits.size());
	    ++planehits[plane];
	  }
	  unsigned int prefplane = 0;
	  for(unsigned int i=0; i<planehits.size(); ++i) {
	    if(planehits[i] > planehits[prefplane])
	      prefplane = i;
	  }
	  log << "Preferred plane = " << prefplane << "\n";

	  // Build and smooth track.

	  KGTrack trg;
	  fKFAlg.setPlane(prefplane);
	  bool ok = fKFAlg.buildTrack(trk, trg, fProp, Propagator::FORWARD, cont);
	  if(ok) {
	    ok = fKFAlg.smoothTrackIter(5, trg, fProp);
	    if(ok) {
	      KETrack tremom;
	      bool pok = fKFAlg.fitMomentum(trg, fProp, tremom);
	      if(pok)
		fKFAlg.updateMomentum(tremom, fProp, trg);
	      ++fNumTrack;
	      kalman_tracks.push_back(trg);
	    }
	  }
	  if(ok)
	    log << "Build track succeeded.\n";
	  else
	    log << "Build track failed.\n";
 	}
      }
    }
  }

  // Fill histograms.

  if(fHist) {

    // First loop over tracks.

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

  // Process Kalman filter tracks into persistent objects.

  tracks->reserve(kalman_tracks.size());
  for(std::deque<KGTrack>::const_iterator k = kalman_tracks.begin();
      k != kalman_tracks.end(); ++k) {
    const KGTrack& kalman_track = *k;

    // Add Track object to collection.

    tracks->push_back(recob::Track());
    kalman_track.fillTrack(tracks->back(), tracks->size() - 1);

    // Make Track to Hit associations.  

    art::PtrVector<recob::Hit> trhits;
    kalman_track.fillHits(hits);
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
  }

  // Add tracks and associations to event.

  evt.put(std::move(tracks));
  evt.put(std::move(spts));
  evt.put(std::move(th_assn));
  evt.put(std::move(tsp_assn));
  evt.put(std::move(sph_assn));
}

//------------------------------------------------------------------------------
/// End job method.
void trkf::TrackKalmanCheater::endJob()
{
  mf::LogInfo("TrackKalmanCheater") 
    << "TrackKalmanCheater statistics:\n"
    << "  Number of events = " << fNumEvent << "\n"
    << "  Number of tracks created = " << fNumTrack;
}
