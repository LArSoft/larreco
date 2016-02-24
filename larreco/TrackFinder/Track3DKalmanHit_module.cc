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
// UsePFParticleSeeds - If true, use seeds associated with PFParticles.
// HitModuleLabel     - Module label for unclustered Hits.
// ClusterModuleLabel - Module label for Clusters.
// PFParticleModuleLabel - Module label for PFParticles.
// StoreNPPlane       - Store nonpreferred planes trajectory points.
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// DoDedx             - Global dE/dx enable flag.
// SelfSeed           - Self seed flag.
// LineSurface        - Hits on line surfaces (true) or plane surfaces (false).
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
// 4.  UsePFParticleSeeds has no effect unless UsePFParticleHits is true.
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Seed.h"
#include "larreco/RecoAlg/Track3DKalmanHitAlg.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "lardata/Utilities/AssociationUtil.h"


#include "TH1F.h"

namespace trkf {
   
   class Track3DKalmanHit : public art::EDProducer {
   public:
      
      // Copnstructors, destructor.
      
      explicit Track3DKalmanHit(fhicl::ParameterSet const & pset);
      virtual ~Track3DKalmanHit();
      
      // Overrides.
      // Put override right next to each function
      virtual void reconfigure(fhicl::ParameterSet const & pset) override;
      virtual void produce(art::Event & e) override;
      virtual void beginJob() override;
      virtual void endJob() override;
      
   private:
      void fillHistograms(std::list<LocalKalmanStruct>& LocalKalmanStructList);
      //Member functions that depend on art::event and use art::Assns
      void getInputfromevent(const art::Event &evt,
                             std::list<LocalKalmanStruct> &LocalKalmanStructList);
      void getClusteredHits(const art::Event & evt,
                            art::PtrVector<recob::Hit>& hits) const;
      void getPFParticleHits(const art::Event & evt,
                             std::list<LocalKalmanStruct> & LocalKalmanStructs) const;
      void getAllHits(const art::Event & evt,
                      art::PtrVector<recob::Hit>& hits) const;
      void persistObjects(const art::Event &evt,
                          std::list<LocalKalmanStruct> const &LocalKalmanStructList,
                          std::vector<recob::Track> &tracks,
                          std::vector<recob::SpacePoint> &spts,
                          art::Assns<recob::Track, recob::Hit> &th_assn,
                          art::Assns<recob::Track, recob::SpacePoint> &tsp_assn,
                          art::Assns<recob::SpacePoint, recob::Hit> &sph_assn,
                          art::Assns<recob::PFParticle, recob::Track> &pfPartTrack_assns);
      
      // Fcl parameters.
      
      bool fHist;                         ///< Make histograms.
      bool fUseClusterHits;               ///< Use clustered hits as input.
      bool fUsePFParticleHits;            ///< Use PFParticle hits as input.
      bool fUsePFParticleSeeds;           ///< Use PFParticle seeds.
      std::string fHitModuleLabel;        ///< Unclustered Hits.
      std::string fClusterModuleLabel;    ///< Clustered Hits.
      std::string fPFParticleModuleLabel; ///< PFParticle label.
      bool fStoreNPPlane;                 ///< Store nonpreferred planes trajectory points.
      Track3DKalmanHitAlg fTKHAlg;
      SpacePointAlg fSpacePointAlg;       ///< Space point algorithm.
      
      // Histograms.
      
      TH1F* fHIncChisq;   ///< Incremental chisquare.
      TH1F* fHPull;       ///< Hit pull.
      // Statistics.
      
      int fNumEvent;    ///< Number of events seen.
      
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
fUsePFParticleSeeds(false),
fStoreNPPlane(true),
fTKHAlg(pset.get<fhicl::ParameterSet>("Track3DKalmanHitAlg")),
fSpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
fHIncChisq(0),
fHPull(0),
fNumEvent(0)
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
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// p - Fcl parameter set.
///
// SS: talk to Kyle about using parameter set validation
void trkf::Track3DKalmanHit::reconfigure(fhicl::ParameterSet const & pset)
{
   fHist = pset.get<bool>("Hist");
   fTKHAlg.reconfigure(pset.get<fhicl::ParameterSet>("Track3DKalmanHitAlg"));
   fSpacePointAlg.reconfigure(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
   fUseClusterHits = pset.get<bool>("UseClusterHits");
   fUsePFParticleHits = pset.get<bool>("UsePFParticleHits");
   fUsePFParticleSeeds = pset.get<bool>("UsePFParticleSeeds");
   fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
   fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
   fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
   fStoreNPPlane = pset.get<bool>("StoreNPPlane");
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
   std::cout << "Track3DKalmanHit_module: in Produce" << "\n";
   // Make a collection of tracks, plus associations, that will
   // eventually be inserted into the event.
   //use auto and make_unique
   auto tracks = std::make_unique<std::vector<recob::Track>>();
   auto th_assn = std::make_unique<art::Assns<recob::Track, recob::Hit>>();
   auto tsp_assn = std::make_unique<art::Assns<recob::Track, recob::SpacePoint>>();
   auto pfPartTrack_assns = std::make_unique<art::Assns<recob::PFParticle, recob::Track>>();
   
   // Make a collection of space points, plus associations, that will
   // be inserted into the event.
   
   auto spts = std::make_unique<std::vector<recob::SpacePoint>>();
   auto sph_assn = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();
   
   // Reset space point algorithm.
   fSpacePointAlg.clearHitMap();
   
   // Get Hits.
   std::list<LocalKalmanStruct> LocalKalmanStructList;
   getInputfromevent(evt, LocalKalmanStructList);
   //SS: LocalKalmanStruct.hits holds the input hits, .tracks will have the redulting Kalmantracks
   fTKHAlg.generateKalmantracks(LocalKalmanStructList);
   
   // Fill histograms.
   if(fHist) {
      fillHistograms(LocalKalmanStructList);
   }
   
   // Process Kalman filter tracks into persistent objects.
   persistObjects(evt, LocalKalmanStructList, *tracks, *spts, *th_assn, *tsp_assn, *sph_assn, *pfPartTrack_assns);
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
   << "  Number of events = " << fNumEvent << "\n";
   //<< "  Number of tracks created = " << fNumTrack;
}


//----------------------------------------------------------------------------
// There are three modes of operation:
// 1.  Clustered hits (produces one hit collection).
// 2.  PFParticle hits (products one hit collection for each PFParticle).
// 3.  All hits (produces one hit collection).

void trkf::Track3DKalmanHit::getInputfromevent(const art::Event &evt,
                                               std::list<LocalKalmanStruct>& LocalKalmanStructList){
   if (fUsePFParticleHits) {
      getPFParticleHits(evt, LocalKalmanStructList);
   }
   else {
      LocalKalmanStructList.emplace_back();
      LocalKalmanStruct& local_kalman_struct = LocalKalmanStructList.back();
      art::PtrVector<recob::Hit>& hits = local_kalman_struct.hits;
      if(fUseClusterHits) {
         getClusteredHits(evt, hits);
      }
      else {
         getAllHits(evt, hits);
      }
   }
}


//----------------------------------------------------------------------------
/// Fill a collection using clustered hits

// SS: after method extraction, now I am working on replace temp with query
void trkf::Track3DKalmanHit::getClusteredHits(const art::Event &evt,
                                              art::PtrVector<recob::Hit>& hits) const{
   // Get clusters.
   //SS: Can we use getValidHandle?
   art::Handle< std::vector<recob::Cluster> > clusterh;
   evt.getByLabel(fClusterModuleLabel, clusterh);
   if (!clusterh.isValid()) return;
   
   // Get hits from all clusters.
   art::FindManyP<recob::Hit> hitsbycluster(clusterh, evt, fClusterModuleLabel);
   
   for(size_t i = 0; i < clusterh->size(); ++i) {
      std::vector< art::Ptr<recob::Hit> > clushits = hitsbycluster.at(i);
      hits.insert(hits.end(), clushits.begin(), clushits.end());
   }
}

//----------------------------------------------------------------------------
/// If both UseClusteredHits and UsePFParticles is false use this method to fill in hits

void trkf::Track3DKalmanHit::getAllHits(const art::Event &evt,
                                        art::PtrVector<recob::Hit>& hits) const{
   // Get unclustered hits.
   art::Handle< std::vector<recob::Hit> > hith;
   evt.getByLabel(fHitModuleLabel, hith);
   if(!hith.isValid()) return;
   size_t nhits = hith->size();
   hits.reserve(nhits);
   
   for(size_t i = 0; i < nhits; ++i) {
      hits.push_back(art::Ptr<recob::Hit>(hith, i));
   }
   
}

//----------------------------------------------------------------------------
/// If UsePFParticles is true use this method to fill in hits

void trkf::Track3DKalmanHit::getPFParticleHits(const art::Event &evt,
                                               std::list<LocalKalmanStruct> & localcoll) const{
   
   // Our program is to drive the track creation/fitting off the PFParticles in the data store
   // We'll use the hits associated to the PFParticles for each track - and only those hits.
   // Without a valid collection of PFParticles there is nothing to do here
   // We need a handle to the collection of clusters in the data store so we can
   // handle associations to hits.
   art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
   evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle);
   if (!pfParticleHandle.isValid()) return;
   
   //std::cout << "Track3DKalmanHit: pfParticleHandle size" << pfParticleHandle->size() << "\n";
   art::Handle<std::vector<recob::Cluster> > clusterHandle;
   evt.getByLabel(fClusterModuleLabel, clusterHandle);
   
   // If there are no clusters then something is really wrong
   if (!clusterHandle.isValid()) return;
   //{
   // Recover the collection of associations between PFParticles and clusters, this will
   // be the mechanism by which we actually deal with clusters
   art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleModuleLabel);
   
   // Associations to seeds.
   art::FindManyP<recob::Seed> seedAssns(pfParticleHandle, evt, fPFParticleModuleLabel);
   
   // Likewise, recover the collection of associations to hits
   art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fClusterModuleLabel);
   
   // While PFParticle describes a hierarchal structure, for now we simply loop over the collection
   for(size_t partIdx = 0; partIdx < pfParticleHandle->size(); partIdx++) {
      
      // Add a new empty hit collection.
      localcoll.emplace_back();
      LocalKalmanStruct& local_kalman_struct = localcoll.back();
      local_kalman_struct.pfPartPtr = art::Ptr<recob::PFParticle>(pfParticleHandle, partIdx);
      art::PtrVector<recob::Hit>& hits = local_kalman_struct.hits;
      
      // Fill this hit vector by looping over associated clusters and finding the
      // hits associated to them
      std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(partIdx);
      
      for(const auto& cluster : clusterVec) {
         std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster.key());
         hits.insert(hits.end(), hitVec.begin(), hitVec.end());
      }
      
      // If requested, fill associated seeds.
      if(!fUsePFParticleSeeds) continue;
      art::PtrVector<recob::Seed>& seeds = local_kalman_struct.seeds;
      std::vector<art::Ptr<recob::Seed> > seedVec = seedAssns.at(partIdx);
      seeds.insert(seeds.end(), seedVec.begin(), seedVec.end());
      art::FindManyP<recob::Hit> seedHitAssns(seedVec, evt, fPFParticleModuleLabel);
      // std::cout << "Track3DKalmanHit: seedHitAssns " << seedHitAssns.size() << ", "<< seedVec.size() << "\n";
      for(size_t seedIdx = 0; seedIdx < seedVec.size(); ++seedIdx) {
         std::vector<art::Ptr<recob::Hit> > seedHitVec;
         //SS: why seedIdx can have an invalid value?
         try {
            seedHitVec = seedHitAssns.at(seedIdx);
         }
         catch(art::Exception x) {
            seedHitVec.clear();
         }
         local_kalman_struct.seedhits.emplace_back();
         art::PtrVector<recob::Hit>& seedhits = local_kalman_struct.seedhits.back();
         seedhits.insert(seedhits.end(), seedHitVec.begin(), seedHitVec.end());
      }
   }
   
   // std::cout << "Track3DKalmanHit: end localcoll size " << localcoll.size() << "\n";
}


//-------------------------------------------------------------------------------------
void trkf::Track3DKalmanHit::persistObjects(const art::Event &evt,
                                            std::list<LocalKalmanStruct> const &LocalKalmanStructList,
                                            std::vector<recob::Track> &tracks,
                                            std::vector<recob::SpacePoint> &spts,
                                            art::Assns<recob::Track, recob::Hit> &th_assn,
                                            art::Assns<recob::Track, recob::SpacePoint> &tsp_assn,
                                            art::Assns<recob::SpacePoint, recob::Hit> &sph_assn,
                                            art::Assns<recob::PFParticle, recob::Track> &pfPartTrack_assns)
{
   size_t tracksSize(0);
   for(const auto& local_kalman_struct : LocalKalmanStructList) {
      tracksSize += local_kalman_struct.tracks.size();
   }
   tracks.reserve(tracksSize);
   
   std::cout << "persist Objects " << LocalKalmanStructList.size() <<  "\n";
   
   auto const tid = getProductID<std::vector<recob::Track> >(evt);
   auto const tidgetter = evt.productGetter(tid);
   
   auto const spacepointId = getProductID<std::vector<recob::SpacePoint> >(evt);
   auto const getter = evt.productGetter(spacepointId);
   
   for(auto& local_kalman_struct : LocalKalmanStructList) {
      // Recover the kalman tracks double ended queue
      const std::deque<KGTrack>& kalman_tracks = local_kalman_struct.tracks;
      //      std::cout << "persist Objects: kalman tracks" << kalman_tracks.size() <<  "\n";
      
      for(auto const& kalman_track:kalman_tracks) {
         
         // Add Track object to collection.
         recob::Track track;
         kalman_track.fillTrack(track, tracks.size(), fStoreNPPlane);
         if(track.NumberTrajectoryPoints() < 2) {
            continue;
         }
         tracks.emplace_back(std::move(track));
         // SS: tracks->size() does not change after this point in each iteration
         
         //fill hits from this track
         art::PtrVector<recob::Hit> trhits;
         kalman_track.fillHits(trhits);
         
         // Make space points from this track.
         auto nspt = spts.size();
         fSpacePointAlg.fillSpacePoints(spts, kalman_track.TrackMap());
         
         
         std::vector<art::Ptr<recob::SpacePoint>> sptvec;
         for(auto ispt = nspt; ispt < spts.size(); ++ispt) {
            sptvec.emplace_back(spacepointId, ispt, getter);
            // Associate newly created space points with hits.
            // Make space point to hit associations.
            const auto& sphits = fSpacePointAlg.getAssociatedHits((spts)[ispt]);
            for(auto const& sphit: sphits) {
               sph_assn.addSingle(sptvec.back(), sphit);
            }
         }
         
         art::Ptr<recob::Track> aptr(tid, tracks.size()-1, tidgetter);
         
         // Make Track to Hit associations.
         for (auto const& trhit: trhits) {
            th_assn.addSingle(aptr, trhit);
         }
         
         // Make track to space point associations
         for (auto const& spt: sptvec) {
            tsp_assn.addSingle(aptr, spt);
         }
         
         // Optionally fill track-to-PFParticle associations.
         if (fUsePFParticleHits) {
            pfPartTrack_assns.addSingle(local_kalman_struct.pfPartPtr, aptr);
         }
      } // end of loop over a given collection
   }
}

//----------------------------------------------------------------------------
/// Fill Histograms method
//fHPull and fHIncChisq are private data members of the class Track3DKalmanHit

void trkf::Track3DKalmanHit::fillHistograms(std::list<LocalKalmanStruct>& LocalKalmanStructList)
{
   for(const auto& local_kalman_struct : LocalKalmanStructList) {
      const std::deque<KGTrack>& kalman_tracks = local_kalman_struct.tracks;
      
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
      //  fHIncChisq->Print("all");
      //  fHPull->Print("all");
   }
}



