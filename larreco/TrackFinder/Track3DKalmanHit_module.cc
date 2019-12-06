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
// Hist                  - Histogram flag (generate histograms if true).
// UseClusterHits        - Use clustered hits if true.
// UsePFParticleHits     - Use PFParticle hits if true.
// UsePFParticleSeeds    - If true, use seeds associated with PFParticles.
// HitModuleLabel        - Module label for unclustered Hits.
// ClusterModuleLabel    - Module label for Clusters.
// PFParticleModuleLabel - Module label for PFParticles.
// Track3DKalmanHitAlg   - Algorithm class for Track3DKalmanHit module
// SpacePointAlg         - Parmaeter set for space points.
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
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/RecoObjects/KHit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larreco/RecoAlg/Track3DKalmanHitAlg.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "larreco/RecoAlg/Track3DKalmanHit.h"

#include "TH1F.h"

namespace trkf {

   class Track3DKalmanHit : public art::EDProducer {
   public:
      // Copnstructors, destructor.
      explicit Track3DKalmanHit(fhicl::ParameterSet const & pset);

   private:
      void produce(art::Event & e) override;
      void beginJob() override;
      void endJob() override;

      //Member functions that depend on art::event and use art::Assns
      //note that return types are move aware
      void prepareForInput();
      KalmanInputs getInput(const art::Event &evt) const;
      Hits getClusteredHits(const art::Event &evt) const;
      KalmanInputs getPFParticleStuff(const art::Event &evt) const;
      Hits getAllHits(const art::Event &evt) const;
      void createOutputs(const art::Event &evt,
                         const std::vector<KalmanOutput> &outputs,
                         const KalmanInputs &inputs,
                         std::vector<recob::Track> &tracks,
                         std::vector<recob::SpacePoint> &spts,
                         art::Assns<recob::Track, recob::Hit> &th_assn,
                         art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> &thm_assn,
                         art::Assns<recob::Track, recob::SpacePoint> &tsp_assn,
                         art::Assns<recob::SpacePoint, recob::Hit> &sph_assn,
                         art::Assns<recob::PFParticle, recob::Track> &pfPartTrack_assns);
      void fillHistograms(std::vector<KalmanOutput>& outputs);


      // Fcl parameters.
      bool fHist;                         ///< Make histograms.
      bool fUseClusterHits;               ///< Use clustered hits as input.
      bool fUsePFParticleHits;            ///< Use PFParticle hits as input.
      bool fUsePFParticleSeeds;           ///< Use PFParticle seeds.
      std::string fHitModuleLabel;        ///< Unclustered Hits.
      std::string fClusterModuleLabel;    ///< Clustered Hits.
      std::string fPFParticleModuleLabel; ///< PFParticle label.
      Track3DKalmanHitAlg fTKHAlg;        ///< Track3DKalmanHit algorithm.
      SpacePointAlg fSpacePointAlg;       ///< Space point algorithm.

      // Histograms.
      TH1F* fHIncChisq;                   ///< Incremental chisquare.
      TH1F* fHPull;                      ///< Hit pull.

      // Statistics.
      int fNumEvent;                    ///< Number of events seen.
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
EDProducer{pset},
fHist(false),
fUseClusterHits(false),
fUsePFParticleHits(false),
fUsePFParticleSeeds(false),
fTKHAlg(pset.get<fhicl::ParameterSet>("Track3DKalmanHitAlg")),
fSpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
fHIncChisq(0),
fHPull(0),
fNumEvent(0)
{
   fHist = pset.get<bool>("Hist");
   fUseClusterHits = pset.get<bool>("UseClusterHits");
   fUsePFParticleHits = pset.get<bool>("UsePFParticleHits");
   fUsePFParticleSeeds = pset.get<bool>("UsePFParticleSeeds");
   fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
   fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
   fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
   if(fUseClusterHits && fUsePFParticleHits) {
      throw cet::exception("Track3DKalmanHit")
      << "Using input from both clustered and PFParticle hits.\n";
   }

   produces<std::vector<recob::Track> >();
   produces<std::vector<recob::SpacePoint>            >();
   produces<art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
   produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();
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
void trkf::Track3DKalmanHit::produce(art::Event & evt)
{
   ++fNumEvent;
   auto tracks = std::make_unique<std::vector<recob::Track>>();
   auto th_assn = std::make_unique<art::Assns<recob::Track, recob::Hit>>();
   auto thm_assn = std::make_unique< art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta > >();
   auto tsp_assn = std::make_unique<art::Assns<recob::Track, recob::SpacePoint>>();
   auto pfPartTrack_assns = std::make_unique<art::Assns<recob::PFParticle, recob::Track>>();
   auto spts = std::make_unique<std::vector<recob::SpacePoint>>();
   auto sph_assn = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();

   prepareForInput();
   auto inputs = getInput(evt);

   auto outputs = fTKHAlg.makeTracks(inputs);

   if(fHist) {
      fillHistograms(outputs);
   }

   createOutputs(evt, outputs, inputs, *tracks, *spts, *th_assn, *thm_assn, *tsp_assn, *sph_assn, *pfPartTrack_assns);

   fSpacePointAlg.clearHitMap();

   evt.put(std::move(tracks));
   evt.put(std::move(spts));
   evt.put(std::move(th_assn));
   evt.put(std::move(thm_assn));
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
}


//----------------------------------------------------------------------------
//A temporary method till we find an alternative
void trkf::Track3DKalmanHit::prepareForInput() {
   fSpacePointAlg.clearHitMap();
}

//----------------------------------------------------------------------------
// There are three modes of operation:
// 1.  Clustered hits (produces one hit collection).
// 2.  PFParticle hits (products one hit collection for each PFParticle).
// 3.  All hits (produces one hit collection).

trkf::KalmanInputs trkf::Track3DKalmanHit::getInput(const art::Event &evt) const{
   if (fUsePFParticleHits) return getPFParticleStuff(evt);
   return KalmanInputs(1, KalmanInput(fUseClusterHits ?
                                      getClusteredHits(evt):
                                      getAllHits(evt)));
}

//----------------------------------------------------------------------------
/// Fill a collection using clustered hits
trkf::Hits trkf::Track3DKalmanHit::getClusteredHits(const art::Event &evt) const{
   Hits hits;
   art::Handle< std::vector<recob::Cluster> > clusterh;
   evt.getByLabel(fClusterModuleLabel, clusterh);

   if (!clusterh.isValid()) return hits;
   // Get hits from all clusters.
   art::FindManyP<recob::Hit> hitsbycluster(clusterh, evt, fClusterModuleLabel);

   for(size_t i = 0; i < clusterh->size(); ++i) {
      std::vector< art::Ptr<recob::Hit> > clushits = hitsbycluster.at(i);
      hits.insert(hits.end(), clushits.begin(), clushits.end());
   }
   return hits;
}

//----------------------------------------------------------------------------
/// If both UseClusteredHits and UsePFParticles is false use this method to fill in hits
trkf::Hits trkf::Track3DKalmanHit::getAllHits(const art::Event &evt) const{
   Hits hits;
   art::Handle< std::vector<recob::Hit> > hith;
   evt.getByLabel(fHitModuleLabel, hith);
   if(!hith.isValid()) return hits;
   size_t nhits = hith->size();
   hits.reserve(nhits);

   for(size_t i = 0; i < nhits; ++i) {
      hits.push_back(art::Ptr<recob::Hit>(hith, i));
   }
   return hits;
}

//----------------------------------------------------------------------------
/// If UsePFParticles is true use this method to fill in hits
trkf::KalmanInputs trkf::Track3DKalmanHit::getPFParticleStuff(const art::Event &evt) const{
   KalmanInputs inputs;
   // Our program is to drive the track creation/fitting off the PFParticles in the data store
   // We'll use the hits associated to the PFParticles for each track - and only those hits.
   // Without a valid collection of PFParticles there is nothing to do here
   // We need a handle to the collection of clusters in the data store so we can
   // handle associations to hits.
   art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
   evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle);
   if (!pfParticleHandle.isValid()) return inputs;

   art::Handle<std::vector<recob::Cluster> > clusterHandle;
   evt.getByLabel(fClusterModuleLabel, clusterHandle);

   // If there are no clusters then something is really wrong
   if (!clusterHandle.isValid()) return inputs;

   // Recover the collection of associations between PFParticles and clusters, this will
   // be the mechanism by which we actually deal with clusters
   art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleModuleLabel);

   // Associations to seeds.
   art::FindManyP<recob::Seed> seedAssns(pfParticleHandle, evt, fPFParticleModuleLabel);

   // Likewise, recover the collection of associations to hits
   art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fClusterModuleLabel);

   inputs.reserve(pfParticleHandle->size());

   // While PFParticle describes a hierarchal structure, for now we simply loop over the collection
   for(size_t partIdx = 0; partIdx < pfParticleHandle->size(); partIdx++) {

      // Add a new empty hit collection.
      inputs.emplace_back();
      KalmanInput& kalman_input = inputs.back();
      kalman_input.pfPartPtr = art::Ptr<recob::PFParticle>(pfParticleHandle, partIdx);
      Hits& hits = kalman_input.hits;

      // Fill this hit vector by looping over associated clusters and finding the
      // hits associated to them
      std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(partIdx);

      for(auto const &cluster : clusterVec) {
         std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster.key());
         hits.insert(hits.end(), hitVec.begin(), hitVec.end());
      }

      // If requested, fill associated seeds.
      if(!fUsePFParticleSeeds) continue;
      art::PtrVector<recob::Seed>& seeds = kalman_input.seeds;
      std::vector<art::Ptr<recob::Seed> > seedVec = seedAssns.at(partIdx);
      seeds.insert(seeds.end(), seedVec.begin(), seedVec.end());
      art::FindManyP<recob::Hit> seedHitAssns(seedVec, evt, fPFParticleModuleLabel);

      for(size_t seedIdx = 0; seedIdx < seedVec.size(); ++seedIdx) {
         std::vector<art::Ptr<recob::Hit> > seedHitVec;
         //SS: why seedIdx can have an invalid value?
         try {
            seedHitVec = seedHitAssns.at(seedIdx);
         }
         catch(art::Exception const&) {
            seedHitVec.clear();
         }
         kalman_input.seedhits.emplace_back();
         Hits& seedhits = kalman_input.seedhits.back();
         seedhits.insert(seedhits.end(), seedHitVec.begin(), seedHitVec.end());
      }
   }
   return inputs;
}

//-------------------------------------------------------------------------------------
void trkf::Track3DKalmanHit::createOutputs(const art::Event &evt,
                                           std::vector<KalmanOutput> const &outputs,
                                           KalmanInputs const &inputs,
                                           std::vector<recob::Track> &tracks,
                                           std::vector<recob::SpacePoint> &spts,
                                           art::Assns<recob::Track, recob::Hit> &th_assn,
                                           art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> &thm_assn,
                                           art::Assns<recob::Track, recob::SpacePoint> &tsp_assn,
                                           art::Assns<recob::SpacePoint, recob::Hit> &sph_assn,
                                           art::Assns<recob::PFParticle, recob::Track> &pfPartTrack_assns)
{
   if(outputs.size()!= inputs.size()) return;

   size_t tracksSize(0);
   for(auto const &kalman_output : outputs) {
      tracksSize += kalman_output.tracks.size();
   }
   tracks.reserve(tracksSize);

   auto const tid = evt.getProductID<std::vector<recob::Track> >();
   auto const tidgetter = evt.productGetter(tid);

   auto const spacepointId = evt.getProductID<std::vector<recob::SpacePoint> >();
   auto const getter = evt.productGetter(spacepointId);
   for (size_t i = 0; i<outputs.size();++i) {
      // Recover the kalman tracks double ended queue
      const std::deque<KGTrack>& kalman_tracks = outputs[i].tracks;

      for(auto const& kalman_track:kalman_tracks) {

         // Add Track object to collection.
         recob::Track track;
         kalman_track.fillTrack(track, tracks.size());
         if(track.NumberTrajectoryPoints() < 2) {
            continue;
         }
         unsigned int numtrajpts = track.NumberTrajectoryPoints();
         tracks.emplace_back(std::move(track));
         // SS: tracks->size() does not change after this point in each iteration

         //fill hits from this track
         Hits trhits;
         std::vector<unsigned int> hittpindex; //hit-trajectory point index
         kalman_track.fillHits(trhits, hittpindex);
         if (hittpindex.back()>=numtrajpts){ //something is wrong
           throw cet::exception("Track3DKalmanHit")
             << "Last hit corresponds to trajectory point index "<<hittpindex.back()<<" while the number of trajectory points is "<<numtrajpts<<'\n';
         }

         // Make space points from this track.
         auto nspt = spts.size();
         fSpacePointAlg.fillSpacePoints(spts, kalman_track.TrackMap());

         std::vector<art::Ptr<recob::SpacePoint>> sptvec;
         for(auto ispt = nspt; ispt < spts.size(); ++ispt) {
            sptvec.emplace_back(spacepointId, ispt, getter);
            // Associate newly created space points with hits.
            // Make space point to hit associations.
            auto const &sphits = fSpacePointAlg.getAssociatedHits((spts)[ispt]);
            for(auto const& sphit: sphits) {
               sph_assn.addSingle(sptvec.back(), sphit);
            }
         }

         art::Ptr<recob::Track> aptr(tid, tracks.size()-1, tidgetter);

         // Make Track to Hit associations.
         for (size_t h = 0; h< trhits.size(); ++h){
            th_assn.addSingle(aptr, trhits[h]);
            recob::TrackHitMeta metadata(hittpindex[h], -1);
            thm_assn.addSingle(aptr, trhits[h], metadata);
         }

         // Make track to space point associations
         for (auto const& spt: sptvec) {
            tsp_assn.addSingle(aptr, spt);
         }

         // Optionally fill track-to-PFParticle associations.
         if (fUsePFParticleHits) {
            pfPartTrack_assns.addSingle(inputs[i].pfPartPtr, aptr);
         }
      } // end of loop over a given collection
   }
}

//----------------------------------------------------------------------------
/// Fill Histograms method
//fHPull and fHIncChisq are private data members of the class Track3DKalmanHit

void trkf::Track3DKalmanHit::fillHistograms(std::vector<KalmanOutput>& outputs)
{
   for(auto const &output : outputs) {
      const std::deque<KGTrack>& kalman_tracks = output.tracks;
      for (size_t i = 0; i < kalman_tracks.size(); ++i) {
         const KGTrack& trg = kalman_tracks[i];
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
