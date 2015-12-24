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

#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Seed.h"
#include "RecoAlg/KalmanFilterAlg.h"
#include "RecoAlg/SeedFinderAlgorithm.h"
#include "RecoObjects/KHitContainerWireLine.h"
#include "RecoObjects/KHitContainerWireX.h"
#include "RecoObjects/SurfXYZPlane.h"
#include "RecoObjects/PropAny.h"
#include "RecoObjects/KHit.h"
#include "Utilities/AssociationUtil.h"

#include "TH1F.h"

// Local functions.

namespace {
    
    // Local hit collection struct.
    struct HitCollection
    {
        art::Ptr<recob::PFParticle> pfPartPtr;
        art::PtrVector<recob::Hit> hits;
        art::PtrVector<recob::Seed> seeds;
        std::vector<art::PtrVector<recob::Hit> > seedhits;
        std::deque<trkf::KGTrack> tracks;
    };
    
    /* unused function
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
     */
    
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
        
        // Private methods.
        
        recob::Seed makeSeed(const art::PtrVector<recob::Hit>& hits) const;
        void FillHistograms(std::list<HitCollection>& hit_collections);
        void FilterHitsOnKalmanTracks(std::deque<KGTrack>& kalman_tracks, art::PtrVector<recob::Hit>& hits, art::PtrVector<recob::Hit>& seederhits, int ntracks);
        std::unique_ptr<KHitContainer> FillHitContainer(art::PtrVector<recob::Hit> &hits);
        
        
        // Fcl parameters.
        
        bool fHist;                         ///< Make histograms.
        bool fUseClusterHits;               ///< Use clustered hits as input.
        bool fUsePFParticleHits;            ///< Use PFParticle hits as input.
        bool fUsePFParticleSeeds;           ///< Use PFParticle seeds.
        std::string fHitModuleLabel;        ///< Unclustered Hits.
        std::string fClusterModuleLabel;    ///< Clustered Hits.
        std::string fPFParticleModuleLabel; ///< PFParticle label.
        bool fStoreNPPlane;                 ///< Store nonpreferred planes trajectory points.
        double fMaxTcut;                    ///< Maximum delta ray energy in MeV for restricted dE/dx.
        bool fDoDedx;                       ///< Global dE/dx enable flag.
        bool fSelfSeed;                     ///< Self seed flag.
        bool fLineSurface;                  ///< Line surface flag.
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
fUsePFParticleSeeds(false),
fStoreNPPlane(true),
fMaxTcut(0.),
fDoDedx(false),
fSelfSeed(false),
fLineSurface(false),
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
    fUsePFParticleSeeds = pset.get<bool>("UsePFParticleSeeds");
    fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
    fStoreNPPlane = pset.get<bool>("StoreNPPlane");
    fMaxTcut = pset.get<double>("MaxTcut");
    fDoDedx = pset.get<bool>("DoDedx");
    fSelfSeed = pset.get<bool>("SelfSeed");
    fLineSurface = pset.get<bool>("LineSurface");
    fMinSeedHits = pset.get<int>("MinSeedHits");
    fMinSeedChopHits = pset.get<int>("MinSeedChopHits");
    fMaxChopHits = pset.get<int>("MaxChopHits");
    fMaxSeedChiDF = pset.get<double>("MaxSeedChiDF");
    fMinSeedSlope = pset.get<double>("MinSeedSlope");
    fInitialMomentum = pset.get<double>("InitialMomentum");
    if(fProp != 0)
        delete fProp;
    fProp = new PropAny(fMaxTcut, fDoDedx);
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
    
    // Make collection of KGTracks where we will save our results, together with
    // corresponding hit collections.
    
    std::list<HitCollection> hit_collections;
    
    // Get Services.
    
    art::ServiceHandle<geo::Geometry> geom;
    
    // Reset space point algorithm.
    
    fSpacePointAlg.clearHitMap();
    
    // Get Hits.
    // There are three modes of operation:
    // 1.  Clustered hits (produces one hit collection).
    // 2.  PFParticle hits (products one hit collection for each PFParticle).
    // 3.  All hits (produces one hit collection).
    
    if(fUseClusterHits) {
        
        // Make one empty hit collection.
        
        hit_collections.emplace_back();
        art::PtrVector<recob::Hit>& hits = hit_collections.back().hits;
        
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
        if (pfParticleHandle.isValid())
        {
            // We need a handle to the collection of clusters in the data store so we can
            // handle associations to hits.
            art::Handle<std::vector<recob::Cluster> > clusterHandle;
            evt.getByLabel(fClusterModuleLabel, clusterHandle);
            
            // If there are no clusters then something is really wrong
            if (clusterHandle.isValid())
            {
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
                    
                    hit_collections.emplace_back();
                    HitCollection& hit_collection = hit_collections.back();
                    hit_collection.pfPartPtr = art::Ptr<recob::PFParticle>(pfParticleHandle, partIdx);
                    art::PtrVector<recob::Hit>& hits = hit_collection.hits;
                    
                    // Fill this hit vector by looping over associated clusters and finding the
                    // hits associated to them
                    std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(partIdx);
                    
                    for(const auto& cluster : clusterVec) {
                        std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster.key());
                        hits.insert(hits.end(), hitVec.begin(), hitVec.end());
                    }
                    
                    // If requested, fill associated seeds.
                    
                    if(fUsePFParticleSeeds) {
                        art::PtrVector<recob::Seed>& seeds = hit_collection.seeds;
                        std::vector<art::Ptr<recob::Seed> > seedVec = seedAssns.at(partIdx);
                        seeds.insert(seeds.end(), seedVec.begin(), seedVec.end());
                        art::FindManyP<recob::Hit> seedHitAssns(seedVec, evt, fPFParticleModuleLabel);
                        for(size_t seedIdx = 0; seedIdx < seedVec.size(); ++seedIdx) {
                            std::vector<art::Ptr<recob::Hit> > seedHitVec;
                            try {
                                seedHitVec = seedHitAssns.at(seedIdx);
                            }
                            catch(art::Exception x) {
                                seedHitVec.clear();
                            }
                            hit_collection.seedhits.emplace_back();
                            art::PtrVector<recob::Hit>& seedhits = hit_collection.seedhits.back();
                            seedhits.insert(seedhits.end(), seedHitVec.begin(), seedHitVec.end());
                        }
                    }
                }
            }
        }
    }
    else {
        
        // Make one empty hit collection.
        
        hit_collections.emplace_back();
        art::PtrVector<recob::Hit>& hits = hit_collections.back().hits;
        
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
    
    // Loop over hit collection / Kalman track combos.
    
    for(auto& hit_collection : hit_collections) {
        
        // Recover the kalman tracks double ended queue
        
        std::deque<KGTrack>& kalman_tracks = hit_collection.tracks;
        art::PtrVector<recob::Hit>& hits = hit_collection.hits;
        
        // The hit collection "hits" (just filled), initially containing all
        // hits, represents hits available for making tracks.  Now we will
        // fill a second hit collection called "seederhits", also initially
        // containing all hits, which will represent hits available for
        // making track seeds.  These collections are not necessarily the
        // same, since hits that are not suitable for seeds may still be
        // suitable for tracks.
        
        art::PtrVector<recob::Hit> seederhits = hits;
        
        // Start of loop.
        
        bool first = true;
        bool done = false;
        while(!done) {
            
            // Use remaining seederhits to make seeds.
            
            std::vector<art::PtrVector<recob::Hit> > hitsperseed;
            std::vector<recob::Seed> seeds;
            
            // On the first trip through this loop, try to use pfparticle-associated seeds.
            // Do this, provided the list of pfparticle-associated seeds and associated
            // hits are not empty.
            
            bool pfseed = false;
            if(first && hit_collection.seeds.size() > 0 && hit_collection.seedhits.size() > 0) {
                pfseed = true;
                seeds.reserve(hit_collection.seeds.size());
                for(const auto& pseed : hit_collection.seeds)
                    seeds.push_back(*pseed);
                hitsperseed.insert(hitsperseed.end(),
                                   hit_collection.seedhits.begin(), hit_collection.seedhits.end());
            }
            else {
                
                // On subsequent trips, or if there were no usable pfparticle-associated seeds,
                // attempt to generate our own seeds.
                
                if(seederhits.size()>0) {
                    if(fSelfSeed) {
                        
                        // Self seed - convert all hits into one big seed.
                        
                        seeds.emplace_back(makeSeed(seederhits));
                        hitsperseed.emplace_back();
                        hitsperseed.back().insert(hitsperseed.back().end(),
                                                  seederhits.begin(), seederhits.end());
                    }
                    else
                        seeds = fSeedFinderAlg.GetSeedsFromUnSortedHits(seederhits, hitsperseed);
                }
            }
            assert(seeds.size() == hitsperseed.size());
            first = false;
            
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
                    
                    // Don't chop pfparticle seeds or self seeds.
                    
                    int nchopmax = std::max(0, int((hpsit->size() - fMinSeedChopHits)/2));
                    if(pfseed || fSelfSeed)
                        nchopmax = 0;
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
                            mf::LogDebug("Track3DKalmanHit")
                            << "Seed found with " << seedhits.size() <<" hits.\n"
                            << "(x,y,z) = " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "\n"
                            << "(dx,dy,dz) = " << dir[0] << ", " << dir[1] << ", " << dir[2] << "\n";
                        } // if debug
                        
                        // Cut on the seed slope dx/ds.
                        
                        double dirlen = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
                        if(std::abs(dir[0]) >= fMinSeedSlope * dirlen) {
                            
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
                                std::unique_ptr<KHitContainer> pseedcont = FillHitContainer(seedhits);
                                
                                // Set the preferred plane to be the one with the most hits.
                                
                                unsigned int prefplane = pseedcont->getPreferredPlane();
                                fKFAlg.setPlane(prefplane);
                                if (mf::isDebugEnabled())
                                    mf::LogDebug("Track3DKalmanHit") << "Preferred plane = " << prefplane << "\n";
                                
                                // Build and smooth seed track.
                                
                                KGTrack trg0(prefplane);
                                bool ok = fKFAlg.buildTrack(trk, trg0, fProp, Propagator::FORWARD, *pseedcont,
                                                            fSelfSeed);
                                if(ok) {
                                    KGTrack trg1(prefplane);
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
                                        double mom0mag = std::sqrt(mom0[0]*mom0[0] + mom0[1]*mom0[1] + mom0[2]*mom0[2]);
                                        double mom1mag = std::sqrt(mom1[0]*mom1[0] + mom1[1]*mom1[1] + mom1[2]*mom1[2]);
                                        double dxds0 = mom0[0] / mom0mag;
                                        double dxds1 = mom1[0] / mom1mag;
                                        ok = ok && (std::abs(dxds0) > fMinSeedSlope &&
                                                    std::abs(dxds1) > fMinSeedSlope);
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
                                                
                                                std::unique_ptr<KHitContainer> ptrackcont = FillHitContainer(trackhits);
                                                
                                                
                                                // Extend the track.  It is not an error for the
                                                // extend operation to fail, meaning that no new hits
                                                // were added.
                                                
                                                bool extendok = fKFAlg.extendTrack(trg1, fProp, *ptrackcont);
                                                if(extendok)
                                                    nfail = 0;
                                                else
                                                    ++nfail;
                                                
                                                // Smooth the extended track, and make a new
                                                // unidirectionally fit track in the opposite
                                                // direction.
                                                
                                                KGTrack trg2(prefplane);
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
                                    mf::LogDebug("Track3DKalmanHit")
                                    << (ok? "Find track succeeded.": "Find track failed.") << "\n";
                                if(ok && !build_all)
                                    break;
                            } // for initial track
                            
                            // Loop over newly added tracks and remove hits contained on
                            // these tracks from hits available for making additional
                            // tracks or track seeds.
                            
                            FilterHitsOnKalmanTracks(kalman_tracks, hits, seederhits, ntracks);
                            
                        }
                    }
                }
            }
        }
    }
    
    // Fill histograms.
    
    if(fHist) {
        //replaced the for loops and body of code with a helper function
        FillHistograms(hit_collections);
    }
    
    // Process Kalman filter tracks into persistent objects.
    
    size_t tracksSize(0);
    for(const auto& hit_collection : hit_collections) {
        tracksSize += hit_collection.tracks.size();
    }
    tracks->reserve(tracksSize);
    
    for(auto& hit_collection : hit_collections) {
        
        // Recover the kalman tracks double ended queue
        const std::deque<KGTrack>& kalman_tracks = hit_collection.tracks;
        
        // Remember how many tracks are already converted
        size_t trackColStartIdx(tracks->size());
        
        for(std::deque<KGTrack>::const_iterator k = kalman_tracks.begin();
            k != kalman_tracks.end(); ++k) {
            const KGTrack& kalman_track = *k;
            
            // Add Track object to collection.
            
            recob::Track track;
            kalman_track.fillTrack(track, tracks->size(), fStoreNPPlane);
            if(track.NumberTrajectoryPoints() >= 2)
                tracks->emplace_back(std::move(track));
            
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
        
        // Optionally fill track-to-PFParticle associations.
        
        if (fUsePFParticleHits) {
            art::Ptr<recob::PFParticle>& pfPartPtr = hit_collection.pfPartPtr;
            
            // Create a vector of art ptrs to the tracks...
            std::vector<art::Ptr<recob::Track> > trackVec;
            
            size_t trackColEndIdx(tracks->size());
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

//----------------------------------------------------------------------------
/// Fill Histograms method
//fHPull and fHIncChisq are private data members of the class Track3DKalmanHit

void trkf::Track3DKalmanHit::FillHistograms(std::list<HitCollection>& hit_collections)
{
    for(const auto& hit_collection : hit_collections) {
        const std::deque<KGTrack>& kalman_tracks = hit_collection.tracks;
        
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
        fHIncChisq->Print("all");
        fHPull->Print("all");
    }
}

//----------------------------------------------------------------------------
/// Filter hits that are on kalman tracks.


void trkf::Track3DKalmanHit::FilterHitsOnKalmanTracks(std::deque<KGTrack>& kalman_tracks,
                                                      art::PtrVector<recob::Hit>& hits,
                                                      art::PtrVector<recob::Hit>& seederhits,
                                                      int ntracks){
    for(unsigned int itrk = ntracks; itrk < kalman_tracks.size(); ++itrk) {
        const KGTrack& trg = kalman_tracks[itrk];
        art::PtrVector<recob::Hit> track_used_hits;
        trg.fillHits(track_used_hits);
        FilterHits(hits, track_used_hits);
        FilterHits(seederhits, track_used_hits);
    }
}

//----------------------------------------------------------------------------
/// Fill hit container with either seedhits or filtered hits i.e. recob::Hit

std::unique_ptr<trkf::KHitContainer> trkf::Track3DKalmanHit::FillHitContainer(art::PtrVector<recob::Hit> &hits) {
    std::unique_ptr<KHitContainer> hitcont;
    if(fLineSurface) {
        KHitContainerWireLine* p = new KHitContainerWireLine;
        p->fill(hits, -1);
        hitcont.reset(p);
    }
    else {
        KHitContainerWireX* p = new KHitContainerWireX;
        p->fill(hits, -1);
        hitcont.reset(p);
    }
    return hitcont;
}



//----------------------------------------------------------------------------
/// Make seed method.
recob::Seed trkf::Track3DKalmanHit::makeSeed(const art::PtrVector<recob::Hit>& hits) const
{
    // Get Services.
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    // Do a linear 3D least squares for of y and z vs. x.
    // y = y0 + ay*(x-x0)
    // z = z0 + az*(x-x0)
    // Here x0 is the global average x coordinate of all hits in all planes.
    // Parameters y0, z0, ay, and az are determined by a simultaneous least squares
    // fit in all planes.
    
    // First, determine x0 by looping over every hit.
    
    double x0 = 0.;
    int n = 0;
    for(auto const& phit : hits) {
        const recob::Hit& hit = *phit;
        geo::WireID wire_id = hit.WireID();
        double time = hit.PeakTime();
        double x = detprop->ConvertTicksToX(time, wire_id);
        x0 += x;
        ++n;
    }
    
    // If there are no hits, return invalid seed.
    
    if(n == 0)
        return recob::Seed();
    
    // Find the average x.
    
    x0 /= n;
    
    // Now do the least squares fit proper.
    
    KSymMatrix<4>::type sm(4);
    KVector<4>::type sv(4);
    sm.clear();
    sv.clear();
    
    // Loop over hits (again).
    
    for(auto const& phit : hits) {
        const recob::Hit& hit = *phit;
        
        // Extract the angle, w and x coordinates from hit.
        
        geo::WireID wire_id = hit.WireID();
        const geo::WireGeo& wgeom = geom->Wire(wire_id);
        double xyz[3];
        wgeom.GetCenter(xyz);
        
        // Phi convention is the one documented in SurfYZPlane.h.
        
        double phi = TMath::PiOver2() - wgeom.ThetaZ();
        double sphi = std::sin(phi);
        double cphi = std::cos(phi);
        double w = -xyz[1]*sphi + xyz[2]*cphi;
        
        double time = hit.PeakTime();
        double x = detprop->ConvertTicksToX(time, wire_id);
        
        // Accumulate data for least squares fit.
        
        double dx = x-x0;
        
        sm(0, 0) += sphi*sphi;
        sm(1, 0) -= sphi*cphi;
        sm(1, 1) += cphi*cphi;
        sm(2, 0) += sphi*sphi * dx;
        sm(2, 1) -= sphi*cphi * dx;
        sm(2, 2) += sphi*sphi * dx*dx;
        sm(3, 0) -= sphi*cphi * dx;
        sm(3, 1) += cphi*cphi * dx;
        sm(3, 2) -= sphi*cphi * dx*dx;
        sm(3, 3) += cphi*cphi * dx*dx;
        
        sv(0) -= sphi * w;
        sv(1) += cphi * w;
        sv(2) -= sphi * w*dx;
        sv(3) += cphi * w*dx;
    }
    
    // Solve.
    
    bool ok = syminvert(sm);
    if(!ok)
        return recob::Seed();
    KVector<4>::type par(4);
    par = prod(sm, sv);
    
    double y0 = par(0);
    double z0 = par(1);
    double dydx = par(2);
    double dzdx = par(3);
    
    // Make seed.
    
    double dsdx = std::hypot(1., std::hypot(dydx, dzdx));
    
    double pos[3] = {x0, y0, z0};
    double dir[3] = {1./dsdx, dydx/dsdx, dzdx/dsdx};
    double poserr[3] = {0., 0., 0.};
    double direrr[3] = {0., 0., 0.};
    
    recob::Seed result(pos, dir, poserr, direrr);
    
    return result;
}
