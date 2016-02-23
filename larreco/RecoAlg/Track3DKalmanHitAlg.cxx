//////////////////////////////////////////////////////////////////////
///
/// \file   Track3DKalmanHitAlg.cxx
///
/// \brief  Track3DKalmanHitAlg.
///
/// \author
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/Track3DKalmanHitAlg.h"

// Local functions.

namespace {
   inline double calcMagnitude(const double *x){
      return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
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

/// Constructor.

trkf::Track3DKalmanHitAlg::Track3DKalmanHitAlg(const fhicl::ParameterSet& pset) :
fDoDedx(false),
fSelfSeed(false),
fMaxTcut(0.),
fLineSurface(false),
fMinSeedHits(0),
fMinSeedChopHits(0),
fMaxChopHits(0),
fMaxSeedChiDF(0.),
fMinSeedSlope(0.),
fInitialMomentum(0.),
fKFAlg(pset.get<fhicl::ParameterSet>("KalmanFilterAlg")),
fSeedFinderAlg(pset.get<fhicl::ParameterSet>("SeedFinderAlg")),
fProp(nullptr),
fNumTrack(0)
{
   mf::LogInfo("Track3DKalmanHitAlg") << "Track3DKalmanHitAlg instantiated.";
   
   // Load fcl parameters.
   
   reconfigure(pset);
}

/// Destructor.
trkf::Track3DKalmanHitAlg::~Track3DKalmanHitAlg()
{
   delete fProp;
}

/// Reconfigure method.
void trkf::Track3DKalmanHitAlg::reconfigure(const fhicl::ParameterSet& pset)
{
   fDoDedx = pset.get<bool>("DoDedx");
   fSelfSeed = pset.get<bool>("SelfSeed");
   fMaxTcut = pset.get<double>("MaxTcut");
   fLineSurface = pset.get<bool>("LineSurface");
   fMinSeedHits = pset.get<size_t>("MinSeedHits");
   fMinSeedChopHits = pset.get<int>("MinSeedChopHits");
   fMaxChopHits = pset.get<int>("MaxChopHits");
   fMaxSeedChiDF = pset.get<double>("MaxSeedChiDF");
   fMinSeedSlope = pset.get<double>("MinSeedSlope");
   fInitialMomentum = pset.get<double>("InitialMomentum");
   fKFAlg.reconfigure(pset.get<fhicl::ParameterSet>("KalmanFilterAlg"));
   fSeedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("SeedFinderAlg"));
   delete fProp;
   fProp = new PropAny(fMaxTcut, fDoDedx);
}


//-----------------
void trkf::Track3DKalmanHitAlg::generateKalmantracks(std::list<LocalKalmanStruct> &LocalKalmanStructList){
   // Loop over hit collection / Kalman track combos.
   for(auto& local_kalman_struct : LocalKalmanStructList) {
      
      // Recover the kalman tracks double ended queue
      std::cout << "Track3DKalmanHitAlg::generateKalmantracks\n";
      std::deque<KGTrack>& kalman_tracks = local_kalman_struct.tracks;
      art::PtrVector<recob::Hit>& hits = local_kalman_struct.hits;
      //std::cout << "Hits num: " << hits.size() << "\n";
      // The hit collection "hits" (just filled), initially containing all
      // hits, represents hits available for making tracks.  Now we will
      // fill a second hit collection called "unusedhits", also initially
      // containing all hits, which will represent hits available for
      // making track seeds.  These collections are not necessarily the
      // same, since hits that are not suitable for seeds may still be
      // suitable for tracks.
      
      art::PtrVector<recob::Hit> unusedhits = hits;
      
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
         auto const seedsize = local_kalman_struct.seeds.size();
         if(first && seedsize > 0 && local_kalman_struct.seedhits.size() > 0) {
            pfseed = true;
            seeds.reserve(seedsize);
            for(const auto& pseed : local_kalman_struct.seeds) {
               seeds.push_back(*pseed);
            }
            hitsperseed.insert(hitsperseed.end(),
                               local_kalman_struct.seedhits.begin(),
                               local_kalman_struct.seedhits.end());
         }
         else {
            // On subsequent trips, or if there were no usable pfparticle-associated seeds,
            // attempt to generate our own seeds.
            if(unusedhits.size()>0) {
               if(fSelfSeed) {
                  // Self seed - convert all hits into one big seed.
                  seeds.emplace_back(makeSeed(unusedhits));
                  hitsperseed.emplace_back();
                  hitsperseed.back().insert(hitsperseed.back().end(),
                                            unusedhits.begin(),
                                            unusedhits.end());
               }
               else
                  seeds = fSeedFinderAlg.GetSeedsFromUnSortedHits(unusedhits, hitsperseed);
            }
            //SS: what if seederhits.size(), that means we have no hits to process
            std::cout << "fSelfSeed: " << fSelfSeed << " unusedhits size " << unusedhits.size() << "\n";
            std::cout << "Seeds: " << seeds.size() << " , "<< hitsperseed.size() << "\n";
         }
         
         assert(seeds.size() == hitsperseed.size());
         first = false;
         
         
         if(seeds.size() == 0) { // Quit loop if we didn't find any new seeds.
            done = true;
            break;
         }
         else {
            processSeeds(pfseed, seeds, hitsperseed, unusedhits, hits, kalman_tracks);
         }
      }
   }
}


//----------------------------------------------------------------------------
/// method to return a seed to surface.

std::shared_ptr<trkf::Surface> trkf::Track3DKalmanHitAlg::makeSurface(const recob::Seed &seed,
                                                                   double *dir)
{
   double xyz[3];
   double err[3];   // Dummy.
   seed.GetPoint(xyz, err);
   seed.GetDirection(dir, err);
   if (mf::isDebugEnabled()) {
      mf::LogDebug("Track3DKalmanHit")
      //<< "Seed found with " << seedhits.size() <<" hits.\n"
      << "(x,y,z) = " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "\n"
      << "(dx,dy,dz) = " << dir[0] << ", " << dir[1] << ", " << dir[2] << "\n";
   } // if debug
   
   return std::shared_ptr<Surface>(new SurfXYZPlane(xyz[0], xyz[1], xyz[2],
                                                    dir[0], dir[1], dir[2]));
}



//----------------------------------------------------------------------------

bool trkf::Track3DKalmanHitAlg::processInitialtracks(const trkf::KTrack &trk,
                                                  art::PtrVector<recob::Hit>& seedhits,
                                                  art::PtrVector<recob::Hit>& hits,
                                                  std::deque<KGTrack>& kalman_tracks){
   
   // Fill hit container with current seed hits.
   std::unique_ptr<KHitContainer> pseedcont = fillHitContainer(seedhits);
   
   // Set the preferred plane to be the one with the most hits.
   unsigned int prefplane = pseedcont->getPreferredPlane();
   fKFAlg.setPlane(prefplane);
   if (mf::isDebugEnabled())
      mf::LogDebug("Track3DKalmanHit") << "Preferred plane = " << prefplane << "\n";
   
   // Build and smooth seed track.
   KGTrack trg0(prefplane);
   bool ok = fKFAlg.buildTrack(trk, trg0, fProp, Propagator::FORWARD, *pseedcont,
                               fSelfSeed);
   if(ok) ok = smoothandextendTrack(trg0, hits, prefplane, kalman_tracks);
   
   if (mf::isDebugEnabled())
      mf::LogDebug("Track3DKalmanHit")
      << (ok? "Find track succeeded.": "Find track failed.") << "\n";
   return ok;
}



//----------------------------------------------------------------------------


void trkf::Track3DKalmanHitAlg::processSeeds(bool pfseed,
                                          std::vector<recob::Seed>& seeds,
                                          std::vector<art::PtrVector<recob::Hit> >& hitsperseed,
                                          art::PtrVector<recob::Hit>& unusedhits,
                                          art::PtrVector<recob::Hit>& hits,
                                          std::deque<KGTrack>& kgtracks){
   std::vector<recob::Seed>::const_iterator sit = seeds.begin();
   std::vector<art::PtrVector<recob::Hit> >::const_iterator hpsit = hitsperseed.begin();
   //SS: use indexing
   //check for size of both containers
   //if (seeds.size() != hitsperseed.size())
   std::cout << "Track3DKalmanHitAlg::processSeeds\n";
   for (;sit != seeds.end() && hpsit != hitsperseed.end(); ++sit, ++hpsit) {
      //method to process Seeds
      const recob::Seed& seed = *sit;
      art::PtrVector<recob::Hit> unchoppedhits;
      
      // Chop a couple of hits off each end of the seed.
      chopHitsOffSeeds(*hpsit, pfseed, unchoppedhits);
      
      // Filter hits used by (chopped) seed from hits available to make future seeds.
      // No matter what, we will never use these hits for another seed.
      // This eliminates the possibility of an infinite loop.
      
      size_t initial_unusedhits = unusedhits.size();
      FilterHits(unusedhits, unchoppedhits);
      
      // Require that this seed be fully disjoint from existing tracks.
      //SS: replace this test with a method with appropriate name
      if(!(unchoppedhits.size() + unusedhits.size() == initial_unusedhits)) continue;
      
      // Convert seed into initial KTracks on surface located at seed point,
      // and normal to seed direction.
      double dir[3];
      std::shared_ptr<Surface> psurf = makeSurface(seed, dir);
      
      // Cut on the seed slope dx/ds.
      //SS: replace test name with a reasonable name
      if (!test1(dir)) continue;
      
      // Make one or two initial KTracks for forward and backward directions.
      //SS: it always makes 2 tracks
      std::vector<KTrack> iKtracks = makeInitialKtracks(psurf);
      
      // Loop over initial tracks.
      auto ntracks = kgtracks.size();   // Remember original track count.
      for(auto const &ktrk: iKtracks) {
         bool ok = processInitialtracks(ktrk, unchoppedhits, hits, kgtracks);
         if(ok && !fDoDedx) break;
      } // for initial track
      
      // Loop over newly added tracks and remove hits contained on
      // these tracks from hits available for making additional
      // tracks or track seeds.
      
      for(unsigned int itrk = ntracks; itrk < kgtracks.size(); ++itrk) {
         const KGTrack& trg = kgtracks[itrk];
         filterHitsOnKalmanTrack(trg, hits, unusedhits);
      }
      
   }
}
//----------------------------------------------------------------------------


std::vector<trkf::KTrack> trkf::Track3DKalmanHitAlg::makeInitialKtracks(const std::shared_ptr<trkf::Surface> psurf){
   // Assume muon (pdgid = 13).
   int pdg = 13; //SS: FIXME another constant?
   
   // SS: FIXME
   // It is a contant value inside a loop so I took it out
   // revisit the linear algebra stuff used here (use of ublas)
   // make a lambda here ... const TrackVector
   //
   TrackVector vec(5);
   vec(0) = 0.;
   vec(1) = 0.;
   vec(2) = 0.;
   vec(3) = 0.;
   vec(4) = (fInitialMomentum != 0. ? 1./fInitialMomentum : 2.);
   
   std::vector<KTrack> initial_tracks;
   
   // The build_all flag specifies whether we should attempt to make
   // tracks from all initial tracks, or alternatively, whether we
   // should declare victory and quit after getting a successful
   // track from one initial track.
   
   int ninit = 2;
   initial_tracks.reserve(ninit);
   initial_tracks.push_back(KTrack(psurf, vec, Surface::FORWARD, pdg));
   initial_tracks.push_back(KTrack(psurf, vec, Surface::BACKWARD, pdg));
   
   return initial_tracks;
}
//----------------------------------------------------------------------------

bool trkf::Track3DKalmanHitAlg::test1(const double *dir) {
   return std::abs(dir[0]) >= fMinSeedSlope * calcMagnitude(dir);
}


//----------------------------------------------------------------------------
/// Filter hits that are on kalman tracks.


void trkf::Track3DKalmanHitAlg::filterHitsOnKalmanTrack(const KGTrack& trg,
                                                     art::PtrVector<recob::Hit>& hits,
                                                     art::PtrVector<recob::Hit>& seederhits) const{
   art::PtrVector<recob::Hit> track_used_hits;
   trg.fillHits(track_used_hits);
   FilterHits(hits, track_used_hits);
   FilterHits(seederhits, track_used_hits);
}

//----------------------------------------------------------------------------
/// Fill hit container with either seedhits or filtered hits i.e. recob::Hit

std::unique_ptr<trkf::KHitContainer> trkf::Track3DKalmanHitAlg::fillHitContainer(const art::PtrVector<recob::Hit> &hits) const{
   // std::cout << "Track3DKalmanHit: fLineSurface " << fLineSurface << "\n";
   std::unique_ptr<KHitContainer> hitcont(fLineSurface ?
                                          static_cast<KHitContainer *>(new KHitContainerWireLine) :
                                          static_cast<KHitContainer *>(new KHitContainerWireX));
   hitcont->fill(hits, -1);
   return hitcont;
}




//----------------------------------------------------------------------------
/// Quality cuts on seed track.
// not sure if this function will be a candidate for generic interface

bool trkf::Track3DKalmanHitAlg::qualityCutsOnSeedTrack(const KGTrack &trg0) const{
   double mom0[3];
   double mom1[3];
   trg0.startTrack().getMomentum(mom0);
   trg0.endTrack().getMomentum(mom1);
   
   double dxds0 = mom0[0] / calcMagnitude(mom0);
   double dxds1 = mom1[0] / calcMagnitude(mom1);
   
   return (std::abs(dxds0) > fMinSeedSlope &&
           std::abs(dxds1) > fMinSeedSlope);
}

//----------------------------------------------------------------------------
/// Chop hits off of the end of seeds.
// Seems like seeds often end at delta rays, Michel electrons,
// or other pathologies.

// Don't chop pfparticle seeds or self seeds.

void trkf::Track3DKalmanHitAlg::chopHitsOffSeeds(art::PtrVector<recob::Hit>const& hpsit,
                                              bool pfseed,
                                              art::PtrVector<recob::Hit> &seedhits) const{
   //int nchopmax = std::max(0, int((hpsit.size() - fMinSeedChopHits)/2));
   //if(pfseed || fSelfSeed) nchopmax = 0;
   
   const int nchopmax = (pfseed || fSelfSeed) ?
   0:
   std::max(0, int((hpsit.size() - fMinSeedChopHits)/2));
   //SS: FIXME
   
   const int nchop = std::min(nchopmax, fMaxChopHits);
   art::PtrVector<recob::Hit>::const_iterator itb = hpsit.begin();
   art::PtrVector<recob::Hit>::const_iterator ite = hpsit.end();
   itb += nchop;
   ite -= nchop;
   seedhits.reserve(hpsit.size());
   for(art::PtrVector<recob::Hit>::const_iterator it = itb; it != ite; ++it)
      seedhits.push_back(*it);
}


bool trkf::Track3DKalmanHitAlg::smoothandextendTrack(KGTrack &trg0,
                                                  const art::PtrVector<recob::Hit> hits,
                                                  unsigned int prefplane,
                                                  std::deque<KGTrack>& kalman_tracks){
   KGTrack trg1(prefplane);
   bool ok = fKFAlg.smoothTrack(trg0, &trg1, fProp);
   if (!ok) return ok;
   
   // Now we have the seed track in the form of a KGTrack.
   // Do additional quality cuts.
   
   auto const n = trg1.numHits();
   auto const chisq = n * fMaxSeedChiDF;
   
   ok = n >= fMinSeedHits &&
   trg1.startTrack().getChisq() <= chisq &&
   trg1.endTrack().getChisq() <= chisq &&
   trg0.startTrack().getChisq() <= chisq &&
   trg0.endTrack().getChisq() <= chisq &&
   qualityCutsOnSeedTrack(trg0);
   
   if(!ok) return ok;
   
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
      std::unique_ptr<KHitContainer> ptrackcont = fillHitContainer(trackhits);
      
      // Extend the track.  It is not an error for the
      // extend operation to fail, meaning that no new hits
      // were added.
      if(fKFAlg.extendTrack(trg1, fProp, *ptrackcont)) nfail = 0;
      else ++nfail;
      
      // Smooth the extended track, and make a new
      // unidirectionally fit track in the opposite
      // direction.
      
      KGTrack trg2(prefplane);
      ok = fKFAlg.smoothTrack(trg1, &trg2, fProp);
      if(ok) {
         // Skip momentum estimate for constant-momentum tracks.
         if(fDoDedx) {
            fitnupdateMomentum(trg1, trg2);
         }
         trg1 = trg2;
      }
   }
   
   // Do a final smooth.
   
   if(!ok) return ok;
   
   ok = fKFAlg.smoothTrack(trg1, 0, fProp);
   if(!ok) return ok;
   
   // Skip momentum estimate for constant-momentum tracks.
   
   if(fDoDedx) {
      fitnupdateMomentum(trg1, trg1);
   }
   // Save this track.
   ++fNumTrack;
   std::cout << "In the alg: " << fNumTrack << "\n";
   kalman_tracks.push_back(trg1);
   return ok;
}


//----------------------------------------------------------------------------
/// fit and update method.
void trkf::Track3DKalmanHitAlg::fitnupdateMomentum(KGTrack& trg1, KGTrack& trg2) {
   KETrack tremom;
   if(fKFAlg.fitMomentum(trg1, fProp, tremom)) {
      fKFAlg.updateMomentum(tremom, fProp, trg2);
   }
}
//----------------------------------------------------------------------------
/// Make seed method.
recob::Seed trkf::Track3DKalmanHitAlg::makeSeed(const art::PtrVector<recob::Hit>& hits) const
{
   
   // Get Services.
   
   art::ServiceHandle<geo::Geometry> geom;
   const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   
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

