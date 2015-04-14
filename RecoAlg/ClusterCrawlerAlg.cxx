//////////////////////////////////////////////////////////////////////
///
/// ClusterCrawlerAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for crawling along a string of hits to make line clusters
/// Technical note in MicroBooNE docdb #2831
///
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm> // std::fill(), std::find(), std::sort()...

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Utilities/Exception.h" 

// LArSoft libraries
#include "Filters/ChannelFilter.h"
#include "SimpleTypesAndConstants/RawTypes.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/ClusterCrawlerAlg.h"


namespace cluster {
  //----------------------------------------------------------------------------
  // These static members need to be defined somewhere. Bah.
  constexpr ClusterCrawlerAlg::HitInCluster_t::ClusterID_t
    ClusterCrawlerAlg::HitInCluster_t::FreeHit;
  constexpr ClusterCrawlerAlg::HitInCluster_t::ClusterID_t
    ClusterCrawlerAlg::HitInCluster_t::NoHit;
  constexpr size_t ClusterCrawlerAlg::HitInCluster_t::InvalidHitIndex;
  
  void ClusterCrawlerAlg::HitInCluster_t::reset(size_t new_hits) {
    ClusterIDs.resize(new_hits);
    std::fill(ClusterIDs.begin(), ClusterIDs.end(), FreeHit);
  } // ClusterCrawlerAlg::HitInCluster_t::reset()
  
  
  
//------------------------------------------------------------------------------
  ClusterCrawlerAlg::ClusterCrawlerAlg(fhicl::ParameterSet const& pset)
  {
    reconfigure(pset);
  }

//------------------------------------------------------------------------------
  void ClusterCrawlerAlg::reconfigure(fhicl::ParameterSet const& pset)
  { 
    fNumPass            = pset.get<             unsigned short  >("NumPass");
    fMaxHitsFit         = pset.get< std::vector<unsigned short> >("MaxHitsFit");
    fMinHits            = pset.get< std::vector<unsigned short> >("MinHits");
    fNHitsAve           = pset.get< std::vector<unsigned short> >("NHitsAve");
    fChgCut             = pset.get< std::vector<float> >("ChgCut");
    fChiCut             = pset.get< std::vector<float> >("ChiCut");
    fMaxWirSkip         = pset.get< std::vector<unsigned short> >("MaxWirSkip");
    fMinWirAfterSkip    = pset.get< std::vector<unsigned short> >("MinWirAfterSkip");
    fKinkChiRat         = pset.get< std::vector<float> >("KinkChiRat");
    fKinkAngCut         = pset.get< std::vector<float> >("KinkAngCut");
    fDoMerge            = pset.get< std::vector<bool>  >("DoMerge");
    fTimeDelta          = pset.get< std::vector<float> >("TimeDelta");
    fMergeChgCut        = pset.get< std::vector<float> >("MergeChgCut");
    fFindVertices       = pset.get< std::vector<bool>  >("FindVertices");
    fLACrawl            = pset.get< std::vector<bool>  >("LACrawl");

    fMergeOverlap       = pset.get< bool   >("MergeOverlap");
    fChkClusterDS       = pset.get< bool   >("ChkClusterDS");
    fVtxClusterSplit    = pset.get< bool   >("VtxClusterSplit");
    fFindStarVertices   = pset.get< bool   >("FindStarVertices");
    fHitErrFac          = pset.get< float  >("HitErrFac");
    fLAClusAngleCut     = pset.get< float  >("LAClusAngleCut");
    fHitMergeChiCut     = pset.get< float  >("HitMergeChiCut");
    fAllowNoHitWire     = pset.get< unsigned short  >("AllowNoHitWire");
    fVertex3DCut        = pset.get< float  >("Vertex3DCut");
    
    fDebugPlane         = pset.get< short  >("DebugPlane");
    fDebugWire          = pset.get< short  >("DebugWire");
    fDebugHit           = pset.get< short  >("DebugHit");

    // some error checking
    bool badinput = false;
    if(fNumPass > fMaxHitsFit.size()) badinput = true;
    if(fNumPass > fMinHits.size()) badinput = true;
    if(fNumPass > fNHitsAve.size()) badinput = true;
    if(fNumPass > fChgCut.size()) badinput = true;
    if(fNumPass > fChiCut.size()) badinput = true;
    if(fNumPass > fMaxWirSkip.size()) badinput = true;
    if(fNumPass > fMinWirAfterSkip.size()) badinput = true;
    if(fNumPass > fKinkChiRat.size()) badinput = true;
    if(fNumPass > fKinkAngCut.size()) badinput = true;
    if(fNumPass > fDoMerge.size()) badinput = true;
    if(fNumPass > fTimeDelta.size()) badinput = true;
    if(fNumPass > fMergeChgCut.size()) badinput = true;
    if(fNumPass > fFindVertices.size()) badinput = true;
    if(fNumPass > fLACrawl.size()) badinput = true;

    if(badinput) throw art::Exception(art::errors::Configuration)
      << "ClusterCrawlerAlg: Bad input from fcl file";

  } // reconfigure
  
  
  // used for sorting hits on wires
  bool SortByLowHit(short i, short j) {return ((i > j));}

  // used for sorting clusters by length
  typedef std::pair<unsigned int, unsigned int> mypair;
  bool SortByLen(const mypair& L, const mypair& R) {return (L.first > R.first);}


  bool ClusterCrawlerAlg::SortByMultiplet
    (recob::Hit const& a, recob::Hit const& b)
  {
    // compare the wire IDs first:
    int cmp_res = a.WireID().cmp(b.WireID());
    if (cmp_res != 0) return cmp_res < 0; // order is decided, unless equal
    // decide by start time
    if (a.StartTick() != b.StartTick()) return a.StartTick() < b.StartTick();
    // if still undecided, resolve by local index
    return a.LocalIndex() < b.LocalIndex(); // if still unresolved, it's a bug!
  } // ClusterCrawlerAlg::SortByMultiplet()
  
  
//------------------------------------------------------------------------------
  void ClusterCrawlerAlg::ClearResults() {
    fHits.clear();
    tcl.clear();
    vtx.clear();
    vtx3.clear();
    fHitInCluster.clear();
  } // ClusterCrawlerAlg::ClearResults()
  
  
//------------------------------------------------------------------------------
  void ClusterCrawlerAlg::CrawlInit() {
    prt = false; vtxprt = false;
    NClusters = 0;  clBeginSlp = 0; clBeginSlpErr = 0; clBeginTim = 0;
    clBeginWir = 0; clBeginChg = 0; clEndSlp = 0;      clEndSlpErr = 0;
    clEndTim = 0;   clEndWir = 0;   clEndChg = 0;      clChisq = 0;
    clStopCode = 0; clProcCode = 0; fFirstWire = 0;
    fLastWire = 0; fAveChg = 0.; fChgSlp = 0.; pass = 0;
    fScaleF = 0; WireHitRange.clear();
    
    ClearResults();
  }


  void ClusterCrawlerAlg::RunCrawler(
    std::vector<recob::Hit> const& srchits
  )
  {
    // Run the ClusterCrawler algorithm - creating seed clusters and crawling upstream.

    CrawlInit();
    
    fHits = srchits; // plain copy of the sources; it's the base of our hit result
    
    // sort it as needed;
    // that is, sorted by wire ID number,
    // then by start of the region of interest in time, then by the multiplet
    std::sort(fHits.begin(), fHits.end(), &SortByMultiplet);
    
    fHitInCluster.reset(fHits.size()); // initialize all the hits as free
    
    // don't make clusters, just hits
    if(fNumPass == 0) return;
    
    if(fHits.size() < 3) return;
    if(fHits.size() > USHRT_MAX) { // not really, but let's keep things sane
      mf::LogWarning("CC")<<"Too many hits for ClusterCrawler "<<fHits.size();
      return;
    }
    
//    std::cout<<"CC: number of hits "<<fHits.size()<<"\n";
    
    for (geo::TPCID const& tpcid: geom->IterateTPCs()) {
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      for(plane = 0; plane < TPC.Nplanes(); ++plane){
        WireHitRange.clear();
        // define a code to ensure clusters are compared within the same plane
        clCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
        // fill the WireHitRange vector with first/last hit on each wire
        // dead wires and wires with no hits are flagged < 0
        GetHitRange(clCTP, WireHitRange, fFirstWire, fLastWire);

// sanity check
/*
  std::cout<<"Plane "<<plane<<" wire range "<<fFirstWire<<" "<<fLastWire;
  unsigned int nhts = 0;
  for(unsigned short wire = fFirstWire; wire < fLastWire; ++wire) {
    unsigned short index = wire - fFirstWire;
    unsigned short fhit = WireHitRange[index].first;
    unsigned short lhit = WireHitRange[index].second;
    for(unsigned short hit = fhit; hit < lhit; ++hit) {
      ++nhts;
      if(fHits[hit].WireID().Wire != wire) {
        std::cout<<"Bad wire "<<hit<<" "<<fHits[hit].WireID().Wire<<" "<<wire<<"\n";
        return;
      } // check wire
      if(fHits[hit].WireID().Plane != plane) {
        std::cout<<"Bad plane "<<hit<<" "<<fHits[hit].WireID().Plane<<" "<<plane<<"\n";
        return;
      } // check plane
      if(fHits[hit].WireID().TPC != tpc) {
        std::cout<<"Bad tpc "<<hit<<" "<<fHits[hit].WireID().TPC<<" "<<tpc<<"\n";
        return;
      } // check tpc
    } // hit
  } // wire
  std::cout<<" is OK. nhits "<<nhts<<"\n";;
*/
	if (WireHitRange.empty()||(fFirstWire == fLastWire)){
	  mf::LogWarning("ClusterCrawlerAlg")<<"No hits in "<<tpcid<<" plane "<<plane;
	  continue;
	}
	else {
	  LOG_DEBUG("ClusterCrawlerAlg")
	    << WireHitRange.size() << " hits in " << tpcid << " plane " << plane;
	}
	if (WireHitRange[0].first<0){
	  throw art::Exception(art::errors::LogicError)<<"ClusterCrawler WireHitRange[0].first = "<<WireHitRange[0].first;
	}
        fFirstHit = WireHitRange[0].first;
        raw::ChannelID_t channel = fHits[fFirstHit].Channel();
        // get the scale factor to convert dTick/dWire to dX/dU. This is used
        // to make the kink and merging cuts
        float wirePitch = geom->WirePitch(geom->View(channel));
        float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
        tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
        fScaleF = tickToDist / wirePitch;
        // convert Large Angle Cluster crawling cut to a slope cut
        if(fLAClusAngleCut > 0) 
          fLAClusSlopeCut = std::tan(3.142 * fLAClusAngleCut / 180.) / fScaleF;
        fMaxTime = detprop->NumberTimeSamples();
        fNumWires = geom->Nwires(plane);
        
        // look for clusters
        ClusterLoop();
      } // plane
      if(fVertex3DCut > 0) {
        // Match vertices in 3 planes
        VtxMatch(tpcid);
        Vtx3ClusterMatch(tpcid);
        // split clusters using 3D vertices
        Vtx3ClusterSplit(tpcid);
      }
      if(vtxprt) {
        mf::LogVerbatim("ClusterCrawlerAlg")<<"Clustering done in TPC ";
        PrintClusters();
      }
    } // for all tpcs

/*
  cstat = 0; tpc = 0;
  for(plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
    WireHitRange.clear();
    // define a code to ensure clusters are compared within the same plane
    clCTP = EncodeCTP(cstat, tpc, plane);
    // fill the WireHitRange vector with first/last hit on each wire
    // dead wires and wires with no hits are flagged < 0
    GetHitRange(clCTP, WireHitRange, fFirstWire, fLastWire);
    unsigned int nhts = 0, nhtsinCls = 0;
    for(unsigned short wire = fFirstWire; wire < fLastWire; ++wire) {
      unsigned short index = wire - fFirstWire;
      unsigned short fhit = WireHitRange[index].first;
      unsigned short lhit = WireHitRange[index].second;
      for(unsigned short hit = fhit; hit < lhit; ++hit) {
        if (isHitPresent(hit)) {
          ++nhts; // there is a hit
          if(!isHitFree(hit)) ++nhtsinCls;
        }
      } // hit
    } // wire
    std::cout<<"After CC: plane "<<plane<<" nhits "<<nhts
      <<" nhits in Clusters "<<nhtsinCls<<"\n";;
  } // plane
*/     
    
    // remove the hits that have become obsolete
    RemoveObsoleteHits();
    
    WireHitRange.clear(); 
    fcl2hits.clear();
    chifits.clear();
    hitnear.clear();
    
  } // RunCrawler
    
////////////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterLoop()
    {
      // looks for seed clusters in a plane and crawls along a trail of hits

      unsigned short nHitsUsed = 0;
      bool AllDone = false;
      for(unsigned short thispass = 0; thispass < fNumPass; ++thispass) {
        pass = thispass;
        // look for a starting cluster that spans a block of wires
        unsigned short span = 3;
        if(fMinHits[pass] < span) span = fMinHits[pass];
        for(unsigned short iwire = fLastWire; iwire > fFirstWire; iwire--) {
          unsigned short index = iwire - fFirstWire;
          // skip bad wires or no hits on the wire
          if(WireHitRange[index].first < 0) continue;
          unsigned short ifirsthit = WireHitRange[index].first;
          unsigned short ilasthit = WireHitRange[index].second;
          for(unsigned short ihit = ifirsthit; ihit < ilasthit; ++ihit) {
            bool ClusterAdded = false;
            recob::Hit const& hit = fHits[ihit];
            // skip used hits
            if(ihit > fHits.size()-1) {
              mf::LogError("ClusterCrawlerAlg")<<"ClusterLoop bad ihit "<<ihit;
              return;
            }
            // skip used and obsolete hits
            if (!isHitFree(ihit)) continue;
            // Check for a hit signal on the next DS wire
            bool SigOK = false;
            ChkSignal(iwire + 1, hit.PeakTime() - 10,
                      iwire + 1, hit.PeakTime() + 10, SigOK);
            // Don't start a seed cluster if there is a hit signal DS. 
            // This is an indicator that we might be trying
            // to start a cluster just US of shower blob.
            if(SigOK && hit.Multiplicity() > 1) continue;
            if((iwire - span + 1) < fFirstWire) continue;
            unsigned short jwire = iwire - span + 1;
            unsigned short jindx = jwire - fFirstWire;
            if(WireHitRange[jindx].first < 0) continue;
            // Find the hit on wire jwire that best matches a line between
            // a nearby vertex and hit ihit. No constraint if useHit < 0
            unsigned short useHit = 0;
            bool doConstrain = false;
            VtxConstraint(iwire, ihit, jwire, useHit, doConstrain);
            unsigned short jfirsthit = WireHitRange[jindx].first;
            unsigned short jlasthit = WireHitRange[jindx].second;
            for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
              recob::Hit const& other_hit = fHits[jhit];
              if(jhit > fHits.size()-1) {
                mf::LogError("ClusterCrawlerAlg")<<"ClusterLoop bad jhit "<<jhit;
                return;
              }
              // skip used and obsolete hits
              if (!isHitFree(jhit)) continue;
              // Vertex constraint
              if(doConstrain && jhit != useHit) continue;
              // start a cluster with these two hits
              fcl2hits.clear();
              chifits.clear();
              hitnear.clear();
              fAveChg = -1.;
              clEndChg = -1.;
              clStopCode = 0;
              clProcCode = pass;
              fcl2hits.push_back(ihit);
              chifits.push_back(0.);
              hitnear.push_back(0);
              fcl2hits.push_back(jhit);
              chifits.push_back(0.);
              hitnear.push_back(0);
              clpar[0] = other_hit.PeakTime();
              clpar[1] = (hit.PeakTime() - other_hit.PeakTime()) / (iwire - jwire);
              // check for hit width consistency. Large angle clusters should have
              // large rms hits or smaller rms hits if they part of a multiplet
              if(std::abs(clpar[1]) > 5) {
                if(pass == 0) continue;
                if(hit.RMS() < 4 && hit.Multiplicity() < 3) continue;
                if(other_hit.RMS() < 4 && other_hit.Multiplicity() < 3) continue;
              } // std::abs(clpar[1]) > 5
              clChisq = 0;
              // now look for hits to add on the intervening wires
              SigOK = false;
              bool HitOK = false;
              bool clok = true;
              for(unsigned short kwire = jwire+1; kwire < iwire; ++kwire) {
                // ensure this cluster doesn't cross a vertex
                if(CrawlVtxChk(kwire)) {
                  clok = false;
                  break;
                }
                AddHit(kwire, HitOK, SigOK);
                // no hit added and no nearby hit either
                if(!HitOK && !SigOK) {
                  clok = false;
                  break;
                }
              }
              // drop it?
              if(fcl2hits.size() < span || !clok) continue;
              // sort them by decreasing wire number
              // assume that this is the same as sorting by decreasing 
              // hit number. This only needs to be done on the starting cluster.
              // Hits will be added in the proper order by CrawlUS
              std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
              // do a real fit
              FitCluster();
              if(clChisq > 10.) continue;
              // check the charge ratio between the DS hit and the next-most
              // DS hit. This ratio is < 2 for a stopping particle. A large
              // ratio indicates that we are starting a cluster just US of a
              // high ionization region
              float chg0 = fHits[fcl2hits[0]].Integral();
              float chg1 = fHits[fcl2hits[1]].Integral();
              if(chg0 > 2 * chg1) continue;
              // check the hit RMS ratio
              chg0 = fHits[fcl2hits[0]].RMS();
              chg1 = fHits[fcl2hits[1]].RMS();
              float rmsrat = std::abs(chg0 - chg1) / (chg0 + chg1);
              if(rmsrat > 0.3) continue;
              // save the cluster begin info
              clBeginWir = iwire;
              clBeginTim = hit.PeakTime();
              clBeginSlp = clpar[1];
              // don't do a small angle crawl if the cluster slope is too large
              // and Large Angle crawling is NOT requested on this pass
              if(!fLACrawl[pass] && std::abs(clBeginSlp) > fLAClusSlopeCut) continue;
              // See if we are trying to start a cluster between a vertex
              // and a cluster that is associated to that vertex. If so, skip it
              if(CrawlVtxChk2()) continue;
              clBeginSlpErr = clparerr[1];
              clBeginChg = 0;
              // decide whether to crawl a large angle cluster. Requirements are:
              // 1) the user has set the LACluster angle cut > 0, AND
              // 2) the cluster slope exceeds the cut
              // Note that if condition 1 is met, normal cluster crawling is done
              // only if the slope is less than the cut
              if(fLACrawl[pass] && fLAClusSlopeCut > 0) {
                // LA cluster crawling requested
                if(std::abs(clBeginSlp) > fLAClusSlopeCut) {
                  LACrawlUS();
                } else {
                  CrawlUS();
                } // std::abs(clBeginSlp) > fLAClusAngleCut
              } else {
                // allow clusters of any angle
                CrawlUS();
              } // fLAClusSlopeCut > 0
              if(fcl2hits.size() >= fMinHits[pass]) {
                // it's long enough so save it
                clEndSlp = clpar[1]; // save the slope at the end
                clEndSlpErr = clparerr[1];
                TmpStore(); // store the cluster
                ClusterAdded = true;
                nHitsUsed += fcl2hits.size();
                AllDone = (nHitsUsed == fHits.size());
                break;
              }
              // kill it
            } // jhit
            if(ClusterAdded || AllDone) break;
          } // ihit
          if(AllDone) break;
        } // iwire

        // try to merge clusters 
        if(fDoMerge[pass]) ChkMerge();
        // form 2D vertices
        if(fFindVertices[pass]) FindVertices();

        if(AllDone) break;
        
      } // pass

      // Kill vertices that are close to long straight tracks. These
      // spurious vertices are likely due to delta rays on muon tracks
//      KillVertices();
      // Merge overlapping clusters
      if(fMergeOverlap) MergeOverlap();
      // Check the DS end of clusters
      if(fChkClusterDS) ChkClusterDS();
      // split clusters using vertices
      if(fVtxClusterSplit) VtxClusterSplit();
      // Look for 2D vertices with star topology - short, back-to-back clusters
      if(fFindStarVertices) FindStarVertices();

      if(fDebugPlane == (short)plane) {
        mf::LogVerbatim("ClusterCrawlerAlg")<<"Clustering done in plane "<<plane;
        PrintClusters();
      }
    
  } // ClusterLoop()

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeOverlap()
    {
      unsigned short icl, jcl, jj, tclsize;
      unsigned short imidwir, jmidwir;
      float ithb, jthe, da, dw, dt;
      
      prt = (fDebugWire == 666);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"MergeOverlap check. clCTP "<<clCTP;

      tclsize = tcl.size();
      for(icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].BeginVtx >= 0) continue;
        if(tcl[icl].tclhits.size() < 10) continue;
        ithb = std::atan(fScaleF * tcl[icl].BeginSlp);
        imidwir = 0.5 * (tcl[icl].BeginWir + tcl[icl].EndWir);
        for(jcl = 0; jcl < tclsize; ++jcl) {
          if(icl == jcl) continue;
          if(tcl[jcl].ID < 0) continue;
          if(tcl[jcl].CTP != clCTP) continue;
          if(tcl[jcl].tclhits.size() < 10) continue;
          jmidwir = 0.5 * (tcl[jcl].BeginWir + tcl[jcl].EndWir);
          // icl Begin is US of jcl End
          if(tcl[icl].BeginWir < tcl[jcl].EndWir) continue;
          // icl Begin is DS more than 50% of the length of jcl
          if(tcl[icl].BeginWir > jmidwir) continue;
          // jcl End is US more than 50% of the length of icl
          if(tcl[jcl].EndWir < imidwir) continue;
          // rough angle check
          jthe = std::atan(fScaleF * tcl[jcl].EndSlp);
          if(std::abs(ithb - jthe) > 0.2) continue;
          // fit hits on jcl near the Begin wire of icl to make
          // tighter cuts. Find the hits on jcl
          unsigned short jht = tcl[jcl].tclhits.size();
          for(jj = 0; jj < tcl[jcl].tclhits.size(); ++jj) {
            jht = tcl[jcl].tclhits[jj];
            if(fHits[jht].WireID().Wire <= tcl[icl].BeginWir) break;
          } // jj
          // GP I don't know how this could happen, but it would surely be bad:
          if (jht == tcl[jcl].tclhits.size()) continue;
          FitClusterMid(jcl, jht, 5);
          // compare the angle
          da = std::atan(fScaleF * clpar[1]) - ithb;
          dw = tcl[icl].BeginWir - fHits[jht].WireID().Wire;
          dt = clpar[0] + dw * clpar[1] - tcl[icl].BeginTim;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" overlap? "<<icl
    <<" and "<<jcl<<" da "<<da<<" dt "<<dt;
          if(std::abs(da) > 0.1) continue;
          if(std::abs(dt) > 10) continue;
          // merge clusters
          DoMerge(icl, jcl, 600);
        } // jcl
      } // icl
      
    } // MergeOverlap()

//////////////////////////////////////////
  void ClusterCrawlerAlg::RemoveObsoleteHits() {
    
    // check that no cluster hosts obsolete hits
    for (ClusterStore const& cluster: tcl) {
      for (unsigned short iHit: cluster.tclhits) {
        if (!isHitPresent(iHit)) {
          mf::LogError("ClusterCrawlerAlg") << "Cluster ID=" << cluster.ID
            << " claim to own hit #" << iHit
            << ", that does not exist any more! [this is called a ***BUG***]";
        }
      } // for cluster hit
    } // for cluster
    
    size_t iDestHit = 0, iSrcHit = 0;
    for (; iSrcHit < fHits.size(); ++iSrcHit) {
      if (!isHitPresent(iSrcHit)) continue; // this hit is going to disappear
      
      if (iDestHit != iSrcHit) { // need to move the hit
        // move the hit; this is easy
        fHits[iDestHit] = std::move(fHits[iSrcHit]);
        fHitInCluster.setCluster(iDestHit, fHitInCluster[iSrcHit]);
        // update the cluster reference
        if (!isHitFree(iSrcHit)) {
          size_t clusterIndex = size_t(fHitInCluster[iSrcHit] - 1);
          if (clusterIndex >= tcl.size()) {
            throw art::Exception(art::errors::LogicError)
              << "Hit #" << iSrcHit << " claims to belong to cluster ID="
              << clusterIndex
              << ", that does not exist! [this is called a ***BUG***]";
          }
          // tclhits a vector of indices of some type; we don't care which one
          auto& hits = tcl[clusterIndex].tclhits;
          auto iHitIndex = std::find(hits.begin(), hits.end(), iSrcHit);
          if (iHitIndex == hits.end()) {
            mf::LogError("ClusterCrawlerAlg")
              << "Hit #" << iSrcHit << " claims to belong to cluster ID="
              << clusterIndex
              << ", that does not admit it! [this is called a ***BUG***]";
          }
          else {
            *iHitIndex = iDestHit; // update the index
          }
        } // if hit is in a cluster
      } // if need to move the hit
      ++iDestHit;
    } // for iSrcHit
    
    // remove for good the invalid hits
    fHits.resize(iDestHit);
    fHitInCluster.resize(iDestHit);
    
    LOG_DEBUG("ClusterCrawlerAlg") << "RemoveObsoleteHits(): removed "
      << (iSrcHit - iDestHit) << "/" << iSrcHit << " hits";
    
  } // ClusterCrawlerAlg::RemoveObsoleteHits()


//////////////////////////////////////////
    void ClusterCrawlerAlg::ChkClusterDS() {
      // Try to extend clusters DS by a few wires. 
      // DS hits may not have been  included in a cluster if they have high 
      // multiplicity or high charge. 
      // Ref ClusterLoop cuts for starting a seed cluster.

  prt = (fDebugPlane == 3);

  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS clCTP "<<clCTP;

      float thhits, prevth, hitrms, rmsrat;
      bool ratOK;
      std::vector<unsigned short> dshits;
      const unsigned short tclsize = tcl.size();
      for(unsigned short icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore clusters that have a Begin vertex
        if(tcl[icl].BeginVtx >= 0) continue;
        // find the angle using the first 2 hits
        unsigned short ih0 = tcl[icl].tclhits[1];
        unsigned short ih1 = tcl[icl].tclhits[0];
        const float slp = (fHits[ih1].PeakTime()    - fHits[ih0].PeakTime()) / 
                          (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
        prevth = std::atan(fScaleF * slp);
        // move the "origin" to the first hit
        ih0 = ih1;
        unsigned short wire = fHits[ih0].WireID().Wire;
        hitrms = fHits[ih0].RMS();
        dshits.clear();
        // follow DS a few wires. Stop and do nothing if any encountered 
        // hit is associated with a cluster
        for(unsigned short ii = 0; ii < 2; ++ii) {
          ++wire;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS: Extend "
    <<tcl[icl].ID<<" to W:T "<<wire;
          if(wire > fLastWire) break;
          unsigned short index = wire - fFirstWire;
          // stop if no hits on this wire
          if(WireHitRange[index].first == -2) break;
        //  if(WireHitRange[index].first == -1) break; // FIXME: trouble; why??
          unsigned short firsthit = WireHitRange[index].first;
          unsigned short lasthit = WireHitRange[index].second;
          bool hitAdded = false;
          for(ih1 = firsthit; ih1 < lasthit; ++ih1) {
            if(!isHitFree(ih1)) continue;
            const float slp = (fHits[ih1].PeakTime() - fHits[ih0].PeakTime()) /
                              (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
            thhits = std::atan(fScaleF * slp);
            rmsrat = fHits[ih1].RMS() / hitrms;
            ratOK = rmsrat > 0.3 && rmsrat < 3;
            // require small angle and not wildly different width compared
            // to the first hit in the cluster
            if(ratOK && std::abs(thhits - prevth) < 0.3) {
              dshits.push_back(ih1);
              hitAdded = true;
              prevth = thhits;
              ih0 = ih1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Add hit "<<fHits[ih1].WireID().Wire
  <<":"<<(int)fHits[ih1].PeakTime()<<" rmsrat "<<rmsrat;
              break;
            }
          } // ih1
          // stop looking if no hit was added on this wire
          if(!hitAdded) break;
        } // ii
        // Found hits not associated with a different cluster
        if(dshits.size() > 0) {
          // put the tcl cluster into the working vectors
          TmpGet(icl);
          // clobber the hits
          fcl2hits.clear();
          // sort the DS hits
          if(dshits.size() > 1) std::sort(dshits.begin(), dshits.end(), SortByLowHit);
          // stuff them into fcl2hits
          fcl2hits = dshits;
          // Append the existing hits
          for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            // un-assign the hits so that TmpStore will re-assign them
            unsigned short iht = tcl[icl].tclhits[ii];
            fHitInCluster.setFree(iht);
            fcl2hits.push_back(iht);
          }
          clProcCode += 5000;
          pass = fNumPass - 1;
          FitClusterChg();
          clBeginChg = fAveChg;
          TmpStore();
          const size_t newcl = tcl.size() -1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Store "<<newcl;
          tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
          tcl[newcl].EndVtx = tcl[icl].EndVtx;
          // declare the old one obsolete
          tcl[icl].ID = -tcl[icl].ID;
        } // dshits.size() > 0
      } // icl
      
      // look for isolated stopping clusters. Check for hit multiplets
      // on the DS end and merge them
      float prTimeLo, prTimeHi;
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore clusters that have a Begin vertex
        if(tcl[icl].BeginVtx >= 0) continue;
        // ignore short clusters
        if(tcl[icl].tclhits.size() < 6) continue;
        const unsigned short nhts = 4;
        unsigned short nmult = 0;
        unsigned short loWire = 9999;
        unsigned short hiWire = 0;
        for(unsigned short ii = 0; ii < nhts; ++ii) {
          unsigned short iht = tcl[icl].tclhits[ii];
          // ignore hits close to a vertex
          if(tcl[icl].EndVtx >= 0) {
            unsigned short ivx = tcl[icl].EndVtx;
            if((fHits[iht].WireID().Wire - vtx[ivx].Wire) < 5) continue;
          } // tcl[icl].EndVtx >= 0
          if(fHits[iht].Multiplicity() > 1) ++nmult;
          if(fHits[iht].WireID().Wire < loWire) loWire = fHits[iht].WireID().Wire;
        //  if(fHits[iht].WireID().Wire > hiWire) hiWire = fHits[iht].WireID().Wire;
        } // ii
        if(nmult == 0) continue;
        // count the number of hits in this region that are not associated
        // with this cluster
        // inspect (at most) the 5 wires US of the cluster end and the 5 wires DS
        hiWire = std::min(tcl[icl].BeginWir + 5, (int) WireHitRange.size()); // or fLastWire+1?
        unsigned short nUsInClus = 0, nUsNotInClus = 0;
        unsigned short nDsInClus = 0, nDsNotInClus = 0;
        for(unsigned short wire = loWire; wire < hiWire; ++wire) {
          unsigned short index = wire - fFirstWire;
          // ignore dead wires and those with no hits
          if(WireHitRange[index].first < 0) continue;
          prTimeLo = tcl[icl].BeginTim + tcl[icl].BeginSlp * (wire - tcl[icl].BeginWir);
          prTimeHi = prTimeLo + 30;
          prTimeLo -= 30;
          for(unsigned short jht = WireHitRange[index].first; jht < WireHitRange[index].second; ++jht) {
            // ignore obsolete hits
            if(!isHitPresent(jht)) continue;
            if(fHits[jht].PeakTime() > prTimeHi) continue;
            if(fHits[jht].PeakTime() < prTimeLo) continue;
            // ignore hits associated with this cluster
            if(fHitInCluster[jht] == tcl[icl].ID) continue;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" not "
    <<fHits[jht].WireID().Wire<<":"<<(int)fHits[jht].PeakTime()
    <<" InClus "<<fHitInCluster[jht];
            if(wire <= tcl[icl].BeginWir) {
              if(isHitInCluster(jht)) {
                ++nUsInClus;
              } else {
                ++nUsNotInClus;
              }
            } else {
              if(isHitInCluster(jht)) {
                ++nDsInClus;
              } else {
                ++nDsNotInClus;
              }
            }
          } // jht
        } // wire
        // Skip if there are nearby US hits that are included in another cluster
        if(nUsInClus > 0) continue;
        // Also skip for DS nearby cluster hits
        if(nDsInClus > 0) continue;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS: icl "
    <<tcl[icl].ID<<" nmult "<<nmult
    <<" nUsNotInClus "<<nUsNotInClus
    <<" nDsNotInClus "<<nDsNotInClus;
        if(nmult > 0) {
          for(unsigned short ii = 0; ii < nhts; ++ii) {
            unsigned short iht = tcl[icl].tclhits[ii];
            // ignore hits close to a vertex
            if(tcl[icl].EndVtx >= 0) {
              unsigned short ivx = tcl[icl].EndVtx;
              if((fHits[iht].WireID().Wire - vtx[ivx].Wire) < 5) continue;
            } // tcl[icl].EndVtx >= 0
            if(fHits[iht].Multiplicity() > 1) MergeHits(iht);
          } // ii
        } // nmult > 0
        if(nDsNotInClus > 0) {
          // DS hits found that are not included. Try to attach
          unsigned short ih0 = tcl[icl].tclhits[1];
          const unsigned short ih1 = tcl[icl].tclhits[0];
          const float slp = (fHits[ih1].PeakTime()    - fHits[ih0].PeakTime()) /
                            (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
          prevth = std::atan(fScaleF * slp);
          // move the "origin" to the first hit
          ih0 = ih1;
          unsigned short wire = fHits[ih0].WireID().Wire;
          hitrms = fHits[ih0].RMS();
          dshits.clear();
          for(unsigned short ii = 0; ii < 2; ++ii) {
            ++wire;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS: Extend "
    <<tcl[icl].ID<<" to W:T "<<wire;
            if(wire > fLastWire) break;
            unsigned short index = wire - fFirstWire;
            if(WireHitRange[index].first == -2) break;
          //  if(WireHitRange[index].first == -1) break; // FIXME trouble, why?
            bool hitAdded = false;
            for(unsigned short ih1 = WireHitRange[index].first; ih1 < WireHitRange[index].second; ++ih1) {
              if(!isHitFree(ih1)) continue;
              const float slp = (fHits[ih1].PeakTime() - fHits[ih0].PeakTime()) /
                                (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
              thhits = std::atan(fScaleF * slp);
              if(std::abs(thhits - prevth) > 0.4) continue;
              rmsrat = fHits[ih1].RMS() / hitrms;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Found hit. Chk last hit RMS "
    <<hitrms<<" w W:T "<<fHits[ih1].WireID().Wire<<":"<<(int)fHits[ih1].PeakTime()
    <<" RMS "<<fHits[ih1].RMS()<<" rmsrat "<<rmsrat;
              if(rmsrat < 0.7 || rmsrat > 1.5) continue;
              dshits.push_back(ih1);
              hitAdded = true;
              prevth = thhits;
              ih0 = ih1;
              break;
            } // ih1
            if(!hitAdded) break;
          } // ii
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" dshits size "<<dshits.size();
          if(dshits.size() > 0) {
            // put the tcl cluster into the working vectors
            TmpGet(icl);
            // clobber the hits
            fcl2hits.clear();
            // sort the DS hits
            if(dshits.size() > 1) std::sort(dshits.begin(), dshits.end(), SortByLowHit); 
            fcl2hits = dshits;
            for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
              // un-assign the hits so that TmpStore will re-assign them
              unsigned short iht = tcl[icl].tclhits[ii];
              fHitInCluster.setFree(iht);
              fcl2hits.push_back(iht);
            } // ii
            clProcCode += 5000;
            pass = fNumPass - 1;
            FitClusterChg();
            clBeginChg = fAveChg;
            TmpStore();
            const size_t newcl = tcl.size() -1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Store "<<newcl;
            tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
            tcl[newcl].EndVtx = tcl[icl].EndVtx;
            // declare the old one obsolete
            tcl[icl].ID = -tcl[icl].ID;
          } // dshits.size() > 0
        } // nDsNotInClus > 0
      } // icl

    } // ChkClusterDS


//////////////////////////////////////////
  bool ClusterCrawlerAlg::CrawlVtxChk2()
  {
    // returns true if the (presumably short) cluster under construction
    // resides between a vertex and another cluster that is associated with
    // that vertex
    
    if(vtx.size() == 0) return false;
    if(fcl2hits.size() == 0) return false;
    
    unsigned short iht = fcl2hits.size() - 1;
    unsigned short icl;
    float wire0 = (0.5 + fHits[fcl2hits[iht]].WireID().Wire);
    float dt;
    float thclus = std::atan(fScaleF * clpar[1]);

    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      if(wire0 < vtx[iv].Wire) continue;
      if(wire0 > vtx[iv].Wire + 10) continue;
      dt = clpar[0] + (vtx[iv].Wire - wire0) * clpar[1] - vtx[iv].Time;
      if(std::abs(dt) > 10) continue;
      // cluster points to an US vertex. See if the angle is similar to
      // cluster associated with this vertex
      for(icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].EndVtx != iv) continue;
        if(std::abs(tcl[icl].EndAng - thclus) < fKinkAngCut[pass]) return true;
      }
    }
    
    return false;
    
  } // CrawlVtxChk2()

//////////////////////////////////////////
  bool ClusterCrawlerAlg::CrawlVtxChk(unsigned short kwire)
  {
    
    // returns true if the cluster is near a vertex on wire kwire
    if(vtx.size() == 0) return false;
    unsigned short iht = fcl2hits.size() - 1;
    float wire0 = (0.5 + fHits[fcl2hits[iht]].WireID().Wire);
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      if((unsigned short)(0.5 + vtx[iv].Wire) != kwire) continue;
      if(std::abs(prtime - vtx[iv].Time) < 10) return true;
    }
    return false;
  } // CrawlVtxChk()
//////////////////////////////////////////
    void ClusterCrawlerAlg::VtxConstraint(
        unsigned short iwire, unsigned short ihit, unsigned short jwire,
        unsigned short& useHit, bool& doConstrain)
    {
      // checks hits on wire jwire to see if one is on a line between a US vertex
      // and the hit ihit on wire iwire. If one is found, doConstrain is set true
      // and the hit index is returned.
      doConstrain = false;
      if(vtx.size() == 0) return;
      // no vertices made yet on the first pass
      if(pass == 0) return;
      // skip if vertices were not requested to be made on the previous pass
      if( !fFindVertices[pass - 1] ) return;
      
      unsigned short jindx = jwire - fFirstWire;
      unsigned short jfirsthit = WireHitRange[jindx].first;
      unsigned short jlasthit = WireHitRange[jindx].second;
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != clCTP) continue;
        // vertex must be US of the cluster
        if(vtx[iv].Wire > jwire) continue;
        // but not too far US
        if(vtx[iv].Wire < jwire - 10) continue;
        recob::Hit const& hit = fHits[ihit];
        clpar[0] = hit.PeakTime();
        clpar[1] = (vtx[iv].Time - hit.PeakTime()) / (vtx[iv].Wire - iwire);
        float prtime = clpar[0] + clpar[1] * (jwire - iwire);
        for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
          if(!isHitFree(jhit)) continue;
          const float tdiff = std::abs(fHits[jhit].TimeDistanceAsRMS(prtime));
          if(tdiff < 2.5) {
            useHit = jhit;
            doConstrain = true;
            return;
          }
        } // jhit
      } // iv
    } // VtxConstraint()


/////////////////////////////////////////
    void ClusterCrawlerAlg::VtxClusterSplit()
    {

      // split clusters that cross vertices

  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"VtxClusterSplit ";

      if(vtx.size() == 0) return;
      unsigned short tclsize = tcl.size();
      if(tclsize < 2) return;
      
      float dth;
      bool didit;

      for(unsigned short icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore short clusters
        if(tcl[icl].tclhits.size() < 5) continue;
        // check vertices
        didit = false;
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].CTP != clCTP) continue;
          // abandoned vertex?
          if(vtx[ivx].Wght <= 0) continue;
          // already assigned to this vertex?
          if(tcl[icl].BeginVtx == ivx) continue;
          if(tcl[icl].EndVtx == ivx) continue;
          // long dead-straight cluster?
          if(tcl[icl].tclhits.size() > 100) { 
            dth = tcl[icl].BeginAng - tcl[icl].EndAng;
            if(std::abs(dth) < 0.1) continue;
          }
          // vertex wire in-between the cluster ends?
          if(vtx[ivx].Wire < tcl[icl].EndWir) continue;
          if(vtx[ivx].Wire > tcl[icl].BeginWir) continue;
          // vertex time in-between the cluster ends?
          float hiTime = tcl[icl].BeginTim;
          if(tcl[icl].EndTim > hiTime) hiTime = tcl[icl].EndTim;
          if(vtx[ivx].Time > hiTime + 5) continue;
          float loTime = tcl[icl].BeginTim;
          if(tcl[icl].EndTim < loTime) loTime = tcl[icl].EndTim;
          if(vtx[ivx].Time < loTime - 5) continue;
          // find the hit on the cluster that is closest to the vertex on the
          // DS side
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Chk cluster "<<icl
    <<" with vertex "<<ivx;
          short ihvx = -99;
          // nSplit is the index of the hit in the cluster where we will
          // split it if all requirements are met
          unsigned short nSplit = 0;
          unsigned short nLop = 0;
          for(unsigned short ii = tcl[icl].tclhits.size()-1; ii > 0 ; --ii) {
            unsigned short iht = tcl[icl].tclhits[ii];
            ++nLop;
            if(fHits[iht].WireID().Wire >= vtx[ivx].Wire) {
              nSplit = ii;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Split cluster "<<tcl[icl].ID<<" at wire "<<fHits[iht].WireID().Wire
    <<" nSplit "<<nSplit;
              ihvx = iht;
              break;
            }
          } // ii
          // found the wire. Now make a rough time cut
          if(ihvx < 0) continue;
          short dtime = std::abs(short(fHits[ihvx].PeakTime() - vtx[ivx].Time));
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Delta time "<<dtime;
          if(dtime > 10) continue;
          // check the angle between the crossing cluster icl and the
          // clusters that comprise the vertex. 
          // First decide which end of cluster icl to use to define the angle
          float iclAng = 0.;
          if(nSplit > tcl[icl].tclhits.size() / 2) {
            iclAng = tcl[icl].EndAng;
          } else {
            iclAng = tcl[icl].BeginAng;
          }
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" iclAng "<<iclAng;
          // check angle wrt the the vertex clusters
          bool angOK = false;
          for(unsigned short jcl = 0; jcl < tclsize; ++jcl) {
            if(tcl[jcl].ID < 0) continue;
            if(tcl[jcl].CTP != clCTP) continue;
            if(tcl[jcl].BeginVtx == ivx) {
              if(std::abs(tcl[jcl].BeginAng - iclAng) > 0.4) {
                // large angle difference. Set the flag
                angOK = true;
                break;
              }
            } // tcl[jcl].BeginVtx == ivx
            if(tcl[jcl].EndVtx == ivx) {
              if(std::abs(tcl[jcl].EndAng - iclAng) > 0.4) {
                // large angle difference. Set the flag
                angOK = true;
                break;
              }
            } // tcl[jcl].EndVtx == ivx
          } // jcl
          // time to split or chop
          if(angOK) {
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Split/Chop at pos "<<nLop;
            if(nLop < 3) {
              // lop off hits at the US end
              // Put the cluster in the local arrays
              TmpGet(icl);
              for(unsigned short ii = 0; ii < nLop; ++ii) {
                unsigned short iht = fcl2hits[fcl2hits.size()-1];
                fHitInCluster.setFree(iht);
                fcl2hits.pop_back();
              }
              // store it
              clProcCode += 1000;
              TmpStore();
              unsigned short newcl = tcl.size() - 1;
              tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
              tcl[newcl].EndVtx = ivx;
              // declare this cluster obsolete
              tcl[icl].ID = -tcl[icl].ID;
            } else {
              // split the cluster into two
              // correct the split position
              ++nSplit;
              SplitCluster(icl, nSplit, ivx);
              tcl[icl].ProcCode += 1000;
              tcl[tcl.size()-1].ProcCode += 1000;
            }
            didit = true;
          } // angOK
          if(didit) break;
        } // ivx
      } // icl

    } // VtxClusterSplit()

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeHits(const unsigned short theHit) {
      // Merge all unused separate hits in the multiplet of which 
      // theHit is a member into one hit (= theHit).
      // Mark the merged hits other than theHit obsolete.
      // Hits in the multiplet that are associated with an existing cluster are
      // not affected.
      // Hit multiplicity is reworked (including all the hits in the multiplet).
      // Used hits have the multiplicity and index corrected too; the local
      // index reflects the peak time.
      // Note that theHit may or may not be marked free (usually, it is not)

      if(theHit > fHits.size() - 1) {
        mf::LogError("ClusterCrawlerAlg")<<"Bad theHit";
        return;
      }
      
      recob::Hit const& hit = fHits[theHit];

  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"MergeHits "
    <<hit.WireID().Wire<<":"<<(int)hit.PeakTime()
    <<" numHits "<<hit.Multiplicity();

      // number of hits in this hit multiplet
      std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(theHit);
      
      // ensure that this is a high multiplicity hit:
      if (MultipletRange.second <= MultipletRange.first) return;
      
      // calculate the Charge normalization factor using the hit information
      // instead of passing CCHitFinder ChgNorms all the way down here
      float chgNorm = 2.507 * hit.PeakAmplitude() * hit.RMS() / hit.Integral();

      short loTime = 9999;
      short hiTime = 0;
      unsigned short nGaus = 1;
      float hitSep;
      // number of hits that are close to theHit
      unsigned short nclose = 0;
      // find the time range for the hit multiplet
      for(size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
        if (!isHitPresent(jht)) continue; // just a ghost
        
        recob::Hit const& other_hit = fHits[jht];

  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<" P:W:T "<<plane<<":"<<other_hit.WireID().Wire<<":"<<(int)other_hit.PeakTime()
    <<" Amp "<<(int)other_hit.PeakAmplitude()
    <<" RMS "<<other_hit.RMS()
    <<" Charge "<<(int)other_hit.Integral()
    <<" InClus "<<fHitInCluster[jht];
  }
        // error checking
        if((other_hit.StartTick() != hit.StartTick())
          || (other_hit.WireID() != hit.WireID()))
        {
          mf::LogError("ClusterCrawlerAlg")<<"Hit multiplet ID error "
            << other_hit.WireID() << " @" << other_hit.StartTick()
              << " " << other_hit.LocalIndex() << "/" << other_hit.Multiplicity()
            << " vs. " << hit.WireID() << " @" << hit.StartTick()
              << " " << hit.LocalIndex() << "/" << hit.Multiplicity()
            ;
          return;
        }
        if (other_hit.Multiplicity() != hit.Multiplicity()) {
          mf::LogError("ClusterCrawlerAlg")
            << " hit #" << jht << " in the same multiplet as #" << theHit
            << " has different multiplicity!"
            << "\n hit #" << theHit << ": " << hit
            << "\n hit #" << jht << ": " << other_hit;
        }
        // hit is not used by another cluster
        if(!isHitFree(jht)) continue;
        short arg = (short)(other_hit.PeakTimeMinusRMS(3));
        if(arg < loTime) loTime = arg;
        arg = (short)(other_hit.PeakTimePlusRMS(3));
        if(arg > hiTime) hiTime = arg;
        if(jht != theHit) ++nGaus;
        hitSep = std::abs(other_hit.PeakTime() - hit.PeakTime()) / other_hit.RMS();
        if(jht != theHit && hitSep < 3) ++nclose;
      } // jht
      // all hits in the multiplet other than this one used?
      if(nGaus < 2) return;
      
      // the hits in this multiplet will have this multiplicity from now on
      const short int NewMultiplicity = hit.Multiplicity() + 1 - nGaus;
      
      // no reasonably close hits?
//      if(nhitnear < 4 && nclose == 0) return;
      if(loTime < 0) loTime = 0;
      ++hiTime;
      // define a signal shape, fill it with zeros
      std::vector<double> signal(hiTime - loTime, 0.);
      // now add the Gaussians for each hit
      double chgsum = 0.;
      for(size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
        recob::Hit const& other_hit = fHits[jht];
        if(jht != theHit) {
          // hit used in another cluster
          if(!isHitFree(jht)) continue;
          // declare this hit obsolete
          fHitInCluster.makeObsolete(jht);
        } // jht != theHit
        // add up the charge
        chgsum += other_hit.Integral();
        for(unsigned short time = loTime; time < hiTime; ++time) {
          unsigned short indx = time - loTime;
          double arg = (other_hit.PeakTime() - (double)time) / other_hit.RMS();
          signal[indx] += other_hit.PeakAmplitude() * exp(-0.5 * arg * arg);
        } // time
      } // jj
      // find the average weighted time
      double sigsum = 0.;
      double sigsumt = 0.;
      for(unsigned short time = loTime; time < hiTime; ++time) {
        sigsum  += signal[time - loTime];
        sigsumt += signal[time - loTime] * time;
      }
      if(sigsum == 0.) {
        mf::LogError("ClusterCrawlerAlg")<<"MergeHits: bad sum";
        return;
      }
      double aveTime = sigsumt / sigsum;
      // find the RMS
      sigsumt = 0.;
      for(unsigned short time = loTime; time < hiTime; ++time) {
        double dtime = time - aveTime;
        sigsumt += signal[time - loTime] * dtime * dtime;
      }
      const float RMS = std::sqrt(sigsumt / sigsum);
      // find the amplitude from the integrated charge and the RMS
      const float amplitude = chgsum * chgNorm/ (2.507 * RMS);
      
      // modify the hit "in place" (actually completely overwrite it...)
      // TODO a lot of these quantities need revamp!!
      fHits[theHit] = recob::Hit(
        hit.Channel(),
        hit.StartTick(),
        hit.EndTick(),
        aveTime,                // peak_time
        hit.SigmaPeakTime(),
        RMS,                    // rms
        amplitude,              // peak_amplitude
        hit.SigmaPeakAmplitude(),
        hit.SummedADC(),
        chgsum,                 // hit_integral
        hit.SigmaIntegral(),
        NewMultiplicity,        // multiplicity
        0,                      // local index
        hit.GoodnessOfFit(),
        hit.DegreesOfFreedom(),
        hit.View(),
        hit.SignalType(),
        hit.WireID()
        );
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<" theHit "<<fHits[theHit].WireID().Wire<<":"<<(int)aveTime
    <<" RMS "<<std::setprecision(1)<<fHits[theHit].RMS()
    <<" chg "<<(int)chgsum<<" Amp "<<(int)fHits[theHit].PeakAmplitude();
  }
      
      FixMultipletLocalIndices(MultipletRange.first, MultipletRange.second);
      
    } // MergeHits()

/////////////////////////////////////////
    void ClusterCrawlerAlg::FixMultipletLocalIndices
      (size_t begin, size_t end, short int multiplicity /* = -1 */)
    {
      //
      // Resets multiplicity and local index of the hits in the range.
      // All hits are assumed to be in the same multiplet.
      // All hits that are not obsolete are given a multiplicity equal to the
      // number of non-obsolete hits in the multiplet, and the local index is
      // assigned as an increasing number starting from 0 with the first
      // non-obsolete hit on.
      //
      
      // first pass: determine the actual number of hits in the multiplet
      if (multiplicity < 0) {
        multiplicity = 0;
        for (size_t iHit = begin; iHit < end; ++iHit) {
          if (!isHitPresent(iHit)) continue;
          ++multiplicity;
        } // for
      } // if no valid multiplicity is given
      
      // second pass: assign the correct multiplicity
      short int local_index = 0;
      for (size_t iHit = begin; iHit < end; ++iHit) {
        if (!isHitPresent(iHit)) continue;
        
        // copy everything but overwrite the local index and multiplicity
        // TODO use a write wrapper!
        recob::Hit const& hit = fHits[iHit];
        fHits[iHit] = recob::Hit(
          hit.Channel(),
          hit.StartTick(),
          hit.EndTick(),
          hit.PeakTime(),
          hit.SigmaPeakTime(),
          hit.RMS(),
          hit.PeakAmplitude(),
          hit.SigmaPeakAmplitude(),
          hit.SummedADC(),
          hit.Integral(),
          hit.SigmaIntegral(),
          multiplicity,        // multiplicity
          local_index,         // local index
          hit.GoodnessOfFit(),
          hit.DegreesOfFreedom(),
          hit.View(),
          hit.SignalType(),
          hit.WireID()
          );
        
        ++local_index;
      } // for
      
    } // FixMultipletLocalIndices()

/////////////////////////////////////////
    void ClusterCrawlerAlg::FindStarVertices()
    {
      // Make 2D vertices with a star topology with short back-to-back
      // clusters. Vertices must reside on the US end of the longest
      // cluster, so vertex finding uses the End information only.
      // This routine is called after all passes are completed
      // in the current CTP
      if(tcl.size() < 2) return;

  vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")<<"FindStarVertices";
    PrintClusters();
  }

      unsigned short vtxSizeIn = vtx.size();

      float fvw = 0.;
      float fvt = 0.;
      float dsl = 0, dth = 0;
      float es1 = 0, es2 = 0;
      float eth1 = 0, eth2 = 0;
      float bt1 = 0, bt2 = 0;
      float et1 = 0, et2 = 0;
      float bw1 = 0, bw2 = 0;
      float ew1 = 0, ew2 = 0;
      float lotime = 0, hitime = 0, nwirecut = 0;
      unsigned short tclsize = tcl.size();
      for(unsigned short it1 = 0; it1 < tclsize - 1; ++it1) {
        // ignore obsolete clusters
        if(tcl[it1].ID < 0) continue;
        if(tcl[it1].CTP != clCTP) continue;
        // long dead-straight cluster?
        if(tcl[it1].tclhits.size() > 100) { 
          dth = tcl[it1].BeginAng - tcl[it1].EndAng;
          if(std::abs(dth) < 0.1) continue;
        }
        es1 = tcl[it1].EndSlp;
        eth1 = tcl[it1].EndAng;
        ew1 = tcl[it1].EndWir;
        et1 = tcl[it1].EndTim;
        bw1 = tcl[it1].BeginWir;
        bt1 = tcl[it1].BeginTim;
        for(unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
          if(tcl[it2].ID < 0) continue;
          if(tcl[it2].CTP != clCTP) continue;
          // long dead-straight cluster?
          if(tcl[it2].tclhits.size() > 100) { 
            dth = tcl[it2].BeginAng - tcl[it2].EndAng;
            if(std::abs(dth) < 0.05) continue;
          }
          es2 = tcl[it2].EndSlp;
          eth2 = tcl[it2].EndAng;
          ew2 = tcl[it2].EndWir;
          et2 = tcl[it2].EndTim;
          bw2 = tcl[it2].BeginWir;
          bt2 = tcl[it2].BeginTim;
          // require angle difference
          dth = std::abs(eth1 - eth2);
          if(dth < 0.3) continue;
          dsl = es2 - es1;
          fvw = (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
          // intersection within the wire boundaries
          if(fvw < ew1 || fvw > bw1) continue;
          if(fvw < ew2 || fvw > bw2) continue;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo5 vtx wire "<<fvw;
          // ensure that the intersection is close to the US end of the longer
          // cluster if it is more than 20 hits long
          if(tcl[it1].tclhits.size() > tcl[it2].tclhits.size()) {
            if(tcl[it1].tclhits.size() > 20) {
              nwirecut = 0.5 * tcl[it1].tclhits.size();
              if((fvw - ew1) > nwirecut) continue;
            }
          } else {
            if(tcl[it2].tclhits.size() > 20) {
              nwirecut = 0.5 * tcl[it2].tclhits.size();
              if((fvw - ew2) > nwirecut) continue;
            }
          }
          fvt = et1 + (fvw - ew1) * es1;
          // and time boundaries
          lotime = 9999;
          if(et1 < lotime) lotime = et1;
          if(bt1 < lotime) lotime = bt1;
          if(et2 < lotime) lotime = et2;
          if(bt2 < lotime) lotime = bt2;
          hitime = 0;
          if(et1 > hitime) hitime = et1;
          if(bt1 > hitime) hitime = bt1;
          if(et2 > hitime) hitime = et2;
          if(bt2 > hitime) hitime = bt2;
          if(fvt < lotime || fvt > hitime) continue;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<" vertex time "<<fvt
    <<" lotime "<<lotime<<" hitime "<<hitime;
          unsigned short vw = (0.5 + fvw);
          // ensure that the vertex is near a hit on both clusters
          unsigned short pos1 = 0;
          for(unsigned short ii = 0; ii < tcl[it1].tclhits.size(); ++ii) {
            if(fHits[tcl[it1].tclhits[ii]].WireID().Wire <= vw) {
              if(std::abs(int(fHits[tcl[it1].tclhits[ii]].PeakTime() - fvt)) < 10) pos1 = ii;
              break;
            }
          } // ii
          // vertex is not near a hit on cluster 1
          if(pos1 == 0) continue;
          unsigned short pos2 = 0;
          for(unsigned short ii = 0; ii < tcl[it2].tclhits.size(); ++ii) {
            if(fHits[tcl[it2].tclhits[ii]].WireID().Wire <= vw) {
              if(std::abs(int(fHits[tcl[it2].tclhits[ii]].PeakTime() - fvt)) < 10) pos2 = ii;
              break;
            }
          } // ii
          // vertex is not near a hit on cluster 2
          if(pos2 == 0) continue;
          // ensure we aren't near an existing vertex
          for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
            if(std::abs(fvw - vtx[iv].Wire) < 2 &&
               std::abs(fvt - vtx[iv].Time) < 10) continue;
          }
          // make a new vertex
          VtxStore newvx;
          newvx.Wire = fvw;
          newvx.Time = fvt;
          newvx.Wght = 1;
          newvx.Topo = 5;
          newvx.CTP = clCTP;
          vtx.push_back(newvx);
          unsigned short ivx = vtx.size() - 1;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<" new vertex "<<ivx
    <<" cluster "<<tcl[it1].ID<<" split pos "<<pos1;
          SplitCluster(it1, pos1, ivx);
          tcl[it1].ProcCode += 1000;
          tcl[tcl.size()-1].ProcCode += 1000;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<" new vertex "<<ivx
    <<" cluster "<<tcl[it2].ID<<" split pos "<<pos2;
          SplitCluster(it2, pos2, ivx);
          tcl[it2].ProcCode += 1000;
          tcl[tcl.size()-1].ProcCode += 1000;
          break;
        } // it2
      } // it1
            
      if(vtx.size() > vtxSizeIn) {
        // try to split other clusters
        VtxClusterSplit();
        // try to attach other clusters to it
        VertexCluster(vtx.size() - 1);
        VtxWghtAndFit(clCTP);
      } // new vertex

  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")<<"Vertices "<<vtx.size();
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      mf::LogVerbatim("ClusterCrawlerAlg")
        <<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time
        <<" wght "<<(int)vtx[iv].Wght<<" topo "<<vtx[iv].Topo;
    }
    PrintClusters();
  }
      
    } // FindStarVertices()

//////////////////////////////////////////
    void ClusterCrawlerAlg::VertexCluster(unsigned short iv)
    {
      // try to attach clusters to the specified vertex
      if(vtx[iv].Wght < 0) return;
      short dwib, dwie;
      float dw, dt;
        
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != vtx[iv].CTP) continue;
        dwib = std::abs(vtx[iv].Wire - tcl[it].BeginWir);
        dwie = std::abs(vtx[iv].Wire - tcl[it].EndWir);
        if(dwib < dwie && dwib < 10) {
          // Cluster Begin is closer
          if(vtx[iv].Wire > tcl[it].BeginWir + 1) continue;
          dw = vtx[iv].Wire - tcl[it].BeginWir;
          dt = tcl[it].BeginTim + tcl[it].BeginSlp * dw - vtx[iv].Time;
          if(std::abs(dt) > 10.) continue;
          tcl[it].BeginVtx = iv;
        } else if(dwie < 10){
          // Cluster End is closer
          if(vtx[iv].Wire < tcl[it].EndWir - 1) continue;
          dw = vtx[iv].Wire - tcl[it].EndWir;
          dt = tcl[it].EndTim + tcl[it].EndSlp * dw - vtx[iv].Time;
          if(std::abs(dt) > 10.) continue;
          tcl[it].EndVtx = iv;
        }
      } // icl
    } // VertexCluster

//////////////////////////////////////////
    void ClusterCrawlerAlg::VtxWghtAndFit(CTP_t inCTP)
      {

      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != inCTP) continue;
        if(vtx[iv].Wght < 0) continue;
        // fit the vertices
        float ChiDOF = 0.;
        FitVtx(iv, ChiDOF);
      }

      // set the vertex weights
      float cw = 0;
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != inCTP) continue;
        // cluster weight = number of hits, saturated at 10
        cw = tcl[it].tclhits.size();
        if(cw > 10) cw = 10;
        if(tcl[it].BeginVtx >=0) vtx[tcl[it].BeginVtx].Wght += cw;
        if(tcl[it].EndVtx >=0) vtx[tcl[it].EndVtx].Wght += cw;
      }

    } // VtxWghtAndFit
/*
/////////////////////////////////////////
    void ClusterCrawlerAlg::KillVertices()
    {
      // kill vertices that are close to long straight clusters. These
      // are most likely due to numerous short delta ray clusters
      if(vtx.size() == 0) return;
      
      unsigned short iht, vWire;
      bool killit;
      
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != clCTP) continue;
        if(tcl[it].tclhits.size() < 100) continue;
        if(std::abs(tcl[it].BeginAng - tcl[it].EndAng) > 0.1) continue;
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].CTP != clCTP) continue;
          if(vtx[ivx].Wght < 0) continue;
          // don't touch vertices near the long cluster ends
          if(vtx[ivx].Wire < tcl[it].EndWir + 10) continue;
          if(vtx[ivx].Wire > tcl[it].BeginWir - 10) continue;
          vWire = (0.5 + vtx[ivx].Wire);
          killit = false;
          for(unsigned short ii = 0; ii < tcl[it].tclhits.size(); ++ii) {
            iht = tcl[it].tclhits[ii];
            if(fHits[iht].WireID().Wire <= vWire) {
              if(std::abs(fHits[iht].PeakTime() - vtx[ivx].Time) < 60) {
                killit = true;
                break;
              }
              break;
            } // fHits[iht].WireID().Wire <= vWire
          } // ii
          if(killit) {
            vtx[ivx].Wght = -1;
            for(unsigned short jt = 0; jt < tcl.size(); ++jt) {
              if(tcl[jt].BeginVtx == ivx) tcl[jt].BeginVtx = -99;
              if(tcl[jt].EndVtx == ivx) tcl[jt].EndVtx = -99;
            } // jt
          } // killit
        } // ivx
      } // it
      
    } // KillVertices
*/
/////////////////////////////////////////
    void ClusterCrawlerAlg::FindVertices()
    {
      // try to make 2D vertices
      
      if(tcl.size() < 2) return;

      // form vertices starting with the longest
      std::map<unsigned short, unsigned short> sortindex;
      SortByLength(tcl, clCTP, sortindex);
      
      float nwires = fNumWires;
      float maxtime = fMaxTime;
      
      unsigned short vtxSizeIn = vtx.size();

  vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")<<"FindVertices plane "<<plane<<" pass "<<pass;
    PrintClusters();
  }

      float es1 = 0, eth1 = 0, ew1 = 0, et1 = 0;
      float bs1 = 0, bth1 = 0, bw1 = 0, bt1 = 0;
      float es2 = 0, eth2 = 0, ew2 = 0, et2 = 0;
      float bs2 = 0, bth2 = 0, bw2 = 0, bt2 = 0;
      float dth = 0,  dsl = 0, fvw = 0, fvt = 0;
//      float ec1 = 0, bc1 = 0, ec2 = 0, bc2 = 0, dchg = 0;
//      bool bothLong;
      float angcut = 0;
      unsigned short vw = 0, pass1, pass2;
      bool SigOK = false;
      for(unsigned short ii1 = 0; ii1 < sortindex.size() - 1; ++ii1) {
        unsigned short it1 = sortindex[ii1];
        // ignore obsolete clusters
        if(tcl[it1].ID < 0) continue;
        // ignore already attached clusters
        if(tcl[it1].BeginVtx >= 0 && tcl[it1].EndVtx >= 0) continue;
        es1 = tcl[it1].EndSlp;
        eth1 = tcl[it1].EndAng;
        ew1 = tcl[it1].EndWir;
        et1 = tcl[it1].EndTim;
//        ec1 = tcl[it1].EndChg;
        bs1 = tcl[it1].BeginSlp;
        bth1 = tcl[it1].BeginAng;
        bw1 = tcl[it1].BeginWir;
        bt1 = tcl[it1].BeginTim;
//        bc1 = tcl[it1].BeginChg;
        pass1 = tcl[it1].ProcCode - 10 *(tcl[it1].ProcCode / 10);
        for(unsigned short ii2 = ii1 + 1; ii2 < sortindex.size(); ++ii2) {
          unsigned short it2 = sortindex[ii2];
          if(tcl[it2].ID < 0) continue;
          // ignore already attached clusters
          if(tcl[it2].BeginVtx >= 0 && tcl[it2].EndVtx >= 0) continue;
          // try to attach cluster to existing vertices at either end
          ClusterVertex(it2);
          // ignore if both clusters are short
          if(tcl[it1].tclhits.size() < 5 &&
             tcl[it2].tclhits.size() < 5) continue;
          es2 = tcl[it2].EndSlp;
          eth2 = tcl[it2].EndAng;
          ew2 = tcl[it2].EndWir;
          et2 = tcl[it2].EndTim;
//          ec2 = tcl[it2].EndChg;
          bs2 = tcl[it2].BeginSlp;
          bth2 = tcl[it2].BeginAng;
          bw2 = tcl[it2].BeginWir;
          bt2 = tcl[it2].BeginTim;
//          bc2 = tcl[it2].BeginChg;
          pass2 = tcl[it2].ProcCode - 10 *(tcl[it2].ProcCode / 10);
          if(pass1 < pass2) {
            angcut = fKinkAngCut[pass2];
          } else {
            angcut = fKinkAngCut[pass1];
          }
//          bothLong = tcl[it1].tclhits.size() > 10 && tcl[it2].tclhits.size() > 10;
  // topo 1: check for vtx US of the ends of both clusters
          dth = std::abs(eth1 - eth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0 && dth > angcut) {
            dsl = es2 - es1;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            fvw = (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
            // vertex wire in the detector?
            if(fvw > 0. && fvw < nwires) {
              // require vtx in the range of wires with hits AND
              vw = (0.5 + fvw);
              // vtx US of both clusters AND
              // vtx not too far US of both clusters
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Chk clusters topo1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" fvw "<<fvw;
              if(vw >= fFirstWire - 1 && 
                 fvw <= ew1 + 1    && fvw <= ew2 + 1 &&
                 fvw > (ew1 - 10)  && fvw > (ew2 - 10) ) {
                fvt = et1 + (fvw - ew1) * es1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<" Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<"  vtx wire "<<vw<<" time "<<(int)fvt;
  }
                if(fvt > 0 && fvt < maxtime) {
                  // Vertex wire US of cluster ends and time in the detector.
                  // Check for signal at the vertex position and adjust the vertex by 1 wire
                  // if necessary
                  ChkSignal(vw, fvt, vw, fvt, SigOK);
                  if(!SigOK) {
                    fvw += 1.;
                    vw = (0.5 + fvw);
                    ChkSignal(vw, fvt, vw, fvt, SigOK);
                  }
                  // Check this against existing vertices and update
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<" SigOK "<<SigOK;
                  if(SigOK) ChkVertex(fvw, fvt, it1, it2, 1);
                } // fvt in detector
              } // vw topo 1 check
            } // fvw in detector
          } // topo 1
  // topo 2: check for vtx US of it1 and DS of it2
          dth = std::abs(eth1 - bth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0) {
            // large angle difference
            if(dth > angcut) {
              dsl = bs2 - es1;
              if(std::abs(ew1 - bw2) < 3 && std::abs(et1 - bt2) < 20) {
                fvw = 0.5 * (ew1 + bw2);
              } else {
                fvw = (et1 - ew1 * es1 - bt2 + bw2 * bs2) / dsl;
              }
              if(fvw > 0 && fvw < nwires) {
                // vertex wire in the detector
                vw = (0.5 + fvw);
                // require vtx US of cluster 1 End AND
                //         vtx DS of cluster 2 Begin 
                if(fvw <= ew1 + 1  && fvw >= bw2 - 1) {
                  fvt = et1 + (vw - ew1) * es1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk clusters topo2 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo2 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                  if(fvt > 0. && fvt < maxtime) {
                    ChkVertex(fvw, fvt, it1, it2, 2);
                  } // fvt in detector
                } // vw topo 2 check
              } // fvw in detector
            } // dth > 0.3
/*
            else if(bothLong && dth > 0.1) {
              // small angle difference. Require end-to-end clusters with
              // significant charge difference
              dchg = std::abs((ec1 - bc2) / (ec1 + bc2));
              if(ew1 > bw2 && ew1 < bw2 + 3 && std::abs(int(et1 - bt2)) < 10 && dchg > 0.5) {
                fvw = 0.5 *(ew1 + bw2);
                fvt = 0.5 *(et1 + bt2);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk clusters topo2 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo2 vtx wire "<<vw<<" time "<<(int)fvt<<" small angle";
  }
                ChkVertex(fvw, fvt, it1, it2, 2);
              }
            } // dth < 0.3
*/
          } // topo 2
  // topo 3: check for vtx DS of it1 and US of it2
          dth = std::abs(bth1 - eth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0) {
            // large angle difference
            if(dth > angcut) {
              dsl = bs1 - es2;
              if(std::abs(bw1 - ew2) < 3 && std::abs(bt1 - et2) < 20) {
                fvw = 0.5 * (bw1 + ew2);
              } else {
                fvw = (et2 - ew2 * es2 - bt1 + bw1 * bs1) / dsl;
              }
              if(fvw > 0 && fvw < nwires) {
                vw = (0.5 + fvw);
                // require vtx US of cluster 2 Begin AND
                //         vtx DS of cluster 1 End
                if(fvw <= ew2 + 1 && fvw >= bw1 - 1) {
                  fvt = et2 + (fvw - ew2) * es2;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk clusters topo3 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo3 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                  if(fvt > 0. && fvt < maxtime) {
                    ChkVertex(fvw, fvt, it1, it2, 3);
                  } // fvt in detector
                } // vw topo 3 check
              } // fvw in detector
            } 
/*
            else if(bothLong && dth > 0.1) {
              // small angle difference
              dchg = std::abs((bc1 - ec2) / (bc1 + ec2));
              if(ew2 > bw1 && ew2 < bw1 + 3 && std::abs(int(et1 - bt2)) < 10 && dchg > 0.6) {
                fvw = 0.5 *(ew2 + bw1);
                fvt = 0.5 *(et2 + bt1);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk clusters topo3 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo3 vtx wire "<<vw<<" time "<<(int)fvt<<" small angle";
  }
                ChkVertex(fvw, fvt, it1, it2, 3);
              }
            } // dth < 0.3
*/
          } // topo 3
  // topo 4: check for vtx DS of it1 and DS of it2
          dth = std::abs(bth1 - bth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0 && dth > angcut) {
            dsl = bs2 - bs1;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            // convert to integer if within the detector (vw, vt)
            fvw = (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              vw = (0.5 + fvw);
              // require vtx in the range of wires with hits AND
              // vtx DS of both clusters AND
              // vtx not too far DS of both clusters
              if(fvw <= fLastWire + 1 && 
                 fvw >= bw1 - 1 && fvw >= bw2 - 1 && 
                 fvw <  bw2 + 10 && fvw <  bw1 + 10) {
                fvt = bt1 + (fvw - bw1) * bs1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo4 vtx wire "<<vw<<" time "<<fvt;
  }
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check for signal at the vertex position and adjust the vertex by 1 wire
                  // if necessary
                  ChkSignal(vw, fvt, vw, fvt, SigOK);
                  if(!SigOK) {
                    fvw -= 1.;
                    vw = (0.5 + fvw);
                    ChkSignal(vw, fvt, vw, fvt, SigOK);
                  }
                  // Check this against existing vertices and update
                  if(SigOK) ChkVertex(fvw, fvt, it1, it2, 4);
                } // fvt in detector
              } // vw topo 4 check
            } // fvw in detector
          } // topo4
        } // it2
      } // it1

      // Fix vertices where both ends of a cluster are assigned to the 
      // same vertex. This can happen with short clusters
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != clCTP) continue;
        if(tcl[it].BeginVtx > -1 && tcl[it].BeginVtx == tcl[it].EndVtx ) {
          unsigned short iv = tcl[it].BeginVtx;
          float dwir = tcl[it].BeginWir - vtx[iv].Wire;
          float dtim = tcl[it].BeginTim - vtx[iv].Time;
          float rbeg = dwir * dwir + dtim * dtim;
          dwir = tcl[it].EndWir - vtx[iv].Wire;
          dtim = tcl[it].EndTim - vtx[iv].Time;
          float rend = dwir * dwir + dtim * dtim;
          if(rend < rbeg) {
            tcl[it].EndVtx = -99;
          } else {
            tcl[it].BeginVtx = -99;
          }
        } // tcl[it].BeginVtx == tcl[it].EndVtx
      } // it

      if(vtx.size() > vtxSizeIn) VtxWghtAndFit(clCTP);
      
    } // FindVertices()

/////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterVertex(unsigned short it)
    {
      // try to attach cluster it to an existing vertex
      
      if(vtx.size() == 0) return;
      
      unsigned short iv, jv;
      short dwib, dwjb, dwie, dwje;
      bool SigOK, matchUS, matchDS;
      float tdiff, ChiDOF;
      
      for(iv = 0; iv < vtx.size(); ++iv) {
        // ignore vertices in the wrong cryostat/TPC/Plane
        if(vtx[iv].CTP != clCTP) continue;
        // ignore deleted vertices
        if(vtx[iv].Wght < 0) continue;
        // determine which end to match - begin or end. Handle short tracks
        matchUS = false; matchDS = false;
        if(tcl[it].tclhits.size() < 6) {
          // See which end is closer to this vertex vs other vertices
          dwib = std::abs(vtx[iv].Wire - tcl[it].BeginWir);
          if(dwib > 2) dwib = 2;
          dwie = std::abs(vtx[iv].Wire - tcl[it].EndWir);
          if(dwie > 2) dwie = 2;
          dwjb = 999; dwje = 999;
          for(jv = 0; jv < vtx.size(); ++jv) {
            if(iv == jv) continue;
            if(std::abs(vtx[jv].Time - tcl[it].BeginTim) < 20) {
              if(std::abs(vtx[jv].Wire - tcl[it].BeginWir) < dwjb) 
                dwjb = std::abs(vtx[jv].Wire - tcl[it].BeginWir);
            } // std::abs(vtx[jv].Time - tcl[it].BeginTim) < 20
            if(std::abs(vtx[jv].Time - tcl[it].EndTim) < 20) {
              if(std::abs(vtx[jv].Wire - tcl[it].EndWir) < dwje) 
                dwje = std::abs(vtx[jv].Wire - tcl[it].EndWir);
            } // std::abs(vtx[jv].Time - tcl[it].EndTim) < 20
          } // jv
          matchUS = tcl[it].EndVtx != iv && 
                    dwie < 3 && dwie < dwje && dwie < dwib;
          matchDS = tcl[it].BeginVtx != iv && 
                    dwib < 3 && dwib < dwjb && dwib < dwie;
        } else {
          matchUS = tcl[it].EndVtx < 0   && vtx[iv].Wire <= tcl[it].EndWir + 2;
          matchDS = tcl[it].BeginVtx < 0 && vtx[iv].Wire >= tcl[it].BeginWir - 2;
        }
        if(matchUS) {
          // project cluster to US vertex
          tdiff = tcl[it].EndTim + (vtx[iv].Wire - tcl[it].EndWir) * tcl[it].EndSlp - vtx[iv].Time;
          if(std::abs(tdiff) < 10) {
            SigOK = false;
            ChkSignal(vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].EndVtx = iv;
              // re-fit it
              ChiDOF = 99.;
              FitVtx(iv, ChiDOF);
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Add End "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<ChiDOF;
              return;
            } // SigOK
          } // tdiff
        } // matchUS
        
        if(matchDS) {
          // project cluster to DS vertex
          tdiff = tcl[it].BeginTim + (vtx[iv].Wire - tcl[it].BeginWir) * tcl[it].BeginSlp - vtx[iv].Time;
          if(std::abs(tdiff) < 10) {
            SigOK = false;
            ChkSignal(vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].BeginVtx = iv;
              // re-fit it
              ChiDOF = 99.;
              FitVtx(iv, ChiDOF);
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Add Begin "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<ChiDOF;
              return;
            } // SigOK
          } // tdiff
        } // matchDS
        
      } // iv
    } // ClusterVertex()



/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkVertex(
        float fvw, float fvt, unsigned short it1, unsigned short it2, short topo)
      {
        // Checks the vertex (vw, fvt) against the existing set of vertices.
        // If there a match, clusters it1 and/or it2 are associated with it
        // if there is signal between the existing vertex and the start of
        // the cluster. The topo flag indicates the type of vertex that was
        // found: 1 = US of it1 && US of it2, 2 = US of it1 && DS of it2,
        // 3 = DS of it1 && US of it2, 4 = DS of it1 and DS of it2.
        // didit is set true if a cluster is attached to a (new or old) vertex
        
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"ChkVertex "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" topo "<<topo
      <<" fvw "<<fvw<<" fvt "<<fvt;
  }

        unsigned short vw = (unsigned short)(0.5 + fvw);
        unsigned short iht;

        // don't make vertices using short tracks that are close to
        // long straight clusters. These are most likely due to numerous
        // short delta ray clusters
        if(tcl[it1].tclhits.size() < 10 && tcl[it2].tclhits.size() < 10) {
          for(unsigned short it = 0; it < tcl.size(); ++it) {
            if(it == it1 || it == it2) continue;
            if(tcl[it].ID < 0) continue;
            if(tcl[it].CTP != clCTP) continue;
            if(tcl[it].tclhits.size() < 100) continue;
            if(std::abs(tcl[it].BeginAng - tcl[it].EndAng) > 0.1) continue;
            // don't reject because it is near the end of a long cluster
            if(vw <  tcl[it].EndWir + 5) continue;
            if(vw >  tcl[it].BeginWir - 5) continue;
            for(unsigned short ii = 0; ii < tcl[it].tclhits.size(); ++ii) {
              iht = tcl[it].tclhits[ii];
              if(fHits[iht].WireID().Wire <= vw) {
                if(std::abs(fHits[iht].PeakTime() - fvt) < 60) return;
              } // fHits[iht].WireID().Wire <= vWire
            } // ii
          } // it
        }

        bool SigOK = false;
        // check vertex and clusters for proximity to existing vertices
        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if(vtx[iv].CTP != clCTP) continue;
          if( std::abs(fvw - vtx[iv].Wire) < 4 &&
             std::abs(fvt - vtx[iv].Time) < 25) {
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Match to vtx "<<iv;
            // got a match. Check the appropriate cluster end and attach
            if( (topo == 1 || topo == 2) && tcl[it1].EndVtx < 0) {
              ChkSignal(vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, SigOK);
              if(SigOK) tcl[it1].EndVtx = iv;
  if(vtxprt)  mf::LogVerbatim("ClusterCrawlerAlg")<<"12 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv;
            } else if( (topo == 3 || topo == 4) && tcl[it1].BeginVtx < 0) {
              ChkSignal(vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, SigOK);
              if(SigOK) tcl[it1].BeginVtx = iv;
  if(vtxprt)  mf::LogVerbatim("ClusterCrawlerAlg")<<"34 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv;
            } // cluster it2
            if( (topo == 1 || topo == 3) && tcl[it2].EndVtx < 0) {
              ChkSignal(vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, SigOK);
              if(SigOK) tcl[it2].EndVtx = iv;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"13 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv;
            } else if( (topo == 2 || topo == 4) && tcl[it2].BeginVtx < 0) {
              ChkSignal(vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, SigOK);
              if(SigOK) tcl[it2].BeginVtx = iv;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"24 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv;
            } // cluster it2
            return;
          } // matched vertex
        } // iv

        // no match to existing vertices. Ensure that there is a wire signal between
        // the vertex and the appropriate ends of the clusters
        bool Sig1OK = false;
        bool Sig2OK = false;
        if(topo == 1 || topo == 2) {
          ChkSignal(vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, Sig1OK);
        } else {
          ChkSignal(vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, Sig1OK);
        }
        if(topo == 1 || topo == 3) {
          ChkSignal(vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, Sig2OK);
        } else {
          ChkSignal(vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, Sig2OK);
        }
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk new "<<tcl[it1].ID<<" OK "<<Sig1OK<<" "<<tcl[it2].ID<<" OK "<<Sig2OK
      <<" Vtx at "<<vw<<" "<<(int)fvt;
  }
        // both clusters must have an OK signal
        if(Sig1OK && Sig2OK) {
          VtxStore newvx;
          newvx.Wire = vw;
          newvx.Time = fvt;
          newvx.Wght = 1;
          newvx.Topo = topo;
          newvx.CTP = clCTP;
          vtx.push_back(newvx);
          unsigned short iv = vtx.size() - 1;
          if(topo == 1 || topo == 2) {
            tcl[it1].EndVtx = iv;
          } else {
            tcl[it1].BeginVtx = iv;
          }
          if(topo == 1 || topo == 3) {
            tcl[it2].EndVtx = iv;
          } else {
            tcl[it2].BeginVtx = iv;
          }
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"New vtx "<<iv<<" in plane "<<plane<<" topo "<<topo<<" cls "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" time "<<(int)fvt<<" wire "<<vw;
  }
        }

      } // ChkVertex()

/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkSignal
      (unsigned short wire1, float time1, unsigned short wire2, float time2, bool& SigOK)
    {
      // returns SigOK true if there is a signal on the line between
      // (wire1, time1) and (wire2, time2).       SigOK = false;
      // Be sure to set time1 < time 2 if you are checking for signals on a single wire
      
      // get the begin and end right
      short wireb = wire1;
      float timeb = time1;
      short wiree = wire2;
      float timee = time2;
      // swap them?
      if(wiree > wireb) {
        wireb = wire2;
        timeb = time2;
        wiree = wire1;
        timee = time1;
      }
      if(wiree < fFirstWire || wiree > fLastWire) return;
      if(wireb < fFirstWire || wireb > fLastWire) return;
      short wire0 = wiree;
      // checking a single wire?
      float slope = 0;
      if(wireb == wiree) {
        slope = 0.;
      } else {
        slope = (timeb - timee) / (wireb - wiree);
      }
      float prTimeLo = 0;
      float prTimeHi = 0;
      for(short wire = wiree; wire < wireb + 1; ++wire) {
        // assume there is no signal on this wire
        bool WireSigOK = false;
        // checking a single wire?
        if(wireb == wiree) {
          prTimeLo = time1;
          prTimeHi = time2;
        } else {
          prTimeLo = timee + (wire - wire0) * slope;
          prTimeHi = prTimeLo;
        }
        unsigned short index = wire - fFirstWire;
        // skip dead wires
        if(WireHitRange[index].first == -1) continue;
        // no hits on this wire
        if(WireHitRange[index].first == -2) {
          SigOK = false;
          return;
        }
        unsigned short firsthit = WireHitRange[index].first;
        unsigned short lasthit = WireHitRange[index].second;
//  mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkSignal "<<wiree<<" "<<wire<<" "<<wireb<<" "<<(int)prtime;
//    <<" first "<<firsthit<<" last "<<lasthit;
        for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
          // outside the hit RAT?
          if(prTimeHi > fHits[khit].EndTick()) continue;
          if(prTimeLo < fHits[khit].StartTick()) continue;
          // found a signal. Skip checking on this wire
          WireSigOK = true;
          break;
        } // khit
//        mf::LogVerbatim("ClusterCrawlerAlg")<<" SigOK "<<WireSigOK;
        if(!WireSigOK) return;
      } // wire
      SigOK = true;
    } // ChkSignal()

/////////////////////////////////////////
    void ClusterCrawlerAlg::SplitCluster
      (unsigned short icl, unsigned short pos, unsigned short ivx)
    {
      // split cluster icl into two clusters

      if(tcl[icl].ID < 0) {
        mf::LogError("ClusterCrawlerAlg")<<"Trying to split obsolete cluster "
          <<icl;
        return;
      }

      if(pos > tcl[icl].tclhits.size()) {
        mf::LogError("ClusterCrawlerAlg")<<"SplitCluster bad split position "
          <<pos<<" in cluster "<<tcl[icl].ID<<" size "<<tcl[icl].tclhits.size();
        return;
      }

      // Create the first cluster (DS) using the Begin info
      clBeginSlp = tcl[icl].BeginSlp;
      clBeginSlpErr = tcl[icl].BeginSlpErr;
      clBeginAng = tcl[icl].BeginAng;
      clBeginWir = tcl[icl].BeginWir;
      clBeginTim = tcl[icl].BeginTim;
      clBeginChg = tcl[icl].BeginChg;
      clStopCode = 5;
      clProcCode = tcl[icl].ProcCode;
      fcl2hits.clear();
      chifits.clear();
      hitnear.clear();
      for(unsigned short ii = 0; ii < pos; ++ii) {
        unsigned short iht = tcl[icl].tclhits[ii];
        fcl2hits.push_back(iht);
      }
      // determine the pass in which this cluster was created
      pass = tcl[icl].ProcCode - 10 * (tcl[icl].ProcCode / 10);
      if(pass > fNumPass-1) pass = fNumPass-1;
      // fit the end hits
      FitCluster();
      clEndSlp = clpar[1];
      clEndSlpErr = clparerr[1];
      clEndAng = std::atan(fScaleF * clEndSlp);
      // find the charge at the end
      FitClusterChg();
      clEndChg = fAveChg;
      TmpStore();
      // associate the End with the supplied vertex
      unsigned short iclnew = tcl.size() - 1;
      tcl[iclnew].EndVtx = ivx;
      tcl[iclnew].BeginVtx = tcl[icl].BeginVtx;

      // now create the second cluster (US)
      clEndSlp = tcl[icl].EndSlp;
      clEndSlpErr = tcl[icl].EndSlpErr;
      clEndAng = std::atan(fScaleF * clEndSlp);
      clEndWir = tcl[icl].EndWir;
      clEndTim = tcl[icl].EndTim;
      clEndChg = tcl[icl].EndChg;
      clStopCode = 5;
      clProcCode = tcl[icl].ProcCode;
      fcl2hits.clear();
      chifits.clear();
      hitnear.clear();
      bool didFit = false;
      for(unsigned short ii = pos; ii < tcl[icl].tclhits.size(); ++ii) {
        unsigned short iht = tcl[icl].tclhits[ii];
  if(fHitInCluster[iht] != tcl[icl].ID) {
    mf::LogError("ClusterCrawlerAlg")
      <<"SplitCluster bad hit "<<iht<<" "<<fHitInCluster[iht]
      <<" "<<tcl[icl].ID<<" ProcCode "<<clProcCode<<std::endl;
  }
        fcl2hits.push_back(iht);
        // define the Begin parameters
        if(fcl2hits.size() == fMaxHitsFit[pass] ||
           fcl2hits.size() == fMinHits[pass]) {
          FitCluster();
          clBeginSlp = clpar[1];
          clBeginAng = std::atan(fScaleF * clBeginSlp);
          clBeginSlpErr = clparerr[1];
          didFit = true;
        }
        if((unsigned short)fcl2hits.size() == fNHitsAve[pass] + 1) {
          FitClusterChg();
          clBeginChg = fAveChg;
          didFit = true;
        }
      } // ii
      // do a fit using all hits if one wasn't done
      if(!didFit) {
        FitCluster();
        FitClusterChg();
        clBeginChg = fAveChg;
      }
      TmpStore();
      // associate the End with the supplied vertex
      iclnew = tcl.size() - 1;
      tcl[iclnew].BeginVtx = ivx;
      tcl[iclnew].EndVtx = tcl[icl].EndVtx;
//  mf::LogVerbatim("ClusterCrawlerAlg")<<"ClusterSplit split "<<tcl[icl].ID;
      // declare icl obsolete
      tcl[icl].ID = -tcl[icl].ID;

    } // SplitCluster()

/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkMerge()
    {
      // Try to merge clusters. Clusters that have been subsumed in other
      // clusters, i.e. no longer valid, have ID < 0
      
      if(tcl.size() < 2) return;
      // The size of the ClusterStore vector will increase while merging
      // is on-going so the upper limit on it1 is fixed tcl.size() - 1 
      // before merging starts

      prt = (fDebugPlane == (short)plane && fDebugWire < 0);
      
      unsigned short it1, it2, nh1, pass1, pass2;
      float bs1, bth1, bt1, bc1, es1, eth1, et1, ec1;
      float bs2, bth2, bt2, bc2, es2, eth2, et2, ec2;
      unsigned short bw1, ew1, bw2, ew2;
      float dth, angcut, chgrat, chgcut, dtim, timecut;
      bool bothLong, NoVtx, SigOK;
      
      unsigned short skipcut, tclsize = tcl.size();
      
      for(it1 = 0; it1 < tclsize - 1; ++it1) {
        // ignore already merged clusters
        if(tcl[it1].ID < 0) continue;
        if(tcl[it1].CTP != clCTP) continue;
        bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        bth1 = std::atan(fScaleF * bs1);
        // more compact notation for begin/end, wire/time/chg/slp/theta, 1/2
        bw1 = tcl[it1].BeginWir;
        bt1 = tcl[it1].BeginTim;
        bc1 = tcl[it1].BeginChg;
        es1 = tcl[it1].EndSlp;
        eth1 = tcl[it1].EndAng;
        ew1 = tcl[it1].EndWir;
        et1 = tcl[it1].EndTim;
        ec1 = tcl[it1].EndChg;
        nh1 = tcl[it1].tclhits.size();
        pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
        if(pass1 > fNumPass) pass1 = fNumPass;
        for(it2 = it1 + 1; it2 < tclsize; ++it2) {
          // ignore already merged clusters
          if(tcl[it1].ID < 0) continue;
          if(tcl[it2].ID < 0) continue;
          // only merge if they are in the right cryostat/TPC/plane
          if(tcl[it2].CTP != clCTP) continue;
          bs2 = tcl[it2].BeginSlp;
          bth2 = std::atan(fScaleF * bs2);
          bw2 = tcl[it2].BeginWir;
          bt2 = tcl[it2].BeginTim;
          bc2 = tcl[it2].BeginChg;
          es2 = tcl[it2].EndSlp;
          eth2 = tcl[it2].EndAng;
          ew2 = tcl[it2].EndWir;
          et2 = tcl[it2].EndTim;
          ec2 = tcl[it2].EndChg;
          pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
          if(pass2 > fNumPass) pass2 = fNumPass;
          // use the more promiscuous pass for cuts
          angcut = fKinkAngCut[pass1];
          if(fKinkAngCut[pass2] > angcut) angcut = fKinkAngCut[pass2];
          skipcut = fMaxWirSkip[pass1];
          if(fMaxWirSkip[pass2] > skipcut) skipcut = fMaxWirSkip[pass2];
          chgcut = fMergeChgCut[pass1];
          if(fMergeChgCut[pass2] > chgcut) chgcut = fMergeChgCut[pass2];
          timecut = fTimeDelta[pass];
          if(fTimeDelta[pass2] > timecut) timecut = fTimeDelta[pass2];
          // increase the time cut for large angle clusters
          timecut *= (2 - 1/(1 + std::abs(clpar[1])));

          // tweak the cuts for long straight tracks
          bothLong = (nh1 > 100 && tcl[it2].tclhits.size() > 100);
          
          // look for US and DS broken clusters w similar angle.
          // US cluster 2 merge with DS cluster 1?
          // This is the most likely occurrence given the order in which
          // clusters are created so put it first.
          dth = std::abs(bth2 - eth1);
          // require no vertex between
          NoVtx = (tcl[it1].EndVtx < 0) && (tcl[it2].BeginVtx < 0);
  if(prt && bw2 < ew1 ) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(ew1 - bw2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          if(NoVtx && bw2 < ew1 && (ew1 - bw2)  < skipcut && dth < angcut) {
            chgrat = 2 * std::abs(bc2 - ec1) / (bc2 + ec1);
            // ignore the charge cut for large angle clusters
            if(std::abs(es1) > fLAClusSlopeCut) chgrat = 0.;
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.05) chgrat = 0.;
            // project bw2,bt2 to ew1
            dtim = std::abs(bt2 + (ew1-bw2)*bs2 - et1);
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<" dtim "<<dtim<<" timecut "<<(int)timecut
    <<" ec1 "<<ec1<<" bc2 "<<bc2
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              // ensure there is a signal between cluster ends
              SigOK = true;
              ChkSignal(ew1,et1,bw2,bt2,SigOK);
              if(SigOK) {
                DoMerge(it2, it1, 10);
                tclsize = tcl.size();
                break;
              }
            }
          } // US cluster 2 merge with DS cluster 1?
          
          // look for US and DS broken clusters w similar angle
          // US cluster 1 merge with DS cluster 2?
          dth = std::abs(bth1 - eth2);
  if(prt && bw1 < ew2 && (ew2 - bw1)  < skipcut) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<"Chk2 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(bw1 - ew2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          // require no vertex between
          NoVtx = (tcl[it2].EndVtx < 0) && (tcl[it1].BeginVtx < 0);
          if(NoVtx && bw1 < ew2 && (ew2 - bw1)  < skipcut && dth < angcut ) {
            chgrat = 2 * std::abs((bc1 - ec2) / (bc1 + ec2));
            // ignore the charge cut for large angle clusters
            if(std::abs(es2) > fLAClusSlopeCut) chgrat = 0.;
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.05) chgrat = 0.;
            // project sw1,st1 to ew2
            dtim = std::abs(bt1 + (ew2-bw1)*bs1 - et2);
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<" dtim "<<dtim<<" err "<<dtim<<" timecut "<<(int)timecut
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              SigOK = true;
              ChkSignal(bw1,bt1,ew2,et2,SigOK);
              if(SigOK) {
                DoMerge(it1, it2, 10);
                tclsize = tcl.size();
                break;
              }
            }
          } // US cluster 1 merge with DS cluster 2
          
          if(bw2 < bw1 && ew2 > ew1) {
            // look for small cl2 within the wire boundary of cl1
            // with similar times and slopes for both clusters
            dth = std::abs(eth2 - eth1);
            dtim = std::abs(et1 +(ew2 - ew1 - 1)*es1 - et2);
            // count the number of wires with no hits on cluster 1
            short nmiss1 = bw1 - ew1 + 1 - tcl[it1].tclhits.size();
            // compare with the number of hits in cluster 2
            short nin2 = tcl[it2].tclhits.size();
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<"cl2: "<<tcl[it2].ID<<" within cl1 "<<tcl[it1].ID
    <<" ? dth "<<dth<<" dtim "<<dtim<<" nmissed "<<nmiss1;
  }
            // make rough cuts before calling ChkMerge12
            // this may not work well for long wandering clusters
            bool didit = false;
            if(dth < 1 && dtim < timecut && nmiss1 >= nin2) 
                ChkMerge12(it1, it2, didit);
            if(didit) {
              tclsize = tcl.size();
              break;
            } //didit
          } // small cl2 within the wire boundary of cl1
          
          if(bw1 < bw2 && ew1 > ew2) {
            // look for small cl1 within the wire boundary of cl2
            // with similar times and slopes for both clusters
            dth = std::abs(eth2 - eth1);
            dtim = std::abs(et2 +(ew1 - ew2 - 1)*es2 - et1);
            // count the number of wires with no hits on cluster 2
            short nmiss2 = bw2 - ew2 + 1 - tcl[it2].tclhits.size();
            // compare with the number of hits in cluster 1
            short nin1 = tcl[it1].tclhits.size();
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<"cl1: "<<tcl[it1].ID<<" within cl2 "<<tcl[it2].ID
    <<" ? dth "<<dth<<" dtim "<<dtim<<" nmissed "<<nmiss2;
  }
            // make rough cuts before calling ChkMerge12
            // this may not work well for long wandering clusters
            bool didit = false;
            if(dth < 1 && dtim < timecut && nmiss2 >= nin1) 
                ChkMerge12(it2, it1, didit);
            if(didit) {
              tclsize = tcl.size();
              break;
            } // didit
          } // small cl1 within the wire boundary of cl2
          
          if(tcl[it1].ID < 0) break;
        } // cluster 2
        if(tcl[it1].ID < 0) continue;
      } // cluster 1
    }

/////////////////////////////////////////
  void ClusterCrawlerAlg::ChkMerge12
    (unsigned short it1, unsigned short it2, bool& didit)
  {
    // Calling routine has done a rough check that cluster it2 is a candidate
    // for merging with cluster it1. The wire range spanned by it2 lies
    // within the wire range of it1 and the clusters are reasonably close
    // together in time.
    
    // assume that no merging was done
    didit = false;
    
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkMerge12 "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    
    ClusterStore& cl1 = tcl[it1];
    // fill a vector spanning the length of cluster 1 and filled with the hit time
    unsigned short ew1 = tcl[it1].EndWir;
    unsigned short bw1 = tcl[it1].BeginWir;
    std::vector<unsigned short> cl1hits;
    // fill the vector with 0s
    for(unsigned short wire = ew1; wire <= bw1; ++wire) {
      cl1hits.push_back(0);
    }
    // put in the hit IDs
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned short hit = cl1.tclhits[iht];
      unsigned short wire = fHits[hit].WireID().Wire;
      if(wire - ew1 < 0 || wire - ew1 > (short)cl1hits.size()) {
        mf::LogError("ClusterCrawlerAlg")<<"ChkMerge12 bad wire "<<(wire-ew1);
        return;
      }
      cl1hits[wire - ew1] = hit;
    }
    unsigned short ew2 = tcl[it2].EndWir;
    float et2 = tcl[it2].EndTim;
    // look for the closest wire with a hit on cluster 1
    unsigned short wiron1 = 0;
    // count the number of missing hits
    short nmiss = 0;
    for(unsigned short wire = ew2 - 1; wire > ew1; --wire) {
      if(cl1hits[wire - ew1] > 0) {
        wiron1 = wire;
        break;
      }
      ++nmiss;
    } // wire
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"chk next US wire "<<wiron1<<" missed "<<nmiss;
    if(wiron1 == 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    
    // compare the wires with hits on cluster 2 with the gap in cluster 1
    // the number of hit wires that fit in the gap
    unsigned short nfit = 0;
    for(unsigned short iht = 0; iht < tcl[it2].tclhits.size(); ++iht) {
      unsigned short hiton2 = tcl[it2].tclhits[iht];
      unsigned short wiron2 = fHits[hiton2].WireID().Wire;
      if(wiron2 < ew1 || wiron2 > bw1) return;
      if(cl1hits[wiron2 - ew1] == 0) ++nfit;
    }
    // require complete filling of the gap
    if(nfit < tcl[it2].tclhits.size()) return;
    
    // decode the pass for both clusters and select the matching cuts
    unsigned short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
    unsigned short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
    unsigned short cpass = pass1;
    // use the tighter cuts
    if(pass2 < pass1) cpass = pass2;
    
    // ***** Check End of Cluster 2 matching with middle of cluster 1
    if(wiron1 - ew1 < 0) return;
    unsigned short hiton1 = cl1hits[wiron1 - ew1];
    if(hiton1 > fHits.size() - 1) {
      mf::LogError("ClusterCrawlerAlg")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    // check the time difference
    float timon1 = fHits[hiton1].PeakTime();
    float dtim = std::abs(et2 + (wiron1 - ew2) * tcl[it2].EndSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    // check the slope difference. First do a local fit on cluster 1 near
    // the matching point
    FitClusterMid(it1, hiton1, 3);
    if(clChisq > 20.) return;
    // fit parameters are now in clpar.
    // check for angle consistency
    float dth = std::abs(std::atan(fScaleF * clpar[1]) - std::atan(fScaleF * tcl[it2].EndSlp));
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"US dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut. fAveChg was calculated in FitClusterMid
    float chgrat = 2 * std::abs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"US chgrat "<<chgrat<<" cut "<<fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    bool SigOK = false;
    ChkSignal(wiron1, timon1, ew2, et2, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"US SigOK? "<<SigOK;
    if( !SigOK ) return;


    // ***** Check Begin of Cluster 2 matching with middle of cluster 1
    unsigned short bw2 = tcl[it2].BeginWir;
    float bt2 = tcl[it2].BeginTim;
    nmiss = 0;
    wiron1 = 0;
    for(unsigned short wire = bw2 + 1; wire < bw1; ++wire) {
      if(cl1hits[wire - ew1] > 0) {
        wiron1 = wire;
        break;
      }
      ++nmiss;
    }
    if(wiron1 == 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    // fit this section of cluster 1 with 4 hits starting at the hit on the
    // closest wire and moving DS
    hiton1 = cl1hits[wiron1 - ew1];
    if(hiton1 > fHits.size() - 1) {
      mf::LogError("ClusterCrawlerAlg")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    timon1 = fHits[hiton1].PeakTime();
    dtim = std::abs(bt2 - (wiron1 - bw2) *tcl[it2].BeginSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    FitClusterMid(it1, hiton1, -3);
    if(clChisq > 20.) return;
    // check for angle consistency
    dth = std::abs(std::atan(fScaleF * clpar[1]) - std::atan(fScaleF * tcl[it2].BeginSlp));
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
      <<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * std::abs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
      <<"DS chgrat "<<chgrat<<" cut "<<fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    SigOK = false;
    ChkSignal(wiron1, timon1, bw2, bt2, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"DS SigOK? "<<SigOK;
    if( !SigOK ) return;

  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Merge em";
    // success. Merge them
    DoMerge(it1, it2, 100);
    didit = true;
  } // ChkMerge12()

/////////////////////////////////////////
  void ClusterCrawlerAlg::DoMerge
    (unsigned short it1, unsigned short it2, short inProcCode)
  {
    // Merge clusters.
    
    ClusterStore& cl1 = tcl[it1];
    ClusterStore& cl2 = tcl[it2];
/*    
    bool myprt = false;
    if(plane == 1 && cl1.ID == 84 && cl2.ID == 63) myprt = true;
*/
    // mark cl1 and cl2 obsolete
    cl1.ID = -cl1.ID;
    cl2.ID = -cl2.ID;
    
    // Find the low and high wire for both clusters.
    // Assume that cluster 1 is DS
    bool cl1DS = true;
    unsigned short hiwire = cl1.BeginWir;
    if(cl2.BeginWir > hiwire) {
      hiwire = cl2.BeginWir;
      cl1DS = false;
    }
    unsigned short lowire = cl1.EndWir;
    if(cl2.EndWir < lowire) lowire = cl2.EndWir;

    // make a vector of wire hits
    std::vector<short> wirehit;
    for(unsigned short wire = lowire; wire < hiwire + 2; ++wire) {
      wirehit.push_back(-1);
    }
    // put in the hit IDs for cluster 2
    for(unsigned short iht = 0; iht < cl2.tclhits.size(); ++iht) {
      unsigned short hit = cl2.tclhits[iht];
      // un-assign the hit
      fHitInCluster.setFree(hit);
      unsigned short wire = fHits[hit].WireID().Wire;
      unsigned short index = wire - lowire;
      wirehit[index] = hit;
/*
  if(myprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Cl2 hit "<<wire<<":"<<(int)fHits[hit].PeakTime()
    <<" wire index "<<index<<" hit index "<<hit;
*/
     } // iht
    // now cluster 1
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned short hit = cl1.tclhits[iht];
      unsigned short wire = fHits[hit].WireID().Wire;
      unsigned short index = wire - lowire;
      fHitInCluster.setFree(hit);
      // TODO: Should merge hits instead of clobbering one of them?
      wirehit[index] = hit;
/*
  if(myprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Cl1 hit "<<wire<<":"<<(int)fHits[hit].PeakTime()
    <<" wire index "<<index<<" hit index "<<hit;
*/
    } // iht
    // make the new cluster
    fcl2hits.clear();
    for(std::vector<short>::reverse_iterator rit = wirehit.rbegin(); 
        rit != wirehit.rend(); ++rit) {
      if(*rit < 0) continue;
      unsigned short hit = *rit;
      fcl2hits.push_back(hit);
    } // rit

    short endVtx = 0;
    short begVtx = 0;
    short del1Vtx = -99;
    short del2Vtx = -99;
    if(cl1DS) {
      // use cluster 1 Begin info
      clBeginSlp = cl1.BeginSlp;
      clBeginSlpErr = cl1.BeginSlpErr;
      clBeginAng = cl1.BeginAng;
      clBeginWir = cl1.BeginWir;
      clBeginTim = cl1.BeginTim;
      clBeginChg = cl1.BeginChg;
      begVtx = cl1.BeginVtx;
      del1Vtx = cl1.EndVtx;
      // and cluster 2 End info
      clEndSlp = cl2.EndSlp;
      clEndSlpErr = cl2.EndSlpErr;
      clEndAng = cl2.EndAng;
      clEndWir = cl2.EndWir;
      clEndTim = cl2.EndTim;
      clEndChg = cl2.EndChg;
      endVtx = cl2.EndVtx;
      del2Vtx = cl2.BeginVtx;
      clStopCode = cl2.StopCode;
    } else {
      // use cluster 2 Begin info
      clBeginSlp = cl2.BeginSlp;
      clBeginSlpErr = cl2.BeginSlpErr;
      clBeginAng = cl2.BeginAng;
      clBeginWir = cl2.BeginWir;
      clBeginTim = cl2.BeginTim;
      clBeginChg = cl2.BeginChg;
      begVtx = cl2.BeginVtx;
      del2Vtx = cl2.EndVtx;
      // and cluster 1 End info
      clEndSlp = cl1.EndSlp;
      clEndSlpErr = cl1.EndSlpErr;
      clEndWir = cl1.EndWir;
      clEndTim = cl1.EndTim;
      clEndChg = cl1.EndChg;
      endVtx = cl1.EndVtx;
      del1Vtx = cl1.BeginVtx;
      clStopCode = cl1.StopCode;
    }

    // append it to the tcl vector
    TmpStore();
    unsigned short itnew = tcl.size()-1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"DoMerge "<<cl1.ID<<" "<<cl2.ID<<" -> "<<tcl[itnew].ID;
    // stuff the processor code with the current pass
    tcl[itnew].ProcCode = inProcCode + pass;
    // transfer the vertex info
    // delete a vertex between these two?
    if(del1Vtx >= 0 && del1Vtx == del2Vtx) vtx[del1Vtx].Wght = -1;
    // preserve the vertex assignments
    tcl[itnew].BeginVtx = begVtx;
    tcl[itnew].EndVtx = endVtx;
  } // DoMerge

/////////////////////////////////////////
  void ClusterCrawlerAlg::PrintClusters()
  {

    float aveRMS = 0;
    // prints clusters to the screen for code development
    mf::LogVerbatim myprt("ClusterCrawlerAlg");
    myprt<<"  ID CTP nht Stop  Proc  beg_W:T  bTheta"
      <<" begChg end_W:T  eTheta eChg  bVx  eVx aveRMS\n";
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      if(fDebugPlane >= 0 && fDebugPlane != (int) DecodeCTP(tcl[ii].CTP).Plane) continue;
      myprt<<std::right<<std::setw(4)<<tcl[ii].ID;
      myprt<<std::right<<std::setw(3)<<tcl[ii].CTP;
      myprt<<std::right<<std::setw(5)<<tcl[ii].tclhits.size();
      myprt<<std::right<<std::setw(4)<<tcl[ii].StopCode;
      myprt<<std::right<<std::setw(6)<<tcl[ii].ProcCode;
      unsigned short iTime = tcl[ii].BeginTim;
      myprt<<std::right<<std::setw(6)<<tcl[ii].BeginWir<<":"<<iTime;
      if(iTime < 10) {
        myprt<<"   ";
      } else if(iTime < 100) {
        myprt<<"  ";
      } else if(iTime < 1000) myprt<<" ";
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].BeginAng;
      myprt<<std::right<<std::setw(5)<<(short)tcl[ii].BeginChg;
      iTime = tcl[ii].EndTim;
      myprt<<std::right<<std::setw(6)<<tcl[ii].EndWir<<":"<<iTime;
      if(iTime < 10) {
        myprt<<"   ";
      } else if(iTime < 100) {
        myprt<<"  ";
      } else if(iTime < 1000) myprt<<" ";
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].EndAng;
      myprt<<std::right<<std::setw(5)<<(short)tcl[ii].EndChg;
      myprt<<std::right<<std::setw(5)<<tcl[ii].BeginVtx;
      myprt<<std::right<<std::setw(5)<<tcl[ii].EndVtx;
      aveRMS = 0;
      unsigned short iht = 0;
      for(unsigned short jj = 0; jj < tcl[ii].tclhits.size(); ++jj) {
        iht = tcl[ii].tclhits[jj];
        aveRMS += fHits[iht].RMS();
      }
      aveRMS /= (float)tcl[ii].tclhits.size();
      myprt<<std::right<<std::setw(5)<<std::fixed<<std::setprecision(1)<<aveRMS;
      myprt<<"\n";
    } // ii
    // print out vertices
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      mf::LogVerbatim("ClusterCrawlerAlg")
        <<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time
        <<" wght "<<(int)vtx[iv].Wght<<" topo "<<vtx[iv].Topo;
    }    
    // Check for incompatible hit->cluster cluster->hit associations
    for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if(tcl[icl].ID < 0) continue;
      for(unsigned short ii = 0; ii< tcl[icl].tclhits.size(); ++ii) {
        unsigned short iht = tcl[icl].tclhits[ii];
        if(fHitInCluster[iht] != tcl[icl].ID) {
          mf::LogVerbatim("PrintClusters")
            <<"Association error: cluster "<<tcl[icl].ID
            <<" Hit "<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()
            <<" InClus is incorrect "<<fHitInCluster[iht]
            <<" Hit index "<<iht;
          return;
        }
      } // ii
    } // icl

  } // PrintClusters()

/////////////////////////////////////////
    void ClusterCrawlerAlg::TmpGet(unsigned short it1)
    {
      // copies temp cluster it1 into the fcl2hits vector, etc. This is 
      // effectively the inverse of cl2TmpStore
      
      if(it1 > tcl.size()) return;


      clBeginSlp = tcl[it1].BeginSlp;
      clBeginSlpErr = tcl[it1].BeginSlpErr;
      clBeginAng = tcl[it1].BeginAng;
      clBeginWir = tcl[it1].BeginWir;
      clBeginTim = tcl[it1].BeginTim;
      clBeginChg = tcl[it1].BeginChg;
      clEndSlp = tcl[it1].EndSlp;
      clEndSlpErr = tcl[it1].EndSlpErr;
      clEndAng = tcl[it1].EndAng;
      clEndWir = tcl[it1].EndWir;
      clEndTim = tcl[it1].EndTim;
      clEndChg = tcl[it1].EndChg;
      clStopCode = tcl[it1].StopCode;
      clProcCode = tcl[it1].ProcCode;
      clCTP = tcl[it1].CTP;
      fcl2hits = tcl[it1].tclhits;
    }


/////////////////////////////////////////
  void ClusterCrawlerAlg::TmpStore()
  {

    if(fcl2hits.size() < 3) return;
    
    if(fcl2hits.size() == USHRT_MAX) return;
    
    ++NClusters;

    // flag all the hits as used
    for(unsigned short it = 0; it < fcl2hits.size(); ++it) {
      unsigned short hit = fcl2hits[it];
      if(!isHitPresent(hit)) {
        mf::LogError("ClusterCrawlerAlg")<<"Trying to use obsolete hit "<<hit
        <<" on wire "<<fHits[hit].WireID().Wire<<" on cluster "<<NClusters
        <<" in plane "<<plane<<" ProcCode "<<clProcCode;
      }
      fHitInCluster.setCluster(hit, NClusters);
    }

    // ensure that the cluster begin/end info is correct

    // define the begin/end charge if it wasn't done already
    if(clEndChg <= 0) {
      // use the next to the last two hits. The End hit may have low charge
      unsigned int ih0 = fcl2hits.size() - 2;
      unsigned int hit = fcl2hits[ih0];
      clEndChg = fHits[hit].Integral();
      hit = fcl2hits[ih0 - 1];
      clEndChg += fHits[hit].Integral();
      clEndChg = clEndChg / 2.;
    }
    if(clBeginChg <= 0) {
      // use the 2nd and third hit. The Begin hit may have low charge
      unsigned int hit = fcl2hits[1];
      clBeginChg = fHits[hit].Integral();
      hit = fcl2hits[2];
      clBeginChg += fHits[hit].Integral();
      clBeginChg = clBeginChg / 2.;
    }
    
    std::vector<unsigned short>::const_iterator ibg = fcl2hits.begin();
    unsigned short hitb = *ibg;
    std::vector<unsigned short>::const_iterator iend = fcl2hits.end() - 1;
    unsigned short hite = *iend;

    // store the cluster in the temporary ClusterStore struct
    ClusterStore clstr;
    
    clstr.ID = NClusters;
    clstr.BeginSlp    = clBeginSlp;
    clstr.BeginSlpErr = clBeginSlpErr;
    clstr.BeginAng    = std::atan(fScaleF * clBeginSlp);
    clstr.BeginWir    = fHits[hitb].WireID().Wire;
    clstr.BeginTim    = fHits[hitb].PeakTime();
    clstr.BeginChg    = clBeginChg;
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndAng      = std::atan(fScaleF * clEndSlp);
    clstr.EndWir      = fHits[hite].WireID().Wire;
    clstr.EndTim      = fHits[hite].PeakTime();
    clstr.EndChg      = clEndChg;
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.BeginVtx    = -99;
    clstr.EndVtx      = -99;
    clstr.CTP         = clCTP;
    clstr.tclhits     = fcl2hits;
    tcl.push_back(clstr);
  } // TmpStore()

/////////////////////////////////////////
  void ClusterCrawlerAlg::LACrawlUS() {
    // Crawl a large angle cluster upstream. Similar to CrawlUS but require
    // that a hit be added on each wire


    unsigned short dhit = fcl2hits[0];
    short dwir = fHits[dhit].WireID().Wire;
    prt = false;
  if(fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
    prt = std::abs(fHits[dhit].PeakTime() - fDebugHit) < 20;

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawlerAlg");
    myprt<<"******************* LACrawlUS PASS "<<pass<<" Hits ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned short iht = fcl2hits[fcl2hits.size() - 1 - ii];
      myprt<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" ";
    }
  }

    // Merge all the hits if this is a very large angle cluster
    if(std::abs(clpar[1]) > 5) {
      unsigned short ii, nmult = 0;
      for(ii = 0; ii < fcl2hits.size(); ++ii) 
        if(fHits[fcl2hits[ii]].Multiplicity() > 1) ++nmult;
      if(nmult == fcl2hits.size()) {
        for(ii = 0; ii < fcl2hits.size(); ++ii) {
          MergeHits(fcl2hits[ii]);
        } // ii
        FitCluster();
        clBeginSlp = clpar[1];
      } // nmult == fcl2hits.size()
    } // pass == fNumPass - 1

    bool SigOK = true;
    bool HitOK = true;
    short nmissed = 0;
    // count the number of kinks encountered. Hits US of the kink are removed
    // and crawling continues unless another kink is encountered
    unsigned short kinkOnWire = 0;
    unsigned short it = fcl2hits.size() - 1;
    unsigned short lasthit = fcl2hits[it];
    unsigned short lastwire = fHits[lasthit].WireID().Wire;
    bool ChkCharge = false;
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"LACrawlUS: next wire "<<nextwire;
      // stop crawling if there is a nearby vertex
      if(CrawlVtxChk(nextwire)) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"LACrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // AddLAHit will merge the hit on nextwire if necessary
      AddLAHit(nextwire, ChkCharge, HitOK, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"LACrawlUS: HitOK "<<HitOK<<" SigOK "<<SigOK;
      if(!SigOK) break;
      if(!HitOK) {
        ++nmissed;
        if(nmissed > fMaxWirSkip[pass]) {
          clStopCode = 1;
          break;
        }
        continue;
      }
      // Merge all of the hit multiplets in the fcl2hits array into single
      // hits when enough hits have been found to call this a credible large
      // angle cluster. The last hit was already merged in AddHit
      if(fcl2hits.size() == 4) {
        for(unsigned short kk = 0; kk< fcl2hits.size()-1; ++kk) {
          unsigned short hit = fcl2hits[kk];
          MergeHits(hit);
        }
        // update the fit
        FitCluster();
        clBeginSlp = clpar[1];
        // start checking the charge ratio when adding new hits
        ChkCharge = true;
        continue;
      } // fcl2hits.size() == 4
      unsigned short chsiz = chifits.size()-1;
      // chsiz is fcl2hits.size() - 1...
      if(chsiz < 6) continue;
      if(fKinkChiRat[pass] <= 0) continue;
      if(chifits.size() != fcl2hits.size()) {
        mf::LogError("ClusterCrawlerAlg")
          <<"LACrawlUS: chifits size error "<<chifits.size()<<" "<<fcl2hits.size();
        return;
      }
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
      if( chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
          chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
        // find the kink angle (crudely) from the 0th and 2nd hit
        unsigned short ih0 = fcl2hits.size() - 1;
        unsigned short hit0 = fcl2hits[ih0];
        unsigned short ih2 = ih0 - 2;
        unsigned short hit2 = fcl2hits[ih2];
        float dt02 = fHits[hit2].PeakTime() - fHits[hit0].PeakTime();
        float dw02 = fHits[hit2].WireID().Wire - fHits[hit0].WireID().Wire;
        float th02 = std::atan( fScaleF * dt02 / dw02);
        // and the 3rd and 5th hit
        unsigned short ih3 = ih2 - 1;
        unsigned short hit3 = fcl2hits[ih3];
        unsigned short ih5 = ih3 - 2;
        unsigned short hit5 = fcl2hits[ih5];
        float dt35 = fHits[hit5].PeakTime() - fHits[hit3].PeakTime();
        float dw35 = fHits[hit5].WireID().Wire - fHits[hit3].WireID().Wire;
        float th35 = std::atan(fScaleF * dt35 / dw35);
        float dth = std::abs(th02 - th35);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Kink angle "<<std::setprecision(3)<<dth<<" cut "<<fKinkAngCut[pass];
        if(dth > fKinkAngCut[pass]) {
          // hit a kink. Lop of the first 3 hits, refit and keep crawling?
          for(short jj = 0; jj < 3; ++jj) {
            fcl2hits.pop_back();
            chifits.pop_back();
            hitnear.pop_back();
          }
          FitCluster();
          // See if this is a second kink and it is close to the first
          // kink (which had hits removed).
          if(kinkOnWire > 0) {
            if(kinkOnWire - nextwire < 4) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Hit a second kink. kinkOnWire = "<<kinkOnWire<<" Stopping";
              // set the kink stop code
              clStopCode = 3;
              break;
            }
          }
          kinkOnWire = nextwire;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Removed kink hits";
        } // kinkang check
      } // chifits test
    } // nextwire

    CheckClusterHitFrac(prt);

    clProcCode += 300;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"LACrawlUS done. Nhits = "<<fcl2hits.size();
    prt = false;
  } // LACrawlUS

/////////////////////////////////////////
  void ClusterCrawlerAlg::CrawlUS()
  {
    // Crawl along a trail of hits moving upstream

    if(fcl2hits.size() < 2) return;

    unsigned short dhit = fcl2hits[0];
    short dwir = fHits[dhit].WireID().Wire;
    prt = false;
  if(fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
    prt = std::abs(fHits[dhit].PeakTime() - fDebugHit) < 20;

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawlerAlg");
    myprt<<"******************* CrawlUS PASS "<<pass<<" Hits: ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned short iht = fcl2hits[fcl2hits.size() - 1 - ii];
      myprt<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" ";
    }
    myprt<<"\n";
  }

    // SigOK = true if there is a ADC signal near the projected cluster position
    bool SigOK = true;
    bool HitOK = true;
    // count the number of missed hits on adjacent wires
    short nmissed = 0;
    // count the number of added hits after skipping
    short nHitAfterSkip = 0;
    bool DidaSkip = false;
    bool PostSkip = false;
    unsigned short it = fcl2hits.size() - 1;
    unsigned short lasthit = fcl2hits[it];
    if(lasthit > fHits.size() - 1) {
      mf::LogError("ClusterCrawlerAlg")<<"CrawlUS bad lasthit "<<lasthit;
    }
    unsigned short lastwire = fHits[lasthit].WireID().Wire;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"CrawlUS: last wire "<<lastwire<<" hit "<<lasthit;
    
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"CrawlUS: next wire "<<nextwire;
      // stop crawling if there is a nearby vertex
      if(CrawlVtxChk(nextwire)) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"CrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // add hits and check for PH and width consistency
      AddHit(nextwire, HitOK, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"CrawlUS: HitOK "<<HitOK<<" SigOK "<<SigOK;
      if(!HitOK) {
        // no hit on this wire. Was there a signal or dead wire?
        if(SigOK) {
        // no hit on the wire but there is a signal
          ++nmissed;
          // stop if too many missed wires
          if(nmissed > fMaxWirSkip[pass]) {
            clStopCode = 1;
            break;
          }
          // see if we are in the PostSkip phase and missed more than 1 wire
// is this an error?
//          if(PostSkip && nmissed > 1) {
          if(PostSkip && nmissed > fMinWirAfterSkip[pass]) {
            // cluster is really short
            if((short)(fcl2hits.size() - nHitAfterSkip) < 4) {
              fcl2hits.clear();
              return;
            }
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" PostSkip && nmissed = "<<nmissed;
            clStopCode = 2;
            for(short jj = 0; jj < nHitAfterSkip; ++jj) {
              fcl2hits.pop_back();
              chifits.pop_back();
              hitnear.pop_back();
            } // pop_back
            FitCluster();
            if(clChisq > 90.) {
              fcl2hits.clear();
              return;
            }
            FitCluster();
            return;
          } // PostSkip && nmissed > 
          if(nmissed > 1) {
            DidaSkip = true;
            PostSkip = false;
          }
        } // SigOK
        else {
          // SigOK is false
          clStopCode = 0;
          if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"No hit or signal on wire "<<nextwire;
          break;
        } // else SigOK false
      } // !HitOK
      else {
        if(clChisq > 99.) {
          if(fcl2hits.size() < 3) return;
          // a fit error occurred. Lop off the leading hit and refit
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Fit failed ";
          fcl2hits.pop_back();
          chifits.pop_back();
          hitnear.pop_back();
          FitCluster();
          if(clChisq > 99.) {
            // something really bad happened. Bail out
            fcl2hits.clear();
            return;
          }
          FitClusterChg();
          continue;
        } // clChisq > 99
        // monitor the onset of a kink. Look for a progressive increase
        // in chisq for the previous 0 - 2 hits.
        if(chifits.size() > 5 && fKinkChiRat[pass] > 0) {
          if(chifits.size() != fcl2hits.size()) {
            mf::LogError("ClusterCrawlerAlg")
              <<"CrawlUS: chifits size error "<<chifits.size()<<" "<<fcl2hits.size();
            return;
          }
          unsigned short chsiz = chifits.size()-1;
  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
          if( chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
              chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
            if(fcl2hits.size() != chifits.size()) {
              mf::LogError("ClusterCrawlerAlg")
              <<"bad kink check size "<<chifits.size()<<" "<<fcl2hits.size()
              <<" plane "<<plane<<" cluster "<<dwir<<":"<<dhit;
              continue;
            }
            // find the kink angle (crudely) from the 0th and 2nd hit
            unsigned short ih0 = fcl2hits.size() - 1;
            unsigned short hit0 = fcl2hits[ih0];
            unsigned short ih2 = ih0 - 2;
            unsigned short hit2 = fcl2hits[ih2];
            float dt02 = fHits[hit2].PeakTime() - fHits[hit0].PeakTime();
            float dw02 = fHits[hit2].WireID().Wire - fHits[hit0].WireID().Wire;
            float th02 = std::atan( fScaleF * dt02 / dw02);
            // and the 3rd and 5th hit
            unsigned short ih3 = ih2 - 1;
            unsigned short hit3 = fcl2hits[ih3];
            unsigned short ih5 = ih3 - 2;
            unsigned short hit5 = fcl2hits[ih5];
            float dt35 = fHits[hit5].PeakTime() - fHits[hit3].PeakTime();
            float dw35 = fHits[hit5].WireID().Wire - fHits[hit3].WireID().Wire;
            float th35 = std::atan(fScaleF * dt35 / dw35);
            float dth = std::abs(th02 - th35);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Kink angle "<<std::setprecision(3)<<dth<<" cut "<<fKinkAngCut[pass];
            // cut on the allowed kink angle
            if(dth > fKinkAngCut[pass]) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Stopped tracking ";
              // kill the last 3 hits, refit and return
              for(short jj = 0; jj < 3; ++jj) {
                fcl2hits.pop_back();
                chifits.pop_back();
                hitnear.pop_back();
              }
              FitCluster();
              FitClusterChg();
              // set the kink stop code but keep looking
              clStopCode = 3;
            } // kinkang check
          } // chifits check
        } // chifits.size() > 5
        // done with kink check
        // update the cluster Begin information?
        if(fcl2hits.size() == fMaxHitsFit[pass] ||
           fcl2hits.size() == fMinHits[pass]) {
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
        }
        // define the begin cluster charge if it's not defined yet
        if(clBeginChg <= 0 && fAveChg > 0) {
          clBeginChg = fAveChg;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Set clBeginChg "<<clBeginChg;
        }
        // reset nmissed
        nmissed = 0;
        // start counting hits added after skipping
        if(DidaSkip) {
          // start PostSkip phase
          PostSkip = true;
          DidaSkip = false;
          nHitAfterSkip = 0;
        } // DidaSkip
        // check for PostSkip phase
        if(PostSkip) {
          // end the PostSkip phase if there are enough hits
          ++nHitAfterSkip;
          if(nHitAfterSkip == fMinWirAfterSkip[pass]) PostSkip = false;
        } 
        // check for bad chisq
        if(clChisq > fChiCut[pass]) {
          // remove the last few hits if there is a systematic increase in chisq and re-fit
          // long tracks only
          for(unsigned short nlop = 0; nlop < 4; ++nlop) {
            unsigned short cfsize = chifits.size() - 1;
            float chirat = chifits[cfsize] / chifits[cfsize - 1];
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"chirat "<<chirat
    <<" last hit "<<fcl2hits[fcl2hits.size()-1];
            if(chirat < 1.2) break;
            fcl2hits.pop_back();
            chifits.pop_back();
            hitnear.pop_back();
            if(fcl2hits.size() < 6) break;
            if(chifits.size() < 6) break;
          } // nlop
          if(fcl2hits.size() < 6) {
            clStopCode = 4;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Bad fit chisq - short cluster. Break";
            break;
          }
          FitCluster();
          FitClusterChg();
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Bad fit chisq - removed hits. Continue...";
        } // clChisq > fChiCut[pass]
      } // !HitOK check
    } // nextwire

    // count the number of hits on adjacent wires at the leading edge and
    // ensure that the count is consistent with fMinWirAfterSkip
    bool reFit = false;
    if((unsigned short)fcl2hits.size() > fMinWirAfterSkip[pass] + 3) {
      unsigned short ih0 = fcl2hits.size() - 1;
      unsigned short hit0 = fcl2hits[ih0];
      unsigned short uswir = fHits[hit0].WireID().Wire;
      unsigned short nAdjHit = 0;
      for(unsigned short ii = ih0 - 1; ii > 0; --ii) {
        unsigned short nxtwir = fHits[ fcl2hits[ii] ].WireID().Wire;
        if(nxtwir != uswir + 1) break;
        ++nAdjHit;
        // break if there are enough hits
        if( nAdjHit == fMinWirAfterSkip[pass] ) break;
        uswir = nxtwir;
      } // ii
      // lop off hits?
      if(nAdjHit < fMinWirAfterSkip[pass]) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Lop off "<<nAdjHit<<" hits ";
        for(unsigned short ii = 0; ii < nAdjHit + 1; ++ii) {
          fcl2hits.pop_back();
          chifits.pop_back();
          hitnear.pop_back();
          reFit = true;
        }
      }
    } // fcl2hits.size() > fMinWirAfterSkip[pass] + 3
    if(reFit) {
      FitCluster();
      FitClusterChg();
    }
    
    CheckClusterHitFrac(prt);

    prt = false;
  } // CrawlUS()

/////////////////////////////////////////
  void ClusterCrawlerAlg::CheckClusterHitFrac(bool prt)
  {


    // Find the fraction of the wires on the cluster that have
    // hits.
    unsigned short iht = fcl2hits[fcl2hits.size() - 1];
    clEndWir = fHits[iht].WireID().Wire;
    float hitFrac = (float)fcl2hits.size() / (float)(clBeginWir - clEndWir + 1);

    if(hitFrac < 0.7) {
      fcl2hits.clear();
      if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
        <<"CheckClusterHitFrac: Poor hit fraction "<<hitFrac;
      return;
    } // hitFrac
    
    // lop off the last hit if it is part of a hit multiplet
    if(fHits[iht].Multiplicity() > 1) {
      fcl2hits.resize(fcl2hits.size() - 1);
    }

    // check for short track ghosts
    if(fcl2hits.size() > 4) return;
    unsigned short nsing = 0;
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      iht = fcl2hits[ii];
      if(fHits[iht].Multiplicity() == 1) ++nsing;
    } // ii
    hitFrac = (float)nsing / (float)fcl2hits.size();
    
    if(hitFrac < 0.7) {
      fcl2hits.clear();
      if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
        <<"CheckClusterHitFrac: Poor short track hit fraction "<<hitFrac;
      return;
    } // hitFrac
        
  } // CheckClusterHitFrac()

/////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterMid
    (unsigned short it1, unsigned short ihtin, short nhit)
  {
    // Fits hits on temp cluster it1 to a line starting at hit ihtin and including
    // nhit hits incrementing towards the hit vector End when nhit > 0 and
    // decrementing towards the hit vector Begin when nhit < 0.
    // The fit params are stashed in the clpar and clparerr arrays. 
    // fAveChg is re-calculated as well.
    
    
    // set chisq bad in case something doesn't work out
    clChisq = 99.;
    
    ClusterStore cls = tcl[it1];
    if(cls.tclhits.size() < 3) return;

    std::vector<float> xwir;
    std::vector<float> ytim;
    std::vector<float> ytimerr2;
    
    short nht = 0;
    unsigned short wire0 = 0;
    if(nhit > 0) {
      nht = nhit;
      // find the first desired hit and move towards the End
      fAveChg = 0.;
      fChgSlp = 0.;
      unsigned short hitcnt = 0;
      bool UseEm = false;
      for(unsigned short it = 0; it < cls.tclhits.size(); ++it) {
        unsigned short ihit = cls.tclhits[it];
        if(ihit > fHits.size()-1) {
          mf::LogError("ClusterCrawlerAlg")<<"FitClusterMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = fHits[ihit].WireID().Wire;
        }
        // use hits after finding the first desired hit
        if(UseEm) {
          unsigned short wire = fHits[ihit].WireID().Wire;
          xwir.push_back(wire - wire0);
          ytim.push_back(fHits[ihit].PeakTime());
          // pass the error^2 to the fitter
          float terr = fHitErrFac * fHits[ihit].RMS();
	  ytimerr2.push_back(terr * terr);
          fAveChg += fHits[ihit].Integral();
          ++hitcnt;
          if(hitcnt == nht) break;
        }
      }    
      nht = hitcnt;
    } else {
      nht = -nhit;
      // find the first desired hit and move towards the Begin
      fAveChg = 0.;
      fChgSlp = 0.;
      unsigned short hitcnt = 0;
      bool UseEm = false;
      for(auto it = cls.tclhits.crbegin(); it != cls.tclhits.crend(); ++it) {
        unsigned short ihit = *it;
        if(ihit > fHits.size()-1) {
          mf::LogVerbatim("ClusterCrawlerAlg")<<"FitClusterMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = fHits[ihit].WireID().Wire;
        }
        // use hits after finding the first desired hit

        if(UseEm) {
          unsigned short wire = fHits[ihit].WireID().Wire;
          xwir.push_back(wire - wire0);
          ytim.push_back(fHits[ihit].PeakTime());
          float terr = fHitErrFac * fHits[ihit].RMS();
	  ytimerr2.push_back(terr * terr);
          fAveChg += fHits[ihit].Integral();
          ++hitcnt;
          if(hitcnt == nht) break;
        }
      }    
      nht = hitcnt;
    }
    
    if(nht < 2) return;
    fAveChg = fAveChg / (float)nht;
    fChgSlp = 0.;
    
    float intcpt = 0.;
    float slope = 0.;
    float intcpterr = 0.;
    float slopeerr = 0.;
    float chidof = 0.;
    fLinFitAlg.LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    clChisq = chidof;
    if(clChisq > 99.) return;
    clpar[0] = intcpt;
    clpar[1] = slope;
    clparerr[0] = intcpterr;
    clparerr[1] = slopeerr;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::FitCluster()
  {
    // Fits the hits on the US end of a cluster. This routine assumes that
    // wires are numbered from lower (upstream) to higher (downstream) and
    // that the hits in the fclhits vector are sorted so that upstream hits
    // are at the end of the vector


    clChisq = 999.;
    
    unsigned short nht = 0;
    // fit all hits or truncate?
    if(fcl2hits.size() < fMaxHitsFit[pass]) {
      nht = fcl2hits.size();
    } else {
      nht = fMaxHitsFit[pass];
    }
    if(nht < 2) return;

    std::vector<float> xwir;
    std::vector<float> ytim;
    std::vector<float> ytimerr2;
    // apply an angle dependent scale factor.
    float angfactor = AngleFactor(clpar[1]);

    // load the hits starting at the End of the fcl2hits vector.
    // These are the most upstream hits
    unsigned short iht = 0;

    bool first = true;
    unsigned short wire0 = 0;
    for(std::vector<unsigned short>::reverse_iterator it = fcl2hits.rbegin();
       it != fcl2hits.rend(); ++it) {
      unsigned short ihit = *it;
      unsigned short wire = fHits[ihit].WireID().Wire;
      if(first) {
        wire0 = wire;
        first = false;
      }
      xwir.push_back(wire - wire0);
      ytim.push_back(fHits[ihit].PeakTime());
      // Scale error by hit multiplicity to account for bias in hit
      // multiplet fitting
      float terr = fHitErrFac * fHits[ihit].RMS() * fHits[ihit].Multiplicity();
      ytimerr2.push_back(angfactor * terr * terr);
      if(iht == nht) break;
      ++iht;
    }

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawlerAlg");
    myprt<<"FitCluster W:T ";
    unsigned short cnt = 0;
    for(std::vector<unsigned short>::reverse_iterator it = fcl2hits.rbegin();
       it != fcl2hits.rend(); ++it) {
      unsigned short ihit = *it;
      unsigned short wire = fHits[ihit].WireID().Wire;
      myprt<<wire<<":"<<(short)fHits[ihit].PeakTime()<<" ";
      ++cnt;
      if(cnt == 8) {
        myprt<<" .... ";
        break;
      }
    }
  } // prt
    
    nht = iht;
    
    if(nht < 2) return;

    float intcpt = 0.;
    float slope = 0.;
    float intcpterr = 0.;
    float slopeerr = 0.;
    float chidof = 0.;
    fLinFitAlg.LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    clChisq = chidof;
    if(chidof > 99.) return;
    clpar[0] = intcpt;
    clpar[1] = slope;
    clparerr[0] = intcpterr;
    clparerr[1] = slopeerr;

  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"nht "<<nht<<" fit par "<<(int)clpar[0]<<"+/-"<<(int)intcpterr
    <<" "<<clpar[1]<<"+/-"<<slopeerr
    <<" clChisq "<<clChisq;
  
  }
/////////////////////////////////////////
  float ClusterCrawlerAlg::AngleFactor(float slope)
  {
    // returns an angle dependent cluster projection error factor for fitting
    // and hit finding
    
    float slp = std::abs(slope);
    if(slp > 15) slp = 15;
    // return a value between 1 and 4
    float angfac = 1 + 0.03 * slp * slp;
    return angfac;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterChg()
  {
    // Fits the charge of hits on the fcl2hits vector to a line, or simply
    // uses the average of 1 or 2 hits as determined by NHitsAve

    unsigned short ih0 = fcl2hits.size() - 1;
    
    if(pass >= fNumPass) {
      mf::LogError("ClusterCrawlerAlg")<<"FitClusterChg bad pass "<<pass;
      return;
    }
    
    // number of hits at the leading edge that we will fit
    unsigned short fitLen = fNHitsAve[pass];
    // start fitting charge when there are at least 6 hits if we are tracking
    // long clusters
    if(fitLen > 5 &&  // Fit 6 hits when tracking long clusters AND
       fcl2hits.size() > 5 && // there are at least 6 hits AND
       fcl2hits.size() < fitLen) // there are less than fNHitsAve[pass]
        fitLen = 5;
    
    // don't find the average charge --> no charge cut is made
    if(fNHitsAve[pass] < 1) return;
    
    if(fNHitsAve[pass] == 1) {
      // simply use the charge and width the last hit
      fAveChg = fHits[fcl2hits[ih0]].Integral();
      fChgSlp = 0.;
    } else if(fNHitsAve[pass] == 2) {
      // average the last two points if requested
      fAveChg = (fHits[fcl2hits[ih0]].Integral() + 
                 fHits[fcl2hits[ih0 - 1]].Integral()) / 2.;
      fChgSlp = 0.;
    } else if((unsigned short)fcl2hits.size() > fitLen){
      // do a real fit
      std::vector<float> xwir;
      std::vector<float> ychg;
      std::vector<float> ychgerr2;
      // origin of the fit
      unsigned short wire0 = fHits[fcl2hits[fcl2hits.size()-1]].WireID().Wire;
      // find the mean and rms of the charge
      unsigned short npt = 0;
      unsigned short imlast = 0;
      float ave = 0.;
      float rms = 0.;
      // this loop intentionally ignores the Begin hit
      for(unsigned int ii = fcl2hits.size() - 1; ii > 0; --ii) {
        ++npt;
        float chg = fHits[fcl2hits[ii]].Integral();
        ave += chg;
        rms += chg * chg;
        if(npt == fitLen) {
          imlast = ii;
          break;
        }
      }
      float fnpt = npt;
      ave /= fnpt;
      rms = std::sqrt((rms - fnpt * ave * ave) / (fnpt - 1));
      float chgcut = ave + rms;
      for(unsigned int ii = fcl2hits.size() - 1; ii > imlast; --ii) {
        unsigned short wire = fHits[fcl2hits[ii]].WireID().Wire;
        float chg = fHits[fcl2hits[ii]].Integral();
        if(chg > chgcut) continue;
        xwir.push_back((float)(wire - wire0));
        ychg.push_back(chg);
        ychgerr2.push_back(chg);
      }
      if(ychg.size() < 3) return;
      float intcpt; float slope; float intcpterr;
      float slopeerr; float chidof;
      fLinFitAlg.LinFit(xwir, ychg, ychgerr2, intcpt, slope, intcpterr, slopeerr, chidof);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"FitClusterChg wire "<<wire0
    <<" chidof "<<(int)chidof<<" npt "<<xwir.size()
    <<" charge = "<<(int)intcpt<<" slope = "<<(int)slope
    <<" first ave "<<(int)ave<<" rms "<<(int)rms;
      if(chidof > 100.) return;
      // fit must have gone wrong if the truncated average is greater than
      // the average using all points
      if(intcpt > ave) return;
      // ensure that change does not exceed 30%
      if(fAveChg > 0) {
        ave = intcpt / fAveChg;
        if(ave > 1.3) return;
        if(ave < 0.77) return;
      }
      fAveChg = intcpt;
      fChgSlp = slope;
    }
  } // fitchg

/////////////////////////////////////////
  void ClusterCrawlerAlg::AddLAHit
    (unsigned short kwire, bool& ChkCharge, bool& HitOK, bool& SigOK)
  {
    // A variant of AddHit for large angle clusters
    
    SigOK = false;
    HitOK = false;
    
    // not in the range of wires with hits
    if(kwire < fFirstWire || kwire > fLastWire) return;
    
    if(fcl2hits.size() == 0) return;

    unsigned short index = kwire - fFirstWire;
    // return if no signal and no hit
    if(WireHitRange[index].first == -2) return;
    // skip bad wire, but assume the track was there
    if(WireHitRange[index].first == -1) {
      SigOK = true;
      return;
    }
    unsigned short firsthit = WireHitRange[index].first;
    unsigned short lasthit = WireHitRange[index].second;
    
    // max allowable time difference between projected cluster and a hit
    float timeDiff = 40 * AngleFactor(clpar[1]);
    float dtime;
    
    // the last hit added to the cluster
    unsigned short lastClHit = fcl2hits[fcl2hits.size()-1];
    unsigned short wire0 = fHits[lastClHit].WireID().Wire;
    // the projected time of the cluster on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    float chgrat;
    
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"AddLAHit: wire "<<kwire<<" prtime "<<(int)prtime
    <<" max time diff "<<timeDiff;
    unsigned short imbest = 0;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(!isHitPresent(khit)) continue;
      chgrat = fHits[khit].Integral() / fHits[lastClHit].Integral();
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" Chk W:T "<<kwire<<":"<<(short)fHits[khit].PeakTime()
    <<" Charge "<<(short)fHits[khit].Integral()
    <<" chgrat "<<std::setw(8)<<std::setprecision(2)<<chgrat
    <<" InClus "<<fHitInCluster[khit]
    <<" mult "<<fHits[khit].Multiplicity()
    <<" RMS "<<std::setprecision(2)<<fHits[khit].RMS()
    <<" Chi2 "<<std::setw(8)<<std::setprecision(2)<<fHits[khit].GoodnessOfFit()
    <<" LoT "<<(int)fHits[khit].StartTick()
    <<" HiT "<<(int)fHits[khit].EndTick();
      // projected time outside the Signal time window?
      if(prtime < fHits[khit].StartTick() - 20) continue;
      if(prtime > fHits[khit].EndTick() + 20) continue;
      SigOK = true;
      // hit used?
      if(isHitInCluster(khit)) continue;
      // ignore very low charge hits
      if(chgrat < 0.1) continue;
      // ignore very high charge hits for intermediate angle clusters
      if(std::abs(clpar[1]) < 5 && chgrat > 2) continue;
      dtime = std::abs(prtime - fHits[khit].PeakTime());
      // Use a tighter requirement for crude hits that may span a very large
      // time range. Crude hits have a large chisq and/or large rms
      if((fHits[khit].GoodnessOfFit() > 500. || fHits[khit].RMS() > 20.) &&
          dtime > 20) continue;
      if(dtime < timeDiff) {
        HitOK = true;
        imbest = khit;
        timeDiff = dtime;
      }
    } // khit
    
  if(prt) {
    if(!HitOK) mf::LogVerbatim("ClusterCrawlerAlg")<<" no hit found ";
  }
    if(!HitOK) return;

  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" Pick hit time "<<(int)fHits[imbest].PeakTime()
    <<" hit index "<<imbest;
    
    // merge the hits in a multiplet?
    bool doMerge = false;
    // assume there is no nearby hit
    short hnear = 0;
    if(fHits[imbest].Multiplicity() > 1) {
      doMerge = true;
      // don't merge if we are close to a vertex
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(vtx[ivx].CTP != clCTP) continue;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" close vtx chk W:T "<<vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time;
        if(std::abs(kwire - vtx[ivx].Wire) < 5 &&
           std::abs(int(fHits[imbest].PeakTime() - vtx[ivx].Time)) < 20 ) {
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" Close to a vertex. Don't merge hits";
          doMerge = false;
        }
      } // ivx
      // Decide which hits in the multiplet to merge. Hits that are well
      // separated from each other should not be merged
      if(doMerge) {
        unsigned short nused = 0;
        // the total charge of the hit multiplet
        float multipletChg = 0.;
        // make an angle dependent chisq cut that varies between
        // 3 sigma at 45 degrees and 6 sigma at 90 degrees
        float chicut = 6 * (1 - 1/(1 + std::abs(clpar[1])));
        // look for a big separation between adjacent hits
        std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(imbest);
        for(size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
          // ignore obsolete hits
          if(!isHitPresent(jht)) continue;
          // count used hits
          if(isHitFree(jht)) multipletChg += fHits[jht].Integral();
          else ++nused;
          // check the neighbor hit separation
          if(jht > MultipletRange.first) {
            // pick the larger RMS of the two hits
            // TODO use std::max()
            float hitRMS = fHits[jht].RMS();
            if(fHits[jht - 1].RMS() > hitRMS) hitRMS = fHits[jht-1].RMS();
            const float tdiff = std::abs(fHits[jht].PeakTime() - fHits[jht-1].PeakTime()) / hitRMS;
            if(tdiff > chicut) doMerge = false;
          } // jht > 0
        } // jht
  if(prt) {
    if(!doMerge) mf::LogVerbatim("ClusterCrawlerAlg")
      <<" Hits are well separated. Don't merge them";
  }
        if(doMerge && nused == 0) {
          // compare the charge with the last hit added?
          if(ChkCharge) {
            // there is a nearby hit
            hnear = 1;
            float chgrat = multipletChg / fHits[lastClHit].Integral();
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" merge hits charge check "
    <<(int)multipletChg<<" Previous hit charge "<<(int)fHits[lastClHit].Integral();
            if(chgrat > 2) doMerge = false;
          }
        } // doMerge && nused == 0
      } // doMerge true
      if(doMerge) {
        // there is a nearby hit and it will be merged
        hnear = -1;
        MergeHits(imbest);
      }
    } // fHits[imbest].Multiplicity() > 1
    
    // attach to the cluster and fit
    fcl2hits.push_back(imbest);
    FitCluster();
    chifits.push_back(clChisq);
    hitnear.push_back(hnear);
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" >>ADD W:T "<<kwire<<":"<<(int)fHits[imbest].PeakTime()
    <<std::setprecision(3)<<" clChisq "<<clChisq
    <<" charge "<<(int)fHits[imbest].Integral();
    // decide what to do with a bad fit
    if(clChisq > fChiCut[pass]) {
      fcl2hits.pop_back();
      chifits.pop_back();
      hitnear.pop_back();
      FitCluster();
      HitOK = false;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" Bad fit. Removed hit. New clChisq "
    <<std::setprecision(3)<<clChisq
    <<" nhits "<<fcl2hits.size();
    }
      // stop tracking if previous chisq values were close to the cut.
      // this is an indicator that the track is wandering too much for this pass
      if(chifits.size() > 2 && 
        chifits[chifits.size()-2] > 0.8 * fChiCut[pass]) SigOK = false;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"  Set SigOK = "<<SigOK;

  } // AddLAHit()


/////////////////////////////////////////
  void ClusterCrawlerAlg::AddHit(unsigned short kwire, bool& HitOK, bool& SigOK)
  {
    // Add a hit to the cluster if it meets several criteria:
    // similar pulse height to the cluster (if fAveChg is defined)
    // closest hit to the project cluster position.
    // Return SigOK if there is a nearby hit that was missed due to the cuts
    
    SigOK = false;
    HitOK = false;
    
    // not in the range of wires with hits
    if(kwire < fFirstWire || kwire > fLastWire) return;

    // the last hit added to the cluster
    unsigned short lastClHit = fcl2hits[fcl2hits.size()-1];
    unsigned short wire0 = fHits[lastClHit].WireID().Wire;

    unsigned short index = kwire - fFirstWire;
    // return if no signal and no hit
    if(fAllowNoHitWire == 0) {
      if(WireHitRange[index].first == -2) return;
    } else {
      // allow a number of wires with no hits
      if(WireHitRange[index].first == -2 && 
        (wire0 - kwire) > fAllowNoHitWire) {
        SigOK = true;
        return;
      }
    }
    // skip bad wire, but assume the track was there
    if(WireHitRange[index].first == -1) {
      SigOK = true;
      return;
    }

    unsigned short firsthit = WireHitRange[index].first;
    unsigned short lasthit = WireHitRange[index].second;

    
    // the projected time of the cluster on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    // Find the projected time error including the projection error and the
    // error from the last hit added
    float prtimerr2 = std::abs(kwire-wire0)*clparerr[1]*clparerr[1];
    // apply an angle dependent scale factor to the hit error
    float hiterr = AngleFactor(clpar[1]) * fHitErrFac * fHits[lastClHit].RMS();
    float err = std::sqrt(prtimerr2 + hiterr * hiterr);
    // Time window for accepting a hit.
    float prtimeLo = prtime - 4 * err;
    float prtimeHi = prtime + 4 * err;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"AddHit: wire "<<kwire
    <<" prtime Lo "<<(int)prtimeLo<<" Hi "<<(int)prtimeHi
    <<" fAveChg "<<(int)fAveChg;

    // loop through the hits
    size_t imbest = HitInCluster_t::InvalidHitIndex; // invalid by default
    float best = 9999.;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(!isHitPresent(khit)) continue;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" Chk W:T "<<kwire<<":"<<(short)fHits[khit].PeakTime()
    <<" InClus "<<fHitInCluster[khit]
    <<" mult "<<fHits[khit].Multiplicity()
    <<" RMS "<<std::setprecision(2)<<fHits[khit].RMS()
    <<" Chi2 "<<std::setprecision(2)<<fHits[khit].GoodnessOfFit()
    <<" Charge "<<(int)fHits[khit].Integral()
    <<" LoT "<<(int)fHits[khit].StartTick()
    <<" HiT "<<(int)fHits[khit].EndTick();
      // check for signal
      if(prtime < fHits[khit].StartTick()) continue;
      if(prtime > fHits[khit].EndTick()) continue;
      SigOK = true;
      // check for good hit
      if(fHits[khit].PeakTime() < prtimeLo) continue;
      if(fHits[khit].PeakTime() > prtimeHi) continue;
      // hit used?
      if(isHitInCluster(khit)) continue;
      float dtime = std::abs(fHits[khit].PeakTime() - prtime);
      if(dtime < best) {
        best = dtime;
        imbest = khit;
      }
    } // khit
    
    if(!SigOK) {
      if(fAllowNoHitWire == 0) return;
      if((wire0 - kwire) > fAllowNoHitWire) return;
      SigOK = true;
    }

    if(imbest == HitInCluster_t::InvalidHitIndex) return;

    recob::Hit const& hit = fHits[imbest];
    
    if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
      <<" Best hit time "<<(int)hit.PeakTime();

    // merge hits in a doublet?
    // assume there is no nearby hit
    short hnear = 0;
    do { // at the end of this fake loop we'll merge (if we get that far)
      
      if (fHitMergeChiCut <= 0. || hit.Multiplicity() != 2) break; // do not merge
      LOG_TRACE("ClusterCrawlerAlg")
        << "Hit #" << imbest << " is in a doublet; maybe we'll merge it.";
      
      // don't merge hits if near a vertex
      bool doMerge = true;
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(std::abs(kwire - vtx[ivx].Wire) < 10 &&
           std::abs(int(hit.PeakTime() - vtx[ivx].Time)) < 20 )
        {
          LOG_TRACE("ClusterCrawlerAlg")
            << "  hit #" << imbest << " not merged: too close to vtx #" << ivx;
          doMerge = false;
          break;
        }
      } // ivx
      if (!doMerge) break; // do not merge
      
      // find the neighbor hit
      size_t imbestn = HitInCluster_t::InvalidHitIndex;
      if(hit.LocalIndex() == 0) {
        imbestn = NextHitPresent(imbest);
        if (imbestn == HitInCluster_t::InvalidHitIndex) {
          mf::LogError("ClusterCrawlerAlg")
            << "Can't find the other hit in the doublet of hits of which #"
            << imbest << " is the first hit";
          break; // do not merge
        }
      } else {
        imbestn = PrevHitPresent(imbest);
        if (imbestn == HitInCluster_t::InvalidHitIndex) {
          mf::LogError("ClusterCrawlerAlg")
            << "Can't find the other hit in the doublet of hits of which #"
            << imbest << " is the last hit";
          break; // do not merge
        }
      }
      
      recob::Hit const& other_hit = fHits[imbestn];
      if (!areInSameMultiplet(hit, other_hit)) {
        mf::LogError("ClusterCrawlerAlg")
          << "Trying to merge hits from different multiplets: "
          << other_hit.WireID() << " @" << other_hit.StartTick()
            << " " << other_hit.LocalIndex() << "/" << other_hit.Multiplicity()
          << " and " << hit.WireID() << " @" << hit.StartTick()
            << " " << hit.LocalIndex() << "/" << hit.Multiplicity()
          ;
        break; // do not merge
      } // if not in the same multiplet
      
      // is the neighbor close and unused?
      // TODO use Hit::TimeDistanceAsRMS()
      float hitSep = std::abs(hit.PeakTime() - other_hit.PeakTime());
      hitSep = hitSep / hit.RMS();
      if (hitSep >= fHitMergeChiCut) {
        LOG_TRACE("ClusterCrawlerAlg")
          << "  hit #" << imbest << " not merged with #" << imbestn
          << ", far " << hitSep << " > " << fHitMergeChiCut << " in time";
        break; // do not merge
      }
      if (!isHitFree(imbestn)) {
        LOG_TRACE("ClusterCrawlerAlg")
          << "  hit #" << imbest << " not merged with #" << imbestn
          << ", that is not free";
        break; // do not merge
      }
      
      // there is a nearby hit
      hnear = 1;
      // Is the charge of the doublet more similar to the charge of the
      // previously added hits than the single hit
      float totChg = hit.Integral() + other_hit.Integral();
      float lastHitChg = fAveChg;
      if(lastHitChg < 0) lastHitChg = fHits[lastClHit].Integral();
      
      // decide whether to merge
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" merge hits charge check: totChg "<<totChg<<" lastHitChg "<<lastHitChg
    <<" hit chg "<<hit.Integral();
      if(std::abs(totChg - lastHitChg) >= std::abs(hit.Integral() - lastHitChg)) {
        LOG_TRACE("ClusterCrawlerAlg")
          << "  hit #" << imbest << " not merged with #" << imbestn
          << ", charge check failed (" << totChg << " " << lastHitChg
          << " " << hit.Integral() << ")";
        break; // do not merge
      }
      // the total charge of both hits is a better match than the 
      // charge of the hit selected
      MergeHits((unsigned short) imbest);
      // the nearby hit was merged
      hnear = -1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" Merging hits "<<imbest<<" and "<<imbestn
    <<" New Time "<<hit.PeakTime()
    <<" New Chg "<<hit.Integral();
      
    } while (false); // fake loop, runs only once

    // Make a charge similarity cut if the average charge is defined
    bool fitChg = true;
    if(fAveChg > 0.) {

      float chgrat = (hit.Integral() - fAveChg) / fAveChg;
    if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
      <<" Chgrat "<<std::setprecision(2)<<chgrat;

      // charge is way too high?
      if(chgrat > 2 * fChgCut[pass]) {
        if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" fails high charge cut";
        return;
      }

      // Determine if the last hit added was a large (low) charge hit
      // This will be used to prevent adding large (low) charge hits on two
      // consecutive fits. This cut is only applied to hits on adjacent wires
      float bigchgcut = 1.5 * fChgCut[pass];
      bool lasthitbig = false;
      bool lasthitlow = false;
      if(std::abs(wire0 - kwire) == 1) {
        float lastchgrat = (fHits[lastClHit].Integral() - fAveChg) / fAveChg;
        lasthitbig = ( lastchgrat > bigchgcut);
        lasthitlow = ( lastchgrat < -fChgCut[pass]);
      }
      
      // the last hit added was low charge and this one is as well
      if(lasthitlow && chgrat < -fChgCut[pass]) {
        if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" fails low charge cut. Stop crawling.";
        SigOK = false;
        return;
      } // lasthitlow
    
      // the last hit was high charge and this one is also
      if(lasthitbig && chgrat > fChgCut[pass]) {
        if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" fails 2nd high charge cut";
        return;
      } // lasthitbig

    
      // require that large charge hits have a very good projection error
      if(chgrat > fChgCut[pass]) {
        if(best > 2 * err) {
          if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" high charge && bad dT= "
            <<best<<" err= "<<err;
          return;
        }
      } // chgrat > fChgCut[pass]

      // decide whether to fit the charge
      fitChg = (chgrat < std::abs(fChgCut[pass]) );
    } // fAveChg > 0
    
    // we now have a hit that meets all the criteria. Fit it
    fcl2hits.push_back(imbest);
    FitCluster();
    chifits.push_back(clChisq);
    hitnear.push_back(hnear);
    // nearby hit check
    ChkClusterNearbyHits(prt);
    HitOK = true;

  if(prt) {
    mf::LogVerbatim("ClusterCrawlerAlg")
      <<" >>ADD W:T "<<kwire<<":"<<(short)hit.PeakTime()<<" dT "<<best
      <<std::setprecision(2)<<" Chisq "<<clChisq
      <<" Chg "<<(int)hit.Integral()<<" AveChg "<<(int)fAveChg
      <<" hitnear "<<hit.Multiplicity()
      <<" fcl2hits size "<<fcl2hits.size();
  }
    if(!fitChg) return;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Fit charge ";
    FitClusterChg();
  } // AddHit()


//////////////////////////////////////
    void ClusterCrawlerAlg::ChkClusterNearbyHits(bool prt)
    {
      // analyze the hitnear vector
      //  0 = no nearby hit exists
      //  1 = a nearby hit exists but was not merged
      // -1 = a nearby hit was merged
      
      if(fHitMergeChiCut <= 0) return;
      
      if(hitnear.size() != fcl2hits.size()) {
        mf::LogWarning("ClusterCrawlerAlg")<<"Coding error: hitnear size != fcl2hits";
        return;
      }
      
      // Analyze the last 6 hits added but don't consider the first few hits
      if(hitnear.size() < 12) return;
      
      // TODO move into loops
      unsigned short ii, indx;
      unsigned short merged = 0;
      unsigned short notmerged = 0;
      for(ii = 0; ii < 6; ++ii) {
        indx = hitnear.size() - 1 - ii;
        if(hitnear[indx] > 0) ++notmerged;
        if(hitnear[indx] < 0) ++merged;
      }
      
//  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
//    <<"ChkClusterNearbyHits: nearby hits merged "<<merged
//    <<" not merged "<<notmerged;

      if(notmerged < 2) return;
      
      // a number of nearby hits were not merged while crawling, so the 
      // average charge is probably wrong. Look at the last 6 hits added
      // and merge them if they are close
      for(ii = 0; ii < 6; ++ii) {
        indx = fcl2hits.size() - 1 - ii;
        const unsigned short iht = fcl2hits[indx];
        recob::Hit const& hit = fHits[iht];
        if(hit.Multiplicity() == 2) {
          // hit doublet. Get the index of the other hit
          unsigned short oht;
          if(hit.LocalIndex() == 0) {
            oht = NextHitPresent(iht);
          } else {
            oht = PrevHitPresent(iht);
          } // hit.LocalIndex() == 0
          recob::Hit const& other_hit = fHits[oht];
          // TODO use Hit::TimeDistanceAsRMS()
          float hitSep = std::abs(hit.PeakTime() - other_hit.PeakTime());
          hitSep /= hit.RMS();
          if(hitSep < fHitMergeChiCut && isHitFree(oht)) {
            // check charge assuming the hits are merged
            float chgRat = hit.Integral() + other_hit.Integral();
            chgRat /= hit.Integral();
            // merge em
            if(chgRat < 4) {
      if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Merging hit doublet "<<iht;
              MergeHits(iht);
              hitnear[indx] = -1;
            }
          } // hitSep OK and not in a cluster
        } // hit doublet
      } // ii
      
      // now re-fit
      FitCluster();
      FitClusterChg();

      if(prt) mf::LogVerbatim("ClusterCrawlerAlg")
        <<"ChkClusterNearbyHits refit cluster. fAveChg= "<<fAveChg;
      
    } // ChkClusterHitNear()

//////////////////////////////////////
    void ClusterCrawlerAlg::FitVtx(unsigned short iv, float& ChiDOF)
    {
      
      std::vector<float> x;
      std::vector<float> y;
      std::vector<float> ey2;
      
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].EndVtx == iv) {
          x.push_back(tcl[icl].EndSlp);
          float arg = tcl[icl].EndSlp * tcl[icl].EndWir - tcl[icl].EndTim;
          y.push_back(arg);
          if(tcl[icl].EndSlpErr > 0.) {
            arg = tcl[icl].EndSlpErr * tcl[icl].EndWir;
          } else {
            arg = .01 * tcl[icl].EndWir;
          }
          ey2.push_back(arg * arg);
        } else if(tcl[icl].BeginVtx == iv) {
          x.push_back(tcl[icl].BeginSlp);
          float arg = tcl[icl].BeginSlp * tcl[icl].BeginWir - tcl[icl].BeginTim;
          y.push_back(arg);
          if(tcl[icl].BeginSlpErr > 0.) {
            arg = tcl[icl].BeginSlpErr * tcl[icl].BeginWir;
          } else {
            arg = .01 * tcl[icl].BeginWir;
          }
          ey2.push_back(arg * arg);
        }
      } // ii
      if(x.size() < 2) {
        vtx[iv].Wght = -1;
        return;
      }
      
      float tv = 0.;
      float tverr = 0.;
      float wv = 0.;
      float wverr = 0.;
      fLinFitAlg.LinFit(x, y, ey2, tv, wv, tverr, wverr, ChiDOF);
      float vtime = -tv;
      float vwire = wv + 0.5;
      if(ChiDOF > 5) return;
      if(std::abs(vwire - vtx[iv].Wire) > 2) return;
      if(std::abs(vtime - vtx[iv].Time) > 10) return;
      vtx[iv].Wire = vwire;
      vtx[iv].Time = vtime;
    } // FitVtx

//////////////////////////////////////
    void ClusterCrawlerAlg::Vtx3ClusterMatch(geo::TPCID const& tpcid)
      {
        // Look for clusters that end/begin near the expected wire/time
        // for incomplete 3D vertices
        if(vtx3.size() == 0) return;
        
        const unsigned int cstat = tpcid.Cryostat;
        const unsigned int tpc = tpcid.TPC;
        
        unsigned short thePlane, theWire;
        float dw, theTime, dt;
        short dwb, dwe;

        for(unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
          // A complete 3D vertex with matching 2D vertices in all planes?
          if(vtx3[ivx].Wire < 0) continue;
	  if(vtx3[ivx].CStat != cstat || vtx3[ivx].TPC != tpc) continue;
          // Find the plane that is missing a 2D vertex
          thePlane = 3;
          theWire = vtx3[ivx].Wire;
          for(plane = 0; plane < 3; ++plane) {
            if(vtx3[ivx].Ptr2D[plane] >= 0) continue;
            thePlane = plane;
            break;
          } // plane
          if(thePlane > 2) continue;
          theTime = detprop->ConvertXToTicks(vtx3[ivx].X, thePlane, tpc, cstat);
          clCTP = EncodeCTP(cstat, tpc, thePlane);
          for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
            if(tcl[icl].ID < 0) continue;
            if(tcl[icl].CTP != clCTP) continue;
            dwb = std::abs(theWire - tcl[icl].BeginWir);
            dwe = std::abs(theWire - tcl[icl].EndWir);
            if(dwb < dwe && dwb < 5) {
              // cluster begin is closer
              if(theWire < tcl[icl].BeginWir) continue;
              dw = theWire - tcl[icl].BeginWir;
              dt = tcl[icl].BeginTim + tcl[icl].BeginSlp * dw - theTime;
              if(std::abs(dt) > 10) continue;
              // create a new 2D vertex
              VtxStore vnew;
              vnew.Wire = tcl[icl].BeginWir;
              vnew.Time = tcl[icl].BeginTim;
              vnew.Wght = 10;
              vnew.Topo = 4;
              vnew.CTP = clCTP;
              vtx.push_back(vnew);
              unsigned short ivnew = vtx.size() -1;
              vtx3[ivx].Ptr2D[thePlane] = ivnew;
              vtx3[ivx].Wire = -1;
              tcl[icl].BeginVtx = ivnew;
  if(vtxprt) mf::LogVerbatim("CC")<<"Vtx3ClusterMatch: Attach cluster index "<<icl
    <<" to new Begin vtx "<<ivnew<<" in plane "<<thePlane;
            } else if(dwe < 5) {
              // cluster end is closer
              if(theWire > tcl[icl].EndWir) continue;
              dw = theWire - tcl[icl].EndWir;
              dt = tcl[icl].EndTim + tcl[icl].EndSlp * dw - theTime;
              if(std::abs(dt) > 10) continue;
              // create a new 2D vertex
              VtxStore vnew;
              vnew.Wire = tcl[icl].EndWir;
              vnew.Time = tcl[icl].EndTim;
              vnew.Wght = 10;
              vnew.Topo = 1;
              vnew.CTP = clCTP;
              vtx.push_back(vnew);
              unsigned short ivnew = vtx.size() -1;
              vtx3[ivx].Ptr2D[thePlane] = ivnew;
              vtx3[ivx].Wire = -1;
              tcl[icl].EndVtx = ivnew;
  if(vtxprt) mf::LogVerbatim("CC")<<"VtxMatch: Attach cluster index "<<icl
    <<" to new End vtx "<<ivnew<<" in plane "<<thePlane;
            } // dwb/dwe check
          } // icl
        } // ivx
      } // Vtx3ClusterMatch

//////////////////////////////////////
    void ClusterCrawlerAlg::Vtx3ClusterSplit(geo::TPCID const& tpcid)
      {
        // Try to split clusters in a view in which there is no 2D vertex
        // assigned to a 3D vertex
        if(vtx3.size() == 0) return;
	const unsigned int cstat = tpcid.Cryostat;
        const unsigned int tpc = tpcid.TPC;
	
	vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);
        
        unsigned short lastplane = 5, kcl, kclID;
        float dth, theTime;
        unsigned short thePlane, theWire, plane;
        unsigned short loWire, hiWire;

        for(unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
	  if(vtx3[ivx].CStat != cstat || vtx3[ivx].TPC != tpc) continue;
	  // Complete 3D vertex with matching 2D vertices in all planes?
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Vtx3ClusterSplit ivx "<<ivx
    <<" Ptr2D "<<vtx3[ivx].Ptr2D[0]<<" "<<vtx3[ivx].Ptr2D[1]<<" "<<vtx3[ivx].Ptr2D[2]
        <<" wire "<<vtx3[ivx].Wire;
          if(vtx3[ivx].Wire < 0) continue;
          // find the plane that needs to be studied
          thePlane = 3;
          theWire = vtx3[ivx].Wire;
          for(plane = 0; plane < 3; ++plane) {
            if(vtx3[ivx].Ptr2D[plane] >= 0) continue;
            thePlane = plane;
            break;
          } // plane
          if(thePlane > 2) continue;
          theTime = detprop->ConvertXToTicks((double)vtx3[ivx].X, 
            (int)thePlane, (int)tpcid.TPC, (int)tpcid.Cryostat);
          // get the hit range if necessary
          if(thePlane != lastplane) {
            clCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, thePlane);
            GetHitRange(clCTP, WireHitRange, fFirstWire, fLastWire);
            lastplane = thePlane;
          }
          // make a list of clusters that have hits near this point on nearby wires xxx
          std::vector<short> clIDs;
          if(theWire > fFirstWire + 5) { loWire = theWire - 5; } else { loWire = fFirstWire; }
          if(theWire < fLastWire  - 5) { hiWire = theWire + 5; } else { hiWire = fLastWire; }
    if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"3DVtx "<<ivx
      <<" look for cluster hits near P:W:T "<<thePlane<<":"<<theWire<<":"<<(int)theTime
      <<" Wire range "<<loWire<<" to "<<hiWire;
          for(unsigned short wire = loWire; wire < hiWire; ++wire) {
            unsigned short index = wire - fFirstWire;
            // ignore dead wires or wires with no hits
            if(WireHitRange[index].first < 0) continue;
            unsigned short firsthit = WireHitRange[index].first;
            unsigned short lasthit = WireHitRange[index].second;
            for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
              // ignore obsolete hits
              if(fHits[khit].Integral() < 0) continue;
              // ignore un-assigned hits
              if(!isHitInCluster(khit)) continue;
              if((unsigned short)fHitInCluster[khit] > tcl.size() + 1) {
                mf::LogError("ClusterCrawlerAlg")<<"Invalid hit InClus. "<<khit
                  <<" "<<fHitInCluster[khit];
                continue;
              }
              // check an expanded time range
              if(theTime < fHits[khit].StartTick() - 10) continue;
              if(theTime > fHits[khit].EndTick() + 10) continue;
              kclID = fHitInCluster[khit];
              kcl = kclID - 1;
              // ignore obsolete clusters
              if(tcl[kcl].ID < 0) continue;
              
              // put the cluster in the list if it's not there already
    if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Bingo "<<ivx<<" plane "<<thePlane
      <<" wire "<<wire<<" hit "<<fHits[khit].WireID().Wire<<":"<<(int)fHits[khit].PeakTime()
        <<" inClus "<<fHitInCluster[khit]
        <<" P:W:T "<<thePlane<<":"<<tcl[kcl].BeginWir<<":"<<(int)tcl[kcl].BeginTim;
              if(std::find(clIDs.begin(), clIDs.end(), kclID) == clIDs.end()) {
                // ignore long straight clusters
                if(tcl[kcl].tclhits.size() > 100 ) {
                  dth = tcl[kcl].BeginAng - tcl[kcl].EndAng;
    if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Long straight check: nhits "
      <<tcl[kcl].tclhits.size()<<" dth "<<dth;
                  if(std::abs(dth) < 0.05) continue;
                } // tcl[kcl].tclhits.size() > 100
                clIDs.push_back(kclID);
              } // std::find
            } // khit
          } // wire
          if(clIDs.size() == 0) continue;
    if(vtxprt) {
      for(unsigned int ii = 0; ii < clIDs.size(); ++ii) {
        mf::LogVerbatim("ClusterCrawlerAlg")<<" clIDs "<<clIDs[ii];
      }
    }

          // do a local fit near the crossing point and make a tighter time cut
          unsigned short ii, icl, jj, iht;
          short nhitfit;
          bool didit;
          for(ii = 0; ii < clIDs.size(); ++ii) {
            icl = clIDs[ii] - 1;
            didit = false;
            for(jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
              iht = tcl[icl].tclhits[jj];
              if(fHits[iht].WireID().Wire < theWire) {
                nhitfit = 3;
                if(jj > 3) nhitfit = -3;
                FitClusterMid(icl, iht, nhitfit);
                dth = clpar[0] + (theWire - fHits[iht].WireID().Wire) * clpar[1] - theTime;
    if(vtxprt) {
      for(unsigned short kk = 0; kk < clIDs.size(); ++kk) {
        mf::LogVerbatim("ClusterCrawlerAlg")<<" cls "<<icl<<" dt "<<dth;
      }
    }
                if(std::abs(dth) > 5) clIDs[ii] = -1;
                didit = true;
                break;
              } // fHits[iht].WireID().Wire < theWire
              if(didit) break;
            } // jj
            if(didit) break;
          } // ii
    if(vtxprt) {
      mf::LogVerbatim("ClusterCrawlerAlg")<<"clIDs size after fit "<<clIDs.size();
      for(ii = 0; ii < clIDs.size(); ++ii) {
        mf::LogVerbatim("ClusterCrawlerAlg")<<" clIDs "<<clIDs[ii];
      }
    }

          // see if any candidates remain
          unsigned short nok = 0;
          for(ii = 0; ii < clIDs.size(); ++ii) if(clIDs[ii] >= 0) ++nok;
          if(nok == 0) continue;

          // make a new 2D vertex
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Make new 2D vtx in plane "<<thePlane
    <<" from 3D vtx "<<ivx;
          unsigned short nvcl = 0;
          VtxStore vnew;
          vnew.Wire = theWire;
          vnew.Time = theTime;
          vnew.Wght = 10;
          vnew.Topo = 1;
          vnew.CTP = clCTP;
          vtx.push_back(vnew);
          // update the 2D -> 3D vertex pointer
          unsigned short ivnew = vtx.size() - 1;
          vtx3[ivx].Ptr2D[thePlane] = ivnew;
          // either split or attach clusters to this vertex 
          for(ii = 0; ii < clIDs.size(); ++ii) {
            if(clIDs[ii] < 0) continue;
            icl = clIDs[ii] - 1;
            // find the split position
            // split the cluster. Find the split position
            unsigned short pos = 0;
            for(unsigned short jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
              iht = tcl[icl].tclhits[jj];
              if(fHits[iht].WireID().Wire < theWire) {
                pos = jj;
                break;
              }
            } // jj
            if(pos == 0) {
              // vertex is DS of the cluster Begin
              tcl[icl].BeginVtx = ivnew;
              ++nvcl;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Attach to Begin "<<icl;
            } else if(pos > tcl[icl].tclhits.size()) {
              // vertex is US of the cluster Eend
              tcl[icl].EndVtx = ivnew;
              ++nvcl;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"Attach to End "<<icl;
            } else {
              // vertex is in the middle of the cluster
              SplitCluster(icl, pos, ivnew);
              tcl[icl].ProcCode += 10000;
              tcl[tcl.size()-1].ProcCode += 10000;
              nvcl += 2;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Split cluster "<<clIDs[ii]<<" at pos "<<pos;
            } // pos check
          } // ii
          // try to attach other clusters
//          VertexCluster(ivnew);
          // Fit the vertex position
          float chisq = 0;
          FitVtx(ivnew, chisq);
          // mark the 3D vertex as complete
          vtx3[ivx].Wire = -1;
        } // ivx
        
      } // Vtx3ClusterSplit()


//////////////////////////////////////
    void ClusterCrawlerAlg::VtxMatch(geo::TPCID const& tpcid)
    {
      // Create 3D vertices from 2D vertices. 3D vertices that are matched
      // in all three planes have Ptr2D >= 0 for all planes
      
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      unsigned int nPln = TPC.Nplanes();
      if(nPln != 3) return;
      
      const unsigned int cstat = tpcid.Cryostat;
      const unsigned int tpc = tpcid.TPC;
      
      vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);
      
      // wire spacing in cm
      float wirePitch = geom->WirePitch(0, 1, 0, tpcid.TPC, tpcid.Cryostat);
            
      // create a vector of vertex indices in each plane
      std::vector<std::vector<unsigned short>> vIndex;
      std::vector<unsigned short> temp;
      for(unsigned short ipl = 0; ipl < 3; ++ipl) {
        temp.clear();
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].Wght < 0) continue;
          geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
	  if (iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
          unsigned int vpl = iplID.Plane;
          if(ipl == vpl) temp.push_back(ivx);
        }
        vIndex.push_back(temp);
      }
      temp.clear();
      
      // vector of 2D vertices -> 3D vertices.
      std::vector<short> vPtr;
      for(unsigned short ii = 0; ii < vtx.size(); ++ii) vPtr.push_back(-1);
      
      // temp vector of all 2D vertex matches
      std::vector<Vtx3Store> v3temp;
      
      double y = 0, z = 0;
      TVector3 WPos = {0, 0, 0};
      // i, j,k indicates 3 different wire planes
      unsigned int ipl = 0, ii = 0, ivx = 0, jpl = 0, jj = 0, jvx = 0;
      unsigned int kpl = 0, kk = 0, kvx = 0;
      float iX = 0, jX = 0, kX = 0;
      float iWire = 0, jWire = 0, kWire = 0;
      float xbest = 0, ybest = 0, zbest = 0;
      // compare vertices in each view
      for(ipl = 0; ipl < 2; ++ipl) {
        for(ii = 0; ii < vIndex[ipl].size(); ++ii) {
          ivx = vIndex[ipl][ii];
          // vertex has been matched already
          if(vPtr[ivx] >= 0) continue;
  
          iX = detprop->ConvertTicksToX((double)vtx[ivx].Time, (int)ipl, 
            (int)tpc, (int)cstat);
          iWire = vtx[ivx].Wire;
          geo::WireID iWireID(tpcid.Cryostat, tpcid.TPC, ipl, iWire);
/*
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"ipl "<<ipl<<" ivx "<<ivx<<" W:T "<<(int)vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time
    <<" iX "<<iX;
*/
          float best = fVertex3DCut;
          // temp array of 2D vertex indices in each plane
          std::array<short, 3> t2dIndex = {-1, -1, -1};
          std::array<short, 3> tmpIndex = {-1, -1, -1};
          for(jpl = ipl + 1; jpl < 3; ++jpl) {
            for(jj = 0; jj < vIndex[jpl].size(); ++jj) {
              jvx = vIndex[jpl][jj];
              if(vPtr[jvx] >= 0) continue;
              jX = detprop->ConvertTicksToX((double)vtx[jvx].Time, (int)jpl, 
                (int)tpc, (int)cstat);
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"jpl "<<jpl<<" jvx "<<jvx<<" W:T "<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time
    <<" jX "<<jX;
              if(std::abs(jX - iX) > fVertex3DCut) continue;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"2DMatchX "<<jX - iX;
              jWire = vtx[jvx].Wire;
            // IntersectionPoint() does not provide any feedback in case of errors:
            //  geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
              geo::WireID jWireID(tpcid.Cryostat, tpcid.TPC, jpl, jWire);
              geo::WireIDIntersection wiresIntersection;
              if (!geom->WireIDsIntersect(iWireID, jWireID, wiresIntersection)) {
                // the wires do not intersect (within the physical TPC)
                continue;
              }
              y = wiresIntersection.y;
              z = wiresIntersection.z;
	      if (!TPC.ContainsYZ(y,z)){
		//if(y < YLo || y > YHi || z < ZLo || z > ZHi) {
                // the previous check is supposed to do all the work
                LOG_WARNING("ClusterCrawlerAlg")
                  << "BUG: this check SHOULD HAVE been unnecessary!";
                continue;
              }
              WPos[1] = y;
              WPos[2] = z;
              // look for the matching vertex in the 3rd plane
              kpl = 3 - ipl - jpl;
              kX = 0.5 * (iX + jX);
              kWire = geom->NearestWire(WPos, kpl, tpc, cstat);
              // save this incomplete 3D vertex
              Vtx3Store v3d;
              tmpIndex[ipl] = ivx;
              tmpIndex[jpl] = jvx;
              tmpIndex[kpl] = -1;
              v3d.Ptr2D = tmpIndex;
              v3d.X = kX;
              v3d.Y = y;
              v3d.Z = z;
              v3d.Wire = kWire;
              v3d.CStat = cstat;
              v3d.TPC = tpc;
              v3temp.push_back(v3d);

  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Match ivx "<<ivx
    <<" P:W:T "<<ipl<<":"<<(int)vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time
    <<" jvx "<<jvx
    <<" P:W:T "<<jpl<<":"<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time;

              for(kk = 0; kk < vIndex[kpl].size(); ++kk) {
                kvx = vIndex[kpl][kk];
                if(vPtr[kvx] >= 0) continue;
                float kvxX = detprop->ConvertTicksToX((double)vtx[kvx].Time, 
                  (int)kpl, (int)tpc, (int)cstat);
                // Wire difference (cm)
                float dW = wirePitch * (vtx[kvx].Wire - kWire);
                // X difference (cm)
                float dX = (kvxX - kX);
                float dr = 0.5 * std::sqrt(dW * dW + dX * dX);
                if(dr < best) {
                  best = dr;
                  xbest = (kvxX + 2 * kX) / 3;
                  ybest = y;
                  zbest = z;
                  t2dIndex[ipl] = ivx;
                  t2dIndex[jpl] = jvx;
                  t2dIndex[kpl] = kvx;
                }

  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<" kvx "<<kvx<<" kpl "<<kpl
    <<" wire "<<(int)vtx[kvx].Wire<<" kTime "<<(int)vtx[kvx].Time
    <<" dr "<<dr;

              } // kk
            } // jj
          } // jpl
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")<<"3DMatch best "<<best;
          if(best < fVertex3DCut) {
            Vtx3Store v3d;
            v3d.Ptr2D = t2dIndex;
            v3d.Wire = -1;
            v3d.X = xbest;
            v3d.Y = ybest;
            v3d.Z = zbest;
            vtx3.push_back(v3d);
            for(unsigned short jj = 0; jj < 3; ++jj) 
              if(t2dIndex[jj] >= 0) vPtr[t2dIndex[jj]] = vtx3.size() - 1;

 if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"New 3D vtx "<<vtx3.size()
    <<" X "<<v3d.X<<" Y "<<v3d.Y<<" Z "<<v3d.Z
    <<" t2dIndex "<<t2dIndex[0]<<" "<<t2dIndex[1]<<" "<<t2dIndex[2];

          } // best < dRCut
        } // ii
      } // ipl
      
      // ignore vertices in the v3temp array that are part of a real 3-plane
      // 3D vertex
      unsigned short vsize = vtx3.size();
      for(unsigned short it = 0; it < v3temp.size(); ++it) {
        bool keepit = true;
        for(unsigned short i3d = 0; i3d < vsize; ++i3d) {
          for(unsigned short plane = 0; plane < 3; ++plane) {
            if(v3temp[it].Ptr2D[plane] == vtx3[i3d].Ptr2D[plane]) {
              keepit = false;
              break;
            }
          } // plane
          if(!keepit) break;
        } // i3d

        if(keepit) vtx3.push_back(v3temp[it]);
      } // it
      

  if(vtxprt) {
    for(unsigned short it = 0; it < vtx3.size(); ++it) {
      mf::LogVerbatim("ClusterCrawlerAlg")
        <<"vtx3 "<<it<<" Ptr2D "<<vtx3[it].Ptr2D[0]<<
        " "<<vtx3[it].Ptr2D[1]<<" "<<vtx3[it].Ptr2D[2]
        <<" wire "<<vtx3[it].Wire;
    }
  }

    } // VtxMatch

//////////////////////////////////
    void ClusterCrawlerAlg::GetHitRange(
      CTP_t CTP, 
      std::vector< std::pair<short, short> >& WireHitRange,
      unsigned short& firstwire, unsigned short& lastwire)
    {
      // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
      firstwire = lastwire = 0;
      WireHitRange.clear();
      bool first = true;
      unsigned short firsthit = 0;
      geo::PlaneID planeID = DecodeCTP(CTP);
      unsigned short lasthit = 0;
      // find the first and last wire with a hit
      for(unsigned short hit = 0; hit < fHits.size(); ++hit) {
        // TODO change to planeID()
        if(fHits[hit].WireID().Plane != planeID.Plane) continue;
        if(fHits[hit].WireID().TPC != planeID.TPC) continue;
        if(fHits[hit].WireID().Cryostat != planeID.Cryostat) continue;
        unsigned short theWireNum = fHits[hit].WireID().Wire;
        if(first) {
          firsthit = hit;
          firstwire = theWireNum;
          first = false;
        }
        lastwire = theWireNum;
        lasthit = hit;
      } //hit
      
      if (first) return; // we collected nothing

      // now we can define the WireHitRange vector.
      // start by defining the "no hits on wire" condition
      short sflag = -2;
      for(unsigned short wire = firstwire; wire <= lastwire; ++wire) {
        WireHitRange.push_back(std::make_pair(sflag, sflag));
      }
      // overwrite with the "dead wires" condition
      filter::ChannelFilter cf;
      sflag = -1;
      for(unsigned short wire = firstwire+1; wire < lastwire; ++wire) {
        raw::ChannelID_t chan = geom->PlaneWireToChannel
          ((int)planeID.Plane,(int)wire,(int)planeID.TPC,(int)planeID.Cryostat);
        // remember to offset references to WireHitRange by the FirstWire
        unsigned short index = wire - firstwire;
        if(cf.BadChannel(chan)) WireHitRange[index] = std::make_pair(sflag, sflag);
      }
          
      lastwire = firstwire;
      unsigned short thishit = firsthit;
      unsigned short lastfirsthit = firsthit;
      // next overwrite with the index of the first/last hit on each wire
      for(unsigned short hit = firsthit; hit <= lasthit; ++hit) {
        recob::Hit const& theHit = fHits[hit];
        if(theHit.WireID().Plane != planeID.Plane) continue;
        if(theHit.WireID().TPC != planeID.TPC) continue;
        if(theHit.WireID().Cryostat != planeID.Cryostat) continue;
	
        unsigned short thiswire = theHit.WireID().Wire;
        if(thiswire > lastwire) {
          unsigned short index = lastwire - firstwire;
          short itmp1 = lastfirsthit;
          short itmp2 = thishit;
          WireHitRange[index] = std::make_pair(itmp1,itmp2);
          lastwire = thiswire;
          lastfirsthit = thishit;
        } else if(thiswire < lastwire) {
          mf::LogError("ClusterCrawlerAlg")<<"ERROR: Hits not sorted!!";
          return;
        }
        ++thishit;
      } //hit
      // define for the last wire
      unsigned short index = lastwire - firstwire;
      short itmp1 = lastfirsthit;
      short itmp2 = thishit;
      WireHitRange[index] = std::make_pair(itmp1,itmp2);
    } // GetHitRange()


//////////////////////////////////////////
    bool ClusterCrawlerAlg::areInSameMultiplet
      (recob::Hit const& first_hit, recob::Hit const& second_hit)
    {
      return (first_hit.StartTick() == second_hit.StartTick())
        && (first_hit.WireID() == second_hit.WireID())
    //   && (first_hit.StopTick() == second_hit.StopTick())
        ;
    } // ClusterCrawlerAlg::areInSameMultiplet()
    
    
    std::pair<size_t, size_t> ClusterCrawlerAlg::FindHitMultiplet
      (size_t iHit) const
    {
      std::pair<size_t, size_t> range{ iHit, iHit };
      
      const raw::TDCtick_t MultipletStartTime = fHits[iHit].StartTick();
      geo::WireID const& wid = fHits[iHit].WireID();
      
      // detect the first hit with the same start time as the reference hit
      while (range.first > 0) {
        if (fHits[range.first - 1].StartTick() != MultipletStartTime) break;
        if (fHits[range.first - 1].WireID() != wid) break;
        --range.first;
      } // while
      
      // detect the last hit with the same start time as the reference hit
      while (++range.second < fHits.size()) {
        if (fHits[range.second].StartTick() != MultipletStartTime) break;
        if (fHits[range.second].WireID() != wid) break;
      } // while
      
      return range;
    } // ClusterCrawlerAlg::FindHitMultiplet()
    
    
/////////////////////////////////////////
    void ClusterCrawlerAlg::SortByLength(
      std::vector<ClusterStore> const& tcl,
      CTP_t inCTP, std::map<unsigned short, unsigned short>& sortindex
      )
    {
      // sorts the temporary cluster vector by decreasing number of hits,
      // while ignoring abandoned clusters. Returns index map with the
      // sort order
      
      // form a vector of pairs of the number of hits and the index
      std::vector< std::pair<unsigned short, unsigned short> > index;
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        if(tcl[ii].ID > 0 && tcl[ii].CTP == inCTP) 
          index.push_back(std::make_pair(tcl[ii].tclhits.size(),ii));
      }
      std::sort(index.begin(), index.end(), SortByLen);
      sortindex.clear();
      for(unsigned short ii = 0; ii < index.size(); ++ii) {
       sortindex[ii]=index[ii].second;
      }
      return; 
    }
/////////////////////////////////////////
    
    size_t ClusterCrawlerAlg::HitInCluster_t::NextPresent
      (size_t iHit, size_t n /* = 1 */) const
    {
      // if we don't actually have to move...
      if (n == 0) return isPresent(iHit)? iHit: InvalidHitIndex;
      
      // start counting...
      while (++iHit < nHits()) {
        if (!isPresent(iHit)) continue;
        if (--n == 0) return iHit;
      } // while
      return InvalidHitIndex;
    } // ClusterCrawlerAlg::HitInCluster_t::NextPresent()
    
    
    size_t ClusterCrawlerAlg::HitInCluster_t::PrevPresent
      (size_t iHit, size_t n /* = 1 */) const
    {
      // if we don't actually have to move...
      if (n == 0) return isPresent(iHit)? iHit: InvalidHitIndex;
      
      // start counting...
      while (iHit-- > 0) {
        if (!isPresent(iHit)) continue;
        if (--n == 0) return iHit;
      } // while
      return InvalidHitIndex;
    } // ClusterCrawlerAlg::HitInCluster_t::PrevPresent()
    
    
    size_t ClusterCrawlerAlg::HitInCluster_t::NextFree
      (size_t iHit, size_t n /* = 1 */) const
    {
      // if we don't actually have to move...
      if (n == 0) return isFree(iHit)? iHit: InvalidHitIndex;
      
      // start counting...
      while (++iHit < nHits()) {
        if (!isFree(iHit)) continue;
        if (--n == 0) return iHit;
      } // while
      return InvalidHitIndex;
    } // ClusterCrawlerAlg::HitInCluster_t::NextFree()
    
    
    size_t ClusterCrawlerAlg::HitInCluster_t::PrevFree
      (size_t iHit, size_t n /* = 1 */) const
    {
      // if we don't actually have to move...
      if (n == 0) return isFree(iHit)? iHit: InvalidHitIndex;
      
      // start counting...
      while (iHit-- > 0) {
        if (!isFree(iHit)) continue;
        if (--n == 0) return iHit;
      } // while
      return InvalidHitIndex;
    } // ClusterCrawlerAlg::HitInCluster_t::PrevFree()
    
    
} // namespace cluster
