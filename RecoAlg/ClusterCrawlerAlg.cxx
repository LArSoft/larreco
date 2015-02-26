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

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include <iostream>
#include <iomanip>

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"

#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoAlg/ClusterCrawlerAlg.h"

namespace cluster {
//------------------------------------------------------------------------------
  ClusterCrawlerAlg::ClusterCrawlerAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    CrawlInit();
  }

//------------------------------------------------------------------------------
  ClusterCrawlerAlg::~ClusterCrawlerAlg()
  {
  }

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

    if(badinput) throw cet::exception("ClusterCrawler")<<"Bad input from fcl file ";

  } // reconfigure

  // used for sorting hits on wires
  bool SortByLowHit(short i, short j) {return ((i > j));}

  // used for sorting clusters by length
  typedef std::pair<unsigned int, unsigned int> mypair;
  bool SortByLen(const mypair& L, const mypair& R) {return (L.first > R.first);}


  void ClusterCrawlerAlg::CrawlInit() {
    prt = false; vtxprt = false;
    NClusters = 0;  clBeginSlp = 0; clBeginSlpErr = 0; clBeginTim = 0;
    clBeginWir = 0; clBeginChg = 0; clEndSlp = 0;      clEndSlpErr = 0;
    clEndTim = 0;   clEndWir = 0;   clEndChg = 0;      clChisq = 0;
    clStopCode = 0; clProcCode = 0; fFirstWire = 0;
    fLastWire = 0; fAveChg = 0.; fChgSlp = 0.; pass = 0;
    fScaleF = 0; tcl.clear(); vtx.clear(); vtx3.clear(); WireHitRange.clear();
  }


  void ClusterCrawlerAlg::RunCrawler(std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // Run the ClusterCrawler algorithm - creating seed clusters and crawling upstream.

    CrawlInit();
    
    // don't make clusters, just hits
    if(fNumPass == 0) return;
    
    if(allhits.size() < 3) return;
    if(allhits.size() > USHRT_MAX) {
      mf::LogWarning("CC")<<"Too many hits for ClusterCrawler "<<allhits.size();
      return;
    }
    
//    std::cout<<"CC: number of hits "<<allhits.size()<<"\n";
    
    for(cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      for(tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
        for(plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
          WireHitRange.clear();
          // define a code to ensure clusters are compared within the same plane
          clCTP = EncodeCTP(cstat, tpc, plane);
          // fill the WireHitRange vector with first/last hit on each wire
          // dead wires and wires with no hits are flagged < 0
          GetHitRange(allhits, clCTP, WireHitRange, fFirstWire, fLastWire);

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
      if(allhits[hit].WirID.Wire != wire) {
        std::cout<<"Bad wire "<<hit<<" "<<allhits[hit].WirID.Wire<<" "<<wire<<"\n";
        return;
      } // check wire
      if(allhits[hit].WirID.Plane != plane) {
        std::cout<<"Bad plane "<<hit<<" "<<allhits[hit].WirID.Plane<<" "<<plane<<"\n";
        return;
      } // check plane
      if(allhits[hit].WirID.TPC != tpc) {
        std::cout<<"Bad tpc "<<hit<<" "<<allhits[hit].WirID.TPC<<" "<<tpc<<"\n";
        return;
      } // check tpc
    } // hit
  } // wire
  std::cout<<" is OK. nhits "<<nhts<<"\n";;
*/
          fFirstHit = WireHitRange[0].first;
          uint32_t channel = allhits[fFirstHit].Wire->Channel();
          // get the scale factor to convert dTick/dWire to dX/dU. This is used
          // to make the kink and merging cuts
          float wirePitch = geom->WirePitch(geom->View(channel));
          float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
          tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
          fScaleF = tickToDist / wirePitch;
          // convert Large Angle Cluster crawling cut to a slope cut
          if(fLAClusAngleCut > 0) 
            fLAClusSlopeCut = tan(3.142 * fLAClusAngleCut / 180.) / fScaleF;
          fMaxTime = detprop->NumberTimeSamples();
          fNumWires = geom->Nwires(plane);
          
          // look for clusters
          ClusterLoop(allhits, tcl, vtx);
        } // plane
        if(fVertex3DCut > 0) {
          // Match vertices in 3 planes
          VtxMatch(allhits, tcl, vtx, vtx3, cstat, tpc);
          Vtx3ClusterMatch(allhits, tcl, vtx, vtx3, cstat, tpc);
          // split clusters using 3D vertices
          Vtx3ClusterSplit(allhits, tcl, vtx, vtx3, cstat, tpc);
        }
        if(vtxprt) {
          mf::LogVerbatim("ClusterCrawler")<<"Clustering done in TPC ";
          PrintClusters(allhits, tcl, vtx);
        }
      } // tpc
    } // cstat
/*
  cstat = 0; tpc = 0;
  for(plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
    WireHitRange.clear();
    // define a code to ensure clusters are compared within the same plane
    clCTP = EncodeCTP(cstat, tpc, plane);
    // fill the WireHitRange vector with first/last hit on each wire
    // dead wires and wires with no hits are flagged < 0
    GetHitRange(allhits, clCTP, WireHitRange, fFirstWire, fLastWire);
    unsigned int nhts = 0, nhtsinCls = 0;
    for(unsigned short wire = fFirstWire; wire < fLastWire; ++wire) {
      unsigned short index = wire - fFirstWire;
      unsigned short fhit = WireHitRange[index].first;
      unsigned short lhit = WireHitRange[index].second;
      for(unsigned short hit = fhit; hit < lhit; ++hit) {
        if(allhits[hit].InClus >= 0) ++nhts;
        if(allhits[hit].InClus > 0) ++nhtsinCls;
      } // hit
    } // wire
    std::cout<<"After CC: plane "<<plane<<" nhits "<<nhts
      <<" nhits in Clusters "<<nhtsinCls<<"\n";;
  } // plane
*/     
    WireHitRange.clear(); 
    fcl2hits.clear();
    chifits.clear();
    hitnear.clear();
    
  } // RunCrawler
    
////////////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterLoop(std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx) {
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
            // skip used hits
            if(ihit > allhits.size()-1) {
              mf::LogError("ClusterCrawler")<<"RunCrawler bad ihit "<<ihit;
              return;
            }
            // skip used and obsolete hits
            if(allhits[ihit].InClus != 0) continue;
            // Check for a hit signal on the next DS wire
            bool SigOK = false;
            ChkSignal(allhits, iwire + 1, allhits[ihit].Time - 10,
                               iwire + 1, allhits[ihit].Time + 10, SigOK);
            // Don't start a seed cluster if there is a hit signal DS. 
            // This is an indicator that we might be trying
            // to start a cluster just US of shower blob.
            if(SigOK && allhits[ihit].numHits > 1) continue;
            if((iwire - span + 1) < fFirstWire) continue;
            unsigned short jwire = iwire - span + 1;
            unsigned short jindx = jwire - fFirstWire;
            if(WireHitRange[jindx].first < 0) continue;
            // Find the hit on wire jwire that best matches a line between
            // a nearby vertex and hit ihit. No constraint if useHit < 0
            unsigned short useHit = 0;
            bool doConstrain = false;
            VtxConstraint(allhits, vtx, iwire, ihit, jwire, useHit, doConstrain);
            unsigned short jfirsthit = WireHitRange[jindx].first;
            unsigned short jlasthit = WireHitRange[jindx].second;
            for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
              if(jhit > allhits.size()-1) {
                mf::LogError("ClusterCrawler")<<"RunCrawler bad jhit "<<jhit;
                return;
              }
              // skip used and obsolete hits
              if(allhits[jhit].InClus != 0) continue;
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
              clpar[0] = allhits[jhit].Time;
              clpar[1] = (allhits[ihit].Time - allhits[jhit].Time) / (iwire - jwire);
              // check for hit width consistency. Large angle clusters should have
              // large rms hits or smaller rms hits if they part of a multiplet
              if(fabs(clpar[1]) > 5) {
                if(pass == 0) continue;
                if(allhits[ihit].RMS < 4 && allhits[ihit].numHits < 3) continue;
                if(allhits[jhit].RMS < 4 && allhits[jhit].numHits < 3) continue;
              } // fabs(clpar[1]) > 5
              clChisq = 0;
              // now look for hits to add on the intervening wires
              SigOK = false;
              bool HitOK = false;
              bool clok = true;
              for(unsigned short kwire = jwire+1; kwire < iwire; ++kwire) {
                // ensure this cluster doesn't cross a vertex
                if(CrawlVtxChk(allhits, vtx, kwire)) {
                  clok = false;
                  break;
                }
                AddHit(allhits, kwire, HitOK, SigOK);
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
              FitCluster(allhits);
              if(clChisq > 10.) continue;
              // check the charge ratio between the DS hit and the next-most
              // DS hit. This ratio is < 2 for a stopping particle. A large
              // ratio indicates that we are starting a cluster just US of a
              // high ionization region
              float chg0 = allhits[fcl2hits[0]].Charge;
              float chg1 = allhits[fcl2hits[1]].Charge;
              if(chg0 > 2 * chg1) continue;
              // check the hit RMS ratio
              chg0 = allhits[fcl2hits[0]].RMS;
              chg1 = allhits[fcl2hits[1]].RMS;
              float rmsrat = fabs(chg0 - chg1) / (chg0 + chg1);
              if(rmsrat > 0.3) continue;
              // save the cluster begin info
              clBeginWir = iwire;
              clBeginTim = allhits[ihit].Time;
              clBeginSlp = clpar[1];
              // don't do a small angle crawl if the cluster slope is too large
              // and Large Angle crawling is NOT requested on this pass
              if(!fLACrawl[pass] && fabs(clBeginSlp) > fLAClusSlopeCut) continue;
              // See if we are trying to start a cluster between a vertex
              // and a cluster that is associated to that vertex. If so, skip it
              if(CrawlVtxChk2(allhits, tcl, vtx)) continue;
              clBeginSlpErr = clparerr[1];
              clBeginChg = 0;
              // decide whether to crawl a large angle cluster. Requirements are:
              // 1) the user has set the LACluster angle cut > 0, AND
              // 2) the cluster slope exceeds the cut
              // Note that if condition 1 is met, normal cluster crawling is done
              // only if the slope is less than the cut
              if(fLACrawl[pass] && fLAClusSlopeCut > 0) {
                // LA cluster crawling requested
                if(fabs(clBeginSlp) > fLAClusSlopeCut) {
                  LACrawlUS(allhits, vtx);
                } else {
                  CrawlUS(allhits, vtx);
                } // fabs(clBeginSlp) > fLAClusAngleCut
              } else {
                // allow clusters of any angle
                CrawlUS(allhits, vtx);
              } // fLAClusSlopeCut > 0
              if(fcl2hits.size() >= fMinHits[pass]) {
                // it's long enough so save it
                clEndSlp = clpar[1]; // save the slope at the end
                clEndSlpErr = clparerr[1];
                TmpStore(allhits, tcl); // store the cluster
                ClusterAdded = true;
                nHitsUsed += fcl2hits.size();
                AllDone = (nHitsUsed == allhits.size());
                break;
              }
              // kill it
            } // jhit
            if(ClusterAdded || AllDone) break;
          } // ihit
          if(AllDone) break;
        } // iwire

        // try to merge clusters 
        if(fDoMerge[pass]) ChkMerge(allhits, tcl);
        // form 2D vertices
        if(fFindVertices[pass]) FindVertices(allhits, tcl, vtx);

        if(AllDone) break;
        
      } // pass

      // Kill vertices that are close to long straight tracks. These
      // spurious vertices are likely due to delta rays on muon tracks
//      KillVertices(allhits, tcl, vtx);
      // Merge overlapping clusters
      if(fMergeOverlap) MergeOverlap(allhits, tcl, vtx);
      // Check the DS end of clusters
      if(fChkClusterDS) ChkClusterDS(allhits, tcl);
      // split clusters using vertices
      if(fVtxClusterSplit) VtxClusterSplit(allhits, tcl, vtx);
      // Look for 2D vertices with star topology - short, back-to-back clusters
      if(fFindStarVertices) FindStarVertices(allhits, tcl, vtx);

      if(fDebugPlane == (short)plane) {
        mf::LogVerbatim("ClusterCrawler")<<"Clustering done in plane "<<plane;
        PrintClusters(allhits, tcl, vtx);
      }
    
  } // ClusterLoop

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeOverlap(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {
      unsigned short icl, jcl, jj, jht, tclsize;
      unsigned short imidwir, jmidwir;
      float ithb, jthe, da, dw, dt;
      
      prt = (fDebugWire == 666);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"MergeOverlap check. clCTP "<<clCTP;

      tclsize = tcl.size();
      for(icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].BeginVtx >= 0) continue;
        if(tcl[icl].tclhits.size() < 10) continue;
        ithb = atan(fScaleF * tcl[icl].BeginSlp);
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
          jthe = atan(fScaleF * tcl[jcl].EndSlp);
          if(fabs(ithb - jthe) > 0.2) continue;
          // fit hits on jcl near the Begin wire of icl to make
          // tighter cuts. Find the hits on jcl
          for(jj = 0; jj < tcl[jcl].tclhits.size(); ++jj) {
            jht = tcl[jcl].tclhits[jj];
            if(allhits[jht].WireNum <= tcl[icl].BeginWir) break;
          } // jj
          FitClusterMid(allhits, tcl, jcl, jht, 5);
          // compare the angle
          da = atan(fScaleF * clpar[1]) - ithb;
          dw = tcl[icl].BeginWir - allhits[jht].WireNum;
          dt = clpar[0] + dw * clpar[1] - tcl[icl].BeginTim;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" overlap? "<<icl
    <<" and "<<jcl<<" da "<<da<<" dt "<<dt;
          if(fabs(da) > 0.1) continue;
          if(fabs(dt) > 10) continue;
          // merge clusters
          DoMerge(allhits, tcl, vtx, icl, jcl, 600);
        } // jcl
      } // icl
      
    } // MergeOverlap

//////////////////////////////////////////
    void ClusterCrawlerAlg::ChkClusterDS(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl)
    {
      // Try to extend clusters DS by a few wires. 
      // DS hits may not have been  included in a cluster if they have high 
      // multiplicity or high charge. 
      // Ref ClusterLoop cuts for starting a seed cluster.

  prt = (fDebugPlane == 3);

  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS clCTP "<<clCTP;

      unsigned short icl, tclsize, ih0, ih1, iht, wire;
      unsigned short ii, newcl;
      float slp, thhits, prevth, hitrms, rmsrat;
      bool ratOK;
      std::vector<unsigned short> dshits;
      tclsize = tcl.size();
      for(icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore clusters that have a Begin vertex
        if(tcl[icl].BeginVtx >= 0) continue;
        // find the angle using the first 2 hits
        ih0 = tcl[icl].tclhits[1];
        ih1 = tcl[icl].tclhits[0];
        slp = (allhits[ih1].Time    - allhits[ih0].Time) / 
              (allhits[ih1].WireNum - allhits[ih0].WireNum);
        prevth = atan(fScaleF * slp);
        // move the "origin" to the first hit
        ih0 = ih1;
        wire = allhits[ih0].WireNum;
        hitrms = allhits[ih0].RMS;
        dshits.clear();
        // follow DS a few wires. Stop and do nothing if any encountered 
        // hit is associated with a cluster
        for(ii = 0; ii < 2; ++ii) {
          ++wire;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS: Extend "
    <<tcl[icl].ID<<" to W:T "<<wire;
          if(wire > fLastWire) break;
          unsigned short index = wire - fFirstWire;
          // stop if no hits on this wire
          if(WireHitRange[index].first == -2) break;
          unsigned short firsthit = WireHitRange[index].first;
          unsigned short lasthit = WireHitRange[index].second;
          bool hitAdded = false;
          for(ih1 = firsthit; ih1 < lasthit; ++ih1) {
            if(allhits[ih1].InClus != 0) continue;
            slp = (allhits[ih1].Time - allhits[ih0].Time) /
                  (allhits[ih1].WireNum - allhits[ih0].WireNum);
            thhits = atan(fScaleF * slp);
            rmsrat = allhits[ih1].RMS / hitrms;
            ratOK = rmsrat > 0.3 && rmsrat < 3;
            // require small angle and not wildly different width compared
            // to the first hit in the cluster
            if(ratOK && fabs(thhits - prevth) < 0.3) {
              dshits.push_back(ih1);
              hitAdded = true;
              prevth = thhits;
              ih0 = ih1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Add hit "<<allhits[ih1].WireNum
  <<":"<<(int)allhits[ih1].Time<<" rmsrat "<<rmsrat;
              break;
            }
          } // ih1
          // stop looking if no hit was added on this wire
          if(!hitAdded) break;
        } // ii
        // Found hits not associated with a different cluster
        if(dshits.size() > 0) {
          // put the tcl cluster into the working vectors
          TmpGet(tcl, icl);
          // clobber the hits
          fcl2hits.clear();
          // sort the DS hits
          if(dshits.size() > 1) std::sort(dshits.begin(), dshits.end(), SortByLowHit);
          // stuff them into fcl2hits
          fcl2hits = dshits;
          // Append the existing hits
          for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            // un-assign the hits so that TmpStore will re-assign them
            iht = tcl[icl].tclhits[ii];
            allhits[iht].InClus = 0;
            fcl2hits.push_back(iht);
          }
          clProcCode += 5000;
          pass = fNumPass - 1;
          FitClusterChg(allhits);
          clBeginChg = fAveChg;
          TmpStore(allhits, tcl);
          newcl = tcl.size() -1;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Store "<<newcl;
          tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
          tcl[newcl].EndVtx = tcl[icl].EndVtx;
          // declare the old one obsolete
          tcl[icl].ID = -tcl[icl].ID;
        } // dshits.size() > 0
      } // icl
      
      // look for isolated stopping clusters. Check for hit multiplets
      // on the DS end and merge them
      unsigned short nhts, nmult, loWire, hiWire, index, jht, ivx;
      unsigned short nUsInClus, nUsNotInClus;
      unsigned short nDsInClus, nDsNotInClus;
      float prTimeLo, prTimeHi;
      for(icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore clusters that have a Begin vertex
        if(tcl[icl].BeginVtx >= 0) continue;
        // ignore short clusters
        if(tcl[icl].tclhits.size() < 6) continue;
        nhts = 4;
        nmult = 0;
        loWire = 9999;
        for(ii = 0; ii < nhts; ++ii) {
          iht = tcl[icl].tclhits[ii];
          // ignore hits close to a vertex
          if(tcl[icl].EndVtx >= 0) {
            ivx = tcl[icl].EndVtx - 1;
            if((allhits[iht].WireNum - vtx[ivx].Wire) < 5) continue;
          } // tcl[icl].EndVtx >= 0
          if(allhits[iht].numHits > 1) ++nmult;
          if(allhits[iht].WireNum < loWire) loWire = allhits[iht].WireNum;
          if(allhits[iht].WireNum > hiWire) hiWire = allhits[iht].WireNum;
        } // ii
        if(nmult == 0) continue;
        // count the number of hits in this region that are not associated
        // with this cluster
        // inspect the 5 wires US of the cluster end and the 5 wires DS
        hiWire = tcl[icl].BeginWir + 5;
        nDsInClus = 0; nDsNotInClus = 0;
        nUsInClus = 0; nUsNotInClus = 0;
        for(wire = loWire; wire < hiWire; ++wire) {
          index = wire - fFirstWire;
          // ignore dead wires and those with no hits
          if(WireHitRange[index].first < 0) continue;
          prTimeLo = tcl[icl].BeginTim + tcl[icl].BeginSlp * (wire - tcl[icl].BeginWir);
          prTimeHi = prTimeLo + 30;
          prTimeLo -= 30;
          for(jht = WireHitRange[index].first; jht < WireHitRange[index].second; ++jht) {
            // ignore obsolete hits
            if(allhits[jht].InClus < 0) continue;
            if(allhits[jht].Time > prTimeHi) continue;
            if(allhits[jht].Time < prTimeLo) continue;
            // ignore hits associated with this cluster
            if(allhits[jht].InClus == tcl[icl].ID) continue;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" not "
    <<allhits[jht].WireNum<<":"<<(int)allhits[jht].Time
    <<" InClus "<<allhits[jht].InClus;
            if(wire <= tcl[icl].BeginWir) {
              if(allhits[jht].InClus > 0) {
                ++nUsInClus;
              } else {
                ++nUsNotInClus;
              }
            } else {
              if(allhits[jht].InClus > 0) {
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
          for(ii = 0; ii < nhts; ++ii) {
            iht = tcl[icl].tclhits[ii];
            // ignore hits close to a vertex
            if(tcl[icl].EndVtx >= 0) {
              ivx = tcl[icl].EndVtx - 1;
              if((allhits[iht].WireNum - vtx[ivx].Wire) < 5) continue;
            } // tcl[icl].EndVtx >= 0
            if(allhits[iht].numHits > 1) MergeHits(allhits, iht);
          } // ii
        } // nmult > 0
        if(nDsNotInClus > 0) {
          // DS hits found that are not included. Try to attach
          ih0 = tcl[icl].tclhits[1];
          ih1 = tcl[icl].tclhits[0];
          slp = (allhits[ih1].Time    - allhits[ih0].Time) /
                (allhits[ih1].WireNum - allhits[ih0].WireNum);
          prevth = atan(fScaleF * slp);
          // move the "origin" to the first hit
          ih0 = ih1;
          wire = allhits[ih0].WireNum;
          hitrms = allhits[ih0].RMS;
          dshits.clear();
          for(ii = 0; ii < 2; ++ii) {
            ++wire;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<"ChkClusterDS: Extend "
    <<tcl[icl].ID<<" to W:T "<<wire;
            if(wire > fLastWire) break;
            index = wire - fFirstWire;
            if(WireHitRange[index].first == -2) break;
            bool hitAdded = false;
            for(ih1 = WireHitRange[index].first; ih1 < WireHitRange[index].second; ++ih1) {
              if(allhits[ih1].InClus != 0) continue;
              slp = (allhits[ih1].Time - allhits[ih0].Time) /
                    (allhits[ih1].WireNum - allhits[ih0].WireNum);
              thhits = atan(fScaleF * slp);
              if(fabs(thhits - prevth) > 0.4) continue;
              rmsrat = allhits[ih1].RMS / hitrms;
  if(prt) mf::LogVerbatim("ClusterCrawlerAlg")<<" Found hit. Chk last hit RMS "
    <<hitrms<<" w W:T "<<allhits[ih1].WireNum<<":"<<(int)allhits[ih1].Time
    <<" RMS "<<allhits[ih1].RMS<<" rmsrat "<<rmsrat;
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
            TmpGet(tcl, icl);
            // clobber the hits
            fcl2hits.clear();
            // sort the DS hits
            if(dshits.size() > 1) std::sort(dshits.begin(), dshits.end(), SortByLowHit); 
            fcl2hits = dshits;
            for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
              // un-assign the hits so that TmpStore will re-assign them
              iht = tcl[icl].tclhits[ii];
              allhits[iht].InClus = 0;
              fcl2hits.push_back(iht);
            } // ii
            clProcCode += 5000;
            pass = fNumPass - 1;
            FitClusterChg(allhits);
            clBeginChg = fAveChg;
            TmpStore(allhits, tcl);
            newcl = tcl.size() -1;
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
  bool ClusterCrawlerAlg::CrawlVtxChk2(
    std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
  {
    // returns true if the (presumably short) cluster under construction
    // resides between a vertex and another cluster that is associated with
    // that vertex
    
    if(vtx.size() == 0) return false;
    if(fcl2hits.size() == 0) return false;
    
    unsigned short iht = fcl2hits.size() - 1;
    unsigned short icl;
    float wire0 = (0.5 + allhits[fcl2hits[iht]].WireNum);
    float dt;
    float thclus = atan(fScaleF * clpar[1]);

    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      if(wire0 < vtx[iv].Wire) continue;
      if(wire0 > vtx[iv].Wire + 10) continue;
      dt = clpar[0] + (vtx[iv].Wire - wire0) * clpar[1] - vtx[iv].Time;
      if(fabs(dt) > 10) continue;
      // cluster points to an US vertex. See if the angle is similar to
      // cluster associated with this vertex
      for(icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].EndVtx != iv) continue;
        if(fabs(tcl[icl].EndAng - thclus) < fKinkAngCut[pass]) return true;
      }
    }
    
    return false;
    
  } // CrawlVtxChk2

//////////////////////////////////////////
  bool ClusterCrawlerAlg::CrawlVtxChk(
    std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<VtxStore>& vtx, unsigned short kwire)
  {
    
    // returns true if the cluster is near a vertex on wire kwire
    if(vtx.size() == 0) return false;
    unsigned short iht = fcl2hits.size() - 1;
    float wire0 = (0.5 + allhits[fcl2hits[iht]].WireNum);
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      if((unsigned short)(0.5 + vtx[iv].Wire) != kwire) continue;
      if(fabs(prtime - vtx[iv].Time) < 10) return true;
    }
    return false;
  }
//////////////////////////////////////////
    void ClusterCrawlerAlg::VtxConstraint(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<VtxStore>& vtx,
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
        clpar[0] = allhits[ihit].Time;
        clpar[1] = (vtx[iv].Time - allhits[ihit].Time) / (vtx[iv].Wire - iwire);
        float prtime = clpar[0] + clpar[1] * (jwire - iwire);
        for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
          if(allhits[jhit].InClus != 0) continue;
          float tdiff = fabs(prtime - allhits[jhit].Time) / allhits[jhit].RMS;
          if(tdiff < 2.5) {
            useHit = jhit;
            doConstrain = true;
            return;
          }
        } // jhit
      } // iv
    } // VtxConstraint


/////////////////////////////////////////
    void ClusterCrawlerAlg::VtxClusterSplit(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
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
            if(fabs(dth) < 0.1) continue;
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
            if(allhits[iht].WireNum >= vtx[ivx].Wire) {
              nSplit = ii;
  if(vtxprt) mf::LogVerbatim("ClusterCrawlerAlg")
    <<"Split cluster "<<tcl[icl].ID<<" at wire "<<allhits[iht].WireNum
    <<" nSplit "<<nSplit;
              ihvx = iht;
              break;
            }
          } // ii
          // found the wire. Now make a rough time cut
          if(ihvx < 0) continue;
          short dtime = abs(allhits[ihvx].Time - vtx[ivx].Time);
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
              if(fabs(tcl[jcl].BeginAng - iclAng) > 0.4) {
                // large angle difference. Set the flag
                angOK = true;
                break;
              }
            } // tcl[jcl].BeginVtx == ivx
            if(tcl[jcl].EndVtx == ivx) {
              if(fabs(tcl[jcl].EndAng - iclAng) > 0.4) {
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
              TmpGet(tcl, icl);
              for(unsigned short ii = 0; ii < nLop; ++ii) {
                unsigned short iht = fcl2hits[fcl2hits.size()-1];
                allhits[iht].InClus = 0;
                fcl2hits.pop_back();
              }
              // store it
              clProcCode += 1000;
              TmpStore(allhits, tcl);
              unsigned short newcl = tcl.size() - 1;
              tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
              tcl[newcl].EndVtx = ivx;
              // declare this cluster obsolete
              tcl[icl].ID = -tcl[icl].ID;
            } else {
              // split the cluster into two
              // correct the split position
              ++nSplit;
              SplitCluster(allhits, tcl, icl, nSplit, ivx);
              tcl[icl].ProcCode += 1000;
              tcl[tcl.size()-1].ProcCode += 1000;
            }
            didit = true;
          } // angOK
          if(didit) break;
        } // ivx
      } // icl

    } // VtxClusterSplit

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeHits(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short theHit)
    {
      // Merge unused separate hits in the multiplet of which 
      // theHit is a member into one hit (= theHit). Set the charge of the 
      // merged hits < 0 to indicate they are obsolete. Hits in the multiplet
      // that are associated with an existing cluster are not affected

      if(theHit > allhits.size() - 1) {
        mf::LogError("ClusterCrawler")<<"Bad theHit";
        return;
      }

  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"MergeHits "
    <<allhits[theHit].WireNum<<":"<<(int)allhits[theHit].Time
    <<" numHits "<<allhits[theHit].numHits;

      // ensure that this is a high multiplicity hit
      if(allhits[theHit].numHits == 1) return;

      // calculate the Charge normalization factor using the hit information
      // instead of passing CCHitFinder ChgNorms all the way down here
      float chgNorm = 2.507 * allhits[theHit].Amplitude * allhits[theHit].RMS
        / allhits[theHit].Charge;

      short loTime = 9999;
      short hiTime = 0;
      unsigned short nGaus = 0;
      // number of hits in this hit multiplet
      unsigned short nhitnear = allhits[theHit].numHits;
      // number of hits that are close to theHit
      float hitSep;
      unsigned short nclose = 0;
      // find the time range for the hit multiplet
      for(unsigned short jj = 0; jj < nhitnear; ++jj) {
        unsigned short jht = allhits[theHit].LoHitID + jj;

  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<" P:W:T "<<plane<<":"<<allhits[jht].WireNum<<":"<<(int)allhits[jht].Time
    <<" Amp "<<(int)allhits[jht].Amplitude
    <<" RMS "<<allhits[jht].RMS
    <<" Charge "<<(int)allhits[jht].Charge
    <<" InClus "<<allhits[jht].InClus;
  }
        // error checking
        if(allhits[jht].LoHitID != allhits[theHit].LoHitID) {
          mf::LogError("ClusterCrawler")<<"Hit multiplet ID error "<<jj<<" "<<allhits[theHit].numHits;
          return;
        }
        // hit is not used by another cluster
        if(allhits[jht].InClus != 0) continue;
        short arg = (short)(allhits[jht].Time - 3 * allhits[jht].RMS);
        if(arg < loTime) loTime = arg;
        arg = (short)(allhits[jht].Time + 3 * allhits[jht].RMS);
        if(arg > hiTime) hiTime = arg;
        ++nGaus;
        hitSep = fabs(allhits[jht].Time - allhits[theHit].Time) / allhits[jht].RMS;
        if(jht != theHit && hitSep < 3) ++nclose;
      } // jj
      // all hits in the multiplet used?
      if(nGaus == 0) return;
      // no reasonably close hits?
//      if(nhitnear < 4 && nclose == 0) return;
      if(loTime < 0) loTime = 0;
      ++hiTime;
      // define a signal shape
      std::vector<double> signal;
      // fill the signal with zeros
      for(unsigned short time = loTime; time < hiTime; ++time) signal.push_back(0.);
      // now add the Gaussians for each hit
      double chgsum = 0.;
      for(unsigned short jj = 0; jj < nhitnear; ++jj) {
        unsigned short jht = allhits[theHit].LoHitID + jj;
        if(jht != theHit) {
          // hit used in another cluster
          if(allhits[jht].InClus != 0) continue;
          // declare this hit obsolete
          allhits[jht].InClus = -1;
        } // jht != theHit
        // add up the charge
        chgsum += allhits[jht].Charge;
        for(unsigned short time = loTime; time < hiTime; ++time) {
          unsigned short indx = time - loTime;
          double arg = (allhits[jht].Time - (double)time) / allhits[jht].RMS;
          signal[indx] += allhits[jht].Amplitude * exp(-0.5 * arg * arg);
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
        mf::LogError("ClusterCrawler")<<"MergeHits: bad sum";
        return;
      }
      double aveTime = sigsumt / sigsum;
      // find the RMS
      sigsumt = 0.;
      for(unsigned short time = loTime; time < hiTime; ++time) {
        double dtime = time - aveTime;
        sigsumt += signal[time - loTime] * dtime * dtime;
      }
      allhits[theHit].RMS = sqrt(sigsumt / sigsum);
      allhits[theHit].Time = aveTime;
      allhits[theHit].Charge = chgsum;
      // find the amplitude from the integrated charge and the RMS
      allhits[theHit].Amplitude = chgsum * chgNorm/ (2.507 * allhits[theHit].RMS);
      allhits[theHit].numHits = 1;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<" theHit "<<allhits[theHit].WireNum<<":"<<(int)aveTime
    <<" RMS "<<std::setprecision(1)<<allhits[theHit].RMS
    <<" chg "<<(int)chgsum<<" Amp "<<(int)allhits[theHit].Amplitude;
  }
  
    } // MergeHits

/////////////////////////////////////////
    void ClusterCrawlerAlg::FindStarVertices(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {
      // Make 2D vertices with a star topology with short back-to-back
      // clusters. Vertices must reside on the US end of the longest
      // cluster, so vertex finding uses the End information only.
      // This routine is called after all passes are completed
      // in the current CTP
      if(tcl.size() < 2) return;

  vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")<<"FindStarVertices";
    PrintClusters(allhits, tcl, vtx);
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
          if(fabs(dth) < 0.1) continue;
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
            if(fabs(dth) < 0.05) continue;
          }
          es2 = tcl[it2].EndSlp;
          eth2 = tcl[it2].EndAng;
          ew2 = tcl[it2].EndWir;
          et2 = tcl[it2].EndTim;
          bw2 = tcl[it2].BeginWir;
          bt2 = tcl[it2].BeginTim;
          // require angle difference
          dth = fabs(eth1 - eth2);
          if(dth < 0.3) continue;
          dsl = es2 - es1;
          fvw = (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
          // intersection within the wire boundaries
          if(fvw < ew1 || fvw > bw1) continue;
          if(fvw < ew2 || fvw > bw2) continue;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
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
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<" vertex time "<<fvt
    <<" lotime "<<lotime<<" hitime "<<hitime;
          unsigned short vw = (0.5 + fvw);
          // ensure that the vertex is near a hit on both clusters
          unsigned short pos1 = 0;
          for(unsigned short ii = 0; ii < tcl[it1].tclhits.size(); ++ii) {
            if(allhits[tcl[it1].tclhits[ii]].WireNum <= vw) {
              if(abs(allhits[tcl[it1].tclhits[ii]].Time - fvt) < 10) pos1 = ii;
              break;
            }
          } // ii
          // vertex is not near a hit on cluster 1
          if(pos1 == 0) continue;
          unsigned short pos2 = 0;
          for(unsigned short ii = 0; ii < tcl[it2].tclhits.size(); ++ii) {
            if(allhits[tcl[it2].tclhits[ii]].WireNum <= vw) {
              if(abs(allhits[tcl[it2].tclhits[ii]].Time - fvt) < 10) pos2 = ii;
              break;
            }
          } // ii
          // vertex is not near a hit on cluster 2
          if(pos2 == 0) continue;
          // ensure we aren't near an existing vertex
          for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
            if(abs(fvw - vtx[iv].Wire) < 2 &&
               abs(fvt - vtx[iv].Time) < 10) continue;
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
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<" new vertex "<<ivx
    <<" cluster "<<tcl[it1].ID<<" split pos "<<pos1;
          SplitCluster(allhits, tcl, it1, pos1, ivx);
          tcl[it1].ProcCode += 1000;
          tcl[tcl.size()-1].ProcCode += 1000;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<" new vertex "<<ivx
    <<" cluster "<<tcl[it2].ID<<" split pos "<<pos2;
          SplitCluster(allhits, tcl, it2, pos2, ivx);
          tcl[it2].ProcCode += 1000;
          tcl[tcl.size()-1].ProcCode += 1000;
          break;
        } // it2
      } // it1
            
      if(vtx.size() > vtxSizeIn) {
        // try to split other clusters
        VtxClusterSplit(allhits, tcl, vtx);
        // try to attach other clusters to it
        VertexCluster(tcl, vtx, vtx.size() - 1);
        VtxWghtAndFit(tcl, vtx, clCTP);
      } // new vertex

  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")<<"Vertices "<<vtx.size();
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      mf::LogVerbatim("ClusterCrawler")
        <<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time
        <<" wght "<<(int)vtx[iv].Wght<<" topo "<<vtx[iv].Topo;
    }
    PrintClusters(allhits, tcl, vtx);
  }
      
    } // FindStarVertices

//////////////////////////////////////////
    void ClusterCrawlerAlg::VertexCluster(std::vector<ClusterStore>& tcl,
        std::vector<VtxStore>& vtx, unsigned short iv)
    {
      // try to attach clusters to the specified vertex
      if(vtx[iv].Wght < 0) return;
      short dwib, dwie;
      float dw, dt;
        
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != vtx[iv].CTP) continue;
        dwib = abs(vtx[iv].Wire - tcl[it].BeginWir);
        dwie = abs(vtx[iv].Wire - tcl[it].EndWir);
        if(dwib < dwie && dwib < 10) {
          // Cluster Begin is closer
          if(vtx[iv].Wire > tcl[it].BeginWir + 1) continue;
          dw = vtx[iv].Wire - tcl[it].BeginWir;
          dt = tcl[it].BeginTim + tcl[it].BeginSlp * dw - vtx[iv].Time;
          if(fabs(dt) > 10) continue;
          tcl[it].BeginVtx = iv;
        } else if(dwie < 10){
          // Cluster End is closer
          if(vtx[iv].Wire < tcl[it].EndWir - 1) continue;
          dw = vtx[iv].Wire - tcl[it].EndWir;
          dt = tcl[it].EndTim + tcl[it].EndSlp * dw - vtx[iv].Time;
          if(fabs(dt) > 10) continue;
          tcl[it].EndVtx = iv;
        }
      } // icl
    } // VertexCluster

//////////////////////////////////////////
    void ClusterCrawlerAlg::VtxWghtAndFit(std::vector<ClusterStore>& tcl,
        std::vector<VtxStore>& vtx, CTP_t inCTP)
      {

      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != inCTP) continue;
        if(vtx[iv].Wght < 0) continue;
        // fit the vertices
        float ChiDOF = 0.;
        FitVtx(tcl, vtx, iv, ChiDOF);
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
    void ClusterCrawlerAlg::KillVertices(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
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
        if(fabs(tcl[it].BeginAng - tcl[it].EndAng) > 0.1) continue;
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
            if(allhits[iht].WireNum <= vWire) {
              if(fabs(allhits[iht].Time - vtx[ivx].Time) < 60) {
                killit = true;
                break;
              }
              break;
            } // allhits[iht].WireNum <= vWire
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
    void ClusterCrawlerAlg::FindVertices(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
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
    mf::LogVerbatim("ClusterCrawler")<<"FindVertices plane "<<plane<<" pass "<<pass;
    PrintClusters(allhits, tcl, vtx);
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
          ClusterVertex(allhits, tcl, vtx, it2);
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
          dth = fabs(eth1 - eth2);
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
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
    <<"Chk clusters topo1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" fvw "<<fvw;
              if(vw >= fFirstWire - 1 && 
                 fvw <= ew1 + 1    && fvw <= ew2 + 1 &&
                 fvw > (ew1 - 10)  && fvw > (ew2 - 10) ) {
                fvt = et1 + (fvw - ew1) * es1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<" Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<"  vtx wire "<<vw<<" time "<<(int)fvt;
  }
                if(fvt > 0 && fvt < maxtime) {
                  // Vertex wire US of cluster ends and time in the detector.
                  // Check for signal at the vertex position and adjust the vertex by 1 wire
                  // if necessary
                  ChkSignal(allhits, vw, fvt, vw, fvt, SigOK);
                  if(!SigOK) {
                    fvw += 1.;
                    vw = (0.5 + fvw);
                    ChkSignal(allhits, vw, fvt, vw, fvt, SigOK);
                  }
                  // Check this against existing vertices and update
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<" SigOK "<<SigOK;
                  if(SigOK) ChkVertex(allhits, tcl, vtx, fvw, fvt, it1, it2, 1);
                } // fvt in detector
              } // vw topo 1 check
            } // fvw in detector
          } // topo 1
  // topo 2: check for vtx US of it1 and DS of it2
          dth = fabs(eth1 - bth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0) {
            // large angle difference
            if(dth > angcut) {
              dsl = bs2 - es1;
              if(fabs(ew1 - bw2) < 3 && fabs(et1 - bt2) < 20) {
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
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters topo2 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo2 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                  if(fvt > 0. && fvt < maxtime) {
                    ChkVertex(allhits, tcl, vtx, fvw, fvt, it1, it2, 2);
                  } // fvt in detector
                } // vw topo 2 check
              } // fvw in detector
            } // dth > 0.3
/*
            else if(bothLong && dth > 0.1) {
              // small angle difference. Require end-to-end clusters with
              // significant charge difference
              dchg = fabs((ec1 - bc2) / (ec1 + bc2));
              if(ew1 > bw2 && ew1 < bw2 + 3 && abs(et1 - bt2) < 10 && dchg > 0.5) {
                fvw = 0.5 *(ew1 + bw2);
                fvt = 0.5 *(et1 + bt2);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters topo2 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo2 vtx wire "<<vw<<" time "<<(int)fvt<<" small angle";
  }
                ChkVertex(allhits, tcl, vtx, fvw, fvt, it1, it2, 2);
              }
            } // dth < 0.3
*/
          } // topo 2
  // topo 3: check for vtx DS of it1 and US of it2
          dth = fabs(bth1 - eth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0) {
            // large angle difference
            if(dth > angcut) {
              dsl = bs1 - es2;
              if(fabs(bw1 - ew2) < 3 && fabs(bt1 - et2) < 20) {
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
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters topo3 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo3 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                  if(fvt > 0. && fvt < maxtime) {
                    ChkVertex(allhits, tcl, vtx, fvw, fvt, it1, it2, 3);
                  } // fvt in detector
                } // vw topo 3 check
              } // fvw in detector
            } 
/*
            else if(bothLong && dth > 0.1) {
              // small angle difference
              dchg = fabs((bc1 - ec2) / (bc1 + ec2));
              if(ew2 > bw1 && ew2 < bw1 + 3 && abs(et1 - bt2) < 10 && dchg > 0.6) {
                fvw = 0.5 *(ew2 + bw1);
                fvt = 0.5 *(et2 + bt1);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters topo3 "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo3 vtx wire "<<vw<<" time "<<(int)fvt<<" small angle";
  }
                ChkVertex(allhits, tcl, vtx, fvw, fvt, it1, it2, 3);
              }
            } // dth < 0.3
*/
          } // topo 3
  // topo 4: check for vtx DS of it1 and DS of it2
          dth = fabs(bth1 - bth2);
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
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo4 vtx wire "<<vw<<" time "<<fvt;
  }
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check for signal at the vertex position and adjust the vertex by 1 wire
                  // if necessary
                  ChkSignal(allhits, vw, fvt, vw, fvt, SigOK);
                  if(!SigOK) {
                    fvw -= 1.;
                    vw = (0.5 + fvw);
                    ChkSignal(allhits, vw, fvt, vw, fvt, SigOK);
                  }
                  // Check this against existing vertices and update
                  if(SigOK) ChkVertex(allhits, tcl, vtx, fvw, fvt, it1, it2, 4);
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

      if(vtx.size() > vtxSizeIn) VtxWghtAndFit(tcl, vtx, clCTP);
      
    } // FindVertices

/////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterVertex(std::vector<CCHitFinderAlg::CCHit>& allhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned short it)
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
          dwib = abs(vtx[iv].Wire - tcl[it].BeginWir);
          if(dwib > 2) dwib = 2;
          dwie = abs(vtx[iv].Wire - tcl[it].EndWir);
          if(dwie > 2) dwie = 2;
          dwjb = 999; dwje = 999;
          for(jv = 0; jv < vtx.size(); ++jv) {
            if(iv == jv) continue;
            if(fabs(vtx[jv].Time - tcl[it].BeginTim) < 20) {
              if(abs(vtx[jv].Wire - tcl[it].BeginWir) < dwjb) 
                dwjb = abs(vtx[jv].Wire - tcl[it].BeginWir);
            } // fabs(vtx[jv].Time - tcl[it].BeginTim) < 20
            if(fabs(vtx[jv].Time - tcl[it].EndTim) < 20) {
              if(abs(vtx[jv].Wire - tcl[it].EndWir) < dwje) 
                dwje = abs(vtx[jv].Wire - tcl[it].EndWir);
            } // fabs(vtx[jv].Time - tcl[it].EndTim) < 20
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
          if(fabs(tdiff) < 10) {
            SigOK = false;
            ChkSignal(allhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].EndVtx = iv;
              // re-fit it
              ChiDOF = 99.;
              FitVtx(tcl, vtx, iv, ChiDOF);
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Add End "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<ChiDOF;
              return;
            } // SigOK
          } // tdiff
        } // matchUS
        
        if(matchDS) {
          // project cluster to DS vertex
          tdiff = tcl[it].BeginTim + (vtx[iv].Wire - tcl[it].BeginWir) * tcl[it].BeginSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            SigOK = false;
            ChkSignal(allhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].BeginVtx = iv;
              // re-fit it
              ChiDOF = 99.;
              FitVtx(tcl, vtx, iv, ChiDOF);
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Add Begin "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<ChiDOF;
              return;
            } // SigOK
          } // tdiff
        } // matchDS
        
      } // iv
    } // ClusterVertex



/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkVertex(std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
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
    mf::LogVerbatim("ClusterCrawler")
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
            if(fabs(tcl[it].BeginAng - tcl[it].EndAng) > 0.1) continue;
            // don't reject because it is near the end of a long cluster
            if(vw <  tcl[it].EndWir + 5) continue;
            if(vw >  tcl[it].BeginWir - 5) continue;
            for(unsigned short ii = 0; ii < tcl[it].tclhits.size(); ++ii) {
              iht = tcl[it].tclhits[ii];
              if(allhits[iht].WireNum <= vw) {
                if(fabs(allhits[iht].Time - fvt) < 60) return;
              } // allhits[iht].WireNum <= vWire
            } // ii
          } // it
        }

        bool SigOK = false;
        // check vertex and clusters for proximity to existing vertices
        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if(vtx[iv].CTP != clCTP) continue;
          if( abs(fvw - vtx[iv].Wire) < 4 &&
             fabs(fvt - vtx[iv].Time) < 25) {
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Match to vtx "<<iv;
            // got a match. Check the appropriate cluster end and attach
            if( (topo == 1 || topo == 2) && tcl[it1].EndVtx < 0) {
              ChkSignal(allhits, vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, SigOK);
              if(SigOK) tcl[it1].EndVtx = iv;
  if(vtxprt)  mf::LogVerbatim("ClusterCrawler")<<"12 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv;
            } else if( (topo == 3 || topo == 4) && tcl[it1].BeginVtx < 0) {
              ChkSignal(allhits, vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, SigOK);
              if(SigOK) tcl[it1].BeginVtx = iv;
  if(vtxprt)  mf::LogVerbatim("ClusterCrawler")<<"34 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv;
            } // cluster it2
            if( (topo == 1 || topo == 3) && tcl[it2].EndVtx < 0) {
              ChkSignal(allhits, vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, SigOK);
              if(SigOK) tcl[it2].EndVtx = iv;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"13 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv;
            } else if( (topo == 2 || topo == 4) && tcl[it2].BeginVtx < 0) {
              ChkSignal(allhits, vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, SigOK);
              if(SigOK) tcl[it2].BeginVtx = iv;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"24 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv;
            } // cluster it2
            return;
          } // matched vertex
        } // iv

        // no match to existing vertices. Ensure that there is a wire signal between
        // the vertex and the appropriate ends of the clusters
        bool Sig1OK = false;
        bool Sig2OK = false;
        if(topo == 1 || topo == 2) {
          ChkSignal(allhits, vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, Sig1OK);
        } else {
          ChkSignal(allhits, vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, Sig1OK);
        }
        if(topo == 1 || topo == 3) {
          ChkSignal(allhits, vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, Sig2OK);
        } else {
          ChkSignal(allhits, vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, Sig2OK);
        }
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
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
    mf::LogVerbatim("ClusterCrawler")
      <<"New vtx "<<iv<<" in plane "<<plane<<" topo "<<topo<<" cls "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" time "<<(int)fvt<<" wire "<<vw;
  }
        }

      }

/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkSignal(std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short wire1, float time1, unsigned short wire2, float time2, bool& SigOK)
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
//  mf::LogVerbatim("ClusterCrawler")<<"ChkSignal "<<wiree<<" "<<wire<<" "<<wireb<<" "<<(int)prtime;
//    <<" first "<<firsthit<<" last "<<lasthit;
        for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
          // outside the hit RAT?
          if(prTimeHi > allhits[khit].HiTime) continue;
          if(prTimeLo < allhits[khit].LoTime) continue;
          // found a signal. Skip checking on this wire
          WireSigOK = true;
          break;
        } // khit
//        mf::LogVerbatim("ClusterCrawler")<<" SigOK "<<WireSigOK;
        if(!WireSigOK) return;
      } // wire
      SigOK = true;
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::SplitCluster(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, 
        unsigned short icl, unsigned short pos, unsigned short ivx)
    {
      // split cluster icl into two clusters

      if(tcl[icl].ID < 0) {
        mf::LogError("ClusterCrawler")<<"Trying to split obsolete cluster "
          <<icl;
        return;
      }

      if(pos > tcl[icl].tclhits.size()) {
        mf::LogError("ClusterCrawler")<<"SplitCluster bad split position "
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
      FitCluster(allhits);
      clEndSlp = clpar[1];
      clEndSlpErr = clparerr[1];
      clEndAng = atan(fScaleF * clEndSlp);
      // find the charge at the end
      FitClusterChg(allhits);
      clEndChg = fAveChg;
      TmpStore(allhits, tcl);
      // associate the End with the supplied vertex
      unsigned short iclnew = tcl.size() - 1;
      tcl[iclnew].EndVtx = ivx;
      tcl[iclnew].BeginVtx = tcl[icl].BeginVtx;

      // now create the second cluster (US)
      clEndSlp = tcl[icl].EndSlp;
      clEndSlpErr = tcl[icl].EndSlpErr;
      clEndAng = atan(fScaleF * clEndSlp);
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
  if(allhits[iht].InClus != tcl[icl].ID) {
    mf::LogError("ClusterCrawler")
      <<"SplitCluster bad hit "<<iht<<" "<<allhits[iht].InClus
      <<" "<<tcl[icl].ID<<" ProcCode "<<clProcCode<<std::endl;
  }
        fcl2hits.push_back(iht);
        // define the Begin parameters
        if(fcl2hits.size() == fMaxHitsFit[pass] ||
           fcl2hits.size() == fMinHits[pass]) {
          FitCluster(allhits);
          clBeginSlp = clpar[1];
          clBeginAng = atan(fScaleF * clBeginSlp);
          clBeginSlpErr = clparerr[1];
          didFit = true;
        }
        if((unsigned short)fcl2hits.size() == fNHitsAve[pass] + 1) {
          FitClusterChg(allhits);
          clBeginChg = fAveChg;
          didFit = true;
        }
      } // ii
      // do a fit using all hits if one wasn't done
      if(!didFit) {
        FitCluster(allhits);
        FitClusterChg(allhits);
        clBeginChg = fAveChg;
      }
      TmpStore(allhits, tcl);
      // associate the End with the supplied vertex
      iclnew = tcl.size() - 1;
      tcl[iclnew].BeginVtx = ivx;
      tcl[iclnew].EndVtx = tcl[icl].EndVtx;
//  mf::LogVerbatim("ClusterCrawler")<<"ClusterSplit split "<<tcl[icl].ID;
      // declare icl obsolete
      tcl[icl].ID = -tcl[icl].ID;

    } // SplitCluster

/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkMerge(std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl)
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
        bth1 = atan(fScaleF * bs1);
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
          bth2 = atan(fScaleF * bs2);
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
          timecut *= (2 - 1/(1 + fabs(clpar[1])));

          // tweak the cuts for long straight tracks
          bothLong = (nh1 > 100 && tcl[it2].tclhits.size() > 100);
          
          // look for US and DS broken clusters w similar angle.
          // US cluster 2 merge with DS cluster 1?
          // This is the most likely occurrence given the order in which
          // clusters are created so put it first.
          dth = fabs(bth2 - eth1);
          // require no vertex between
          NoVtx = (tcl[it1].EndVtx < 0) && (tcl[it2].BeginVtx < 0);
  if(prt && bw2 < ew1 ) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(ew1 - bw2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          if(NoVtx && bw2 < ew1 && (ew1 - bw2)  < skipcut && dth < angcut) {
            chgrat = 2 * fabs(bc2 - ec1) / (bc2 + ec1);
            // ignore the charge cut for large angle clusters
            if(fabs(es1) > fLAClusSlopeCut) chgrat = 0.;
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.05) chgrat = 0.;
            // project bw2,bt2 to ew1
            dtim = fabs(bt2 + (ew1-bw2)*bs2 - et1);
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<" dtim "<<dtim<<" timecut "<<(int)timecut
    <<" ec1 "<<ec1<<" bc2 "<<bc2
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              // ensure there is a signal between cluster ends
              SigOK = true;
              ChkSignal(allhits,ew1,et1,bw2,bt2,SigOK);
              if(SigOK) {
                DoMerge(allhits, tcl, vtx, it2, it1, 10);
                tclsize = tcl.size();
                break;
              }
            }
          } // US cluster 2 merge with DS cluster 1?
          
          // look for US and DS broken clusters w similar angle
          // US cluster 1 merge with DS cluster 2?
          dth = fabs(bth1 - eth2);
  if(prt && bw1 < ew2 && (ew2 - bw1)  < skipcut) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk2 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(bw1 - ew2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          // require no vertex between
          NoVtx = (tcl[it2].EndVtx < 0) && (tcl[it1].BeginVtx < 0);
          if(NoVtx && bw1 < ew2 && (ew2 - bw1)  < skipcut && dth < angcut ) {
            chgrat = 2 * fabs((bc1 - ec2) / (bc1 + ec2));
            // ignore the charge cut for large angle clusters
            if(fabs(es2) > fLAClusSlopeCut) chgrat = 0.;
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.05) chgrat = 0.;
            // project sw1,st1 to ew2
            dtim = fabs(bt1 + (ew2-bw1)*bs1 - et2);
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<" dtim "<<dtim<<" err "<<dtim<<" timecut "<<(int)timecut
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              SigOK = true;
              ChkSignal(allhits,bw1,bt1,ew2,et2,SigOK);
              if(SigOK) {
                DoMerge(allhits, tcl, vtx, it1, it2, 10);
                tclsize = tcl.size();
                break;
              }
            }
          } // US cluster 1 merge with DS cluster 2
          
          if(bw2 < bw1 && ew2 > ew1) {
            // look for small cl2 within the wire boundary of cl1
            // with similar times and slopes for both clusters
            dth = fabs(eth2 - eth1);
            dtim = fabs(et1 +(ew2 - ew1 - 1)*es1 - et2);
            // count the number of wires with no hits on cluster 1
            short nmiss1 = bw1 - ew1 + 1 - tcl[it1].tclhits.size();
            // compare with the number of hits in cluster 2
            short nin2 = tcl[it2].tclhits.size();
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<"cl2: "<<tcl[it2].ID<<" within cl1 "<<tcl[it1].ID
    <<" ? dth "<<dth<<" dtim "<<dtim<<" nmissed "<<nmiss1;
  }
            // make rough cuts before calling ChkMerge12
            // this may not work well for long wandering clusters
            bool didit = false;
            if(dth < 1 && dtim < timecut && nmiss1 >= nin2) 
                ChkMerge12(allhits, tcl, it1, it2, didit);
            if(didit) {
              tclsize = tcl.size();
              break;
            } //didit
          } // small cl2 within the wire boundary of cl1
          
          if(bw1 < bw2 && ew1 > ew2) {
            // look for small cl1 within the wire boundary of cl2
            // with similar times and slopes for both clusters
            dth = fabs(eth2 - eth1);
            dtim = fabs(et2 +(ew1 - ew2 - 1)*es2 - et1);
            // count the number of wires with no hits on cluster 2
            short nmiss2 = bw2 - ew2 + 1 - tcl[it2].tclhits.size();
            // compare with the number of hits in cluster 1
            short nin1 = tcl[it1].tclhits.size();
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<"cl1: "<<tcl[it1].ID<<" within cl2 "<<tcl[it2].ID
    <<" ? dth "<<dth<<" dtim "<<dtim<<" nmissed "<<nmiss2;
  }
            // make rough cuts before calling ChkMerge12
            // this may not work well for long wandering clusters
            bool didit = false;
            if(dth < 1 && dtim < timecut && nmiss2 >= nin1) 
                ChkMerge12(allhits, tcl, it2, it1, didit);
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
  void ClusterCrawlerAlg::ChkMerge12(std::vector<CCHitFinderAlg::CCHit>& allhits, 
     std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short it2, bool& didit)
  {
    // Calling routine has done a rough check that cluster it2 is a candidate
    // for merging with cluster it1. The wire range spanned by it2 lies
    // within the wire range of it1 and the clusters are reasonably close
    // together in time.
    
    // assume that no merging was done
    didit = false;
    
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"ChkMerge12 "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    
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
      unsigned short wire = allhits[hit].WireNum;
      if(wire - ew1 < 0 || wire - ew1 > (short)cl1hits.size()) {
        mf::LogError("ClusterCrawler")<<"ChkMerge12 bad wire "<<(wire-ew1);
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"chk next US wire "<<wiron1<<" missed "<<nmiss;
    if(wiron1 == 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    
    // compare the wires with hits on cluster 2 with the gap in cluster 1
    // the number of hit wires that fit in the gap
    unsigned short nfit = 0;
    for(unsigned short iht = 0; iht < tcl[it2].tclhits.size(); ++iht) {
      unsigned short hiton2 = tcl[it2].tclhits[iht];
      unsigned short wiron2 = allhits[hiton2].WireNum;
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
    if(hiton1 > allhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    // check the time difference
    float timon1 = allhits[hiton1].Time;
    float dtim = fabs(et2 + (wiron1 - ew2) * tcl[it2].EndSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    // check the slope difference. First do a local fit on cluster 1 near
    // the matching point
    FitClusterMid(allhits, tcl, it1, hiton1, 3);
    if(clChisq > 20.) return;
    // fit parameters are now in clpar.
    // check for angle consistency
    float dth = fabs(atan(fScaleF * clpar[1]) - atan(fScaleF * tcl[it2].EndSlp));
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"US dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut. fAveChg was calculated in FitClusterMid
    float chgrat = 2 * fabs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"US chgrat "<<chgrat<<" cut "<<fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    bool SigOK = false;
    ChkSignal(allhits, wiron1, timon1, ew2, et2, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"US SigOK? "<<SigOK;
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
    if(hiton1 > allhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    timon1 = allhits[hiton1].Time;
    dtim = fabs(bt2 - (wiron1 - bw2) *tcl[it2].BeginSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    FitClusterMid(allhits, tcl, it1, hiton1, -3);
    if(clChisq > 20.) return;
    // check for angle consistency
    dth = fabs(atan(fScaleF * clpar[1]) - atan(fScaleF * tcl[it2].BeginSlp));
  if(prt) mf::LogVerbatim("ClusterCrawler")
      <<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * fabs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(prt) mf::LogVerbatim("ClusterCrawler")
      <<"DS chgrat "<<chgrat<<" cut "<<fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    SigOK = false;
    ChkSignal(allhits, wiron1, timon1, bw2, bt2, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"DS SigOK? "<<SigOK;
    if( !SigOK ) return;

  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Merge em";
    // success. Merge them
    DoMerge(allhits, tcl, vtx, it1, it2, 100);
    didit = true;
  } // ChkMerge12

/////////////////////////////////////////
  void ClusterCrawlerAlg::DoMerge(
      std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      unsigned short it1, unsigned short it2, short inProcCode)
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
      allhits[hit].InClus = 0;
      unsigned short wire = allhits[hit].WireNum;
      unsigned short index = wire - lowire;
      wirehit[index] = hit;
/*
  if(myprt) mf::LogVerbatim("ClusterCrawler")
    <<"Cl2 hit "<<wire<<":"<<(int)allhits[hit].Time
    <<" wire index "<<index<<" hit index "<<hit;
*/
     } // iht
    // now cluster 1
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned short hit = cl1.tclhits[iht];
      unsigned short wire = allhits[hit].WireNum;
      unsigned short index = wire - lowire;
      allhits[hit].InClus = 0;
      // TODO: Should merge hits instead of clobbering one of them?
      wirehit[index] = hit;
/*
  if(myprt) mf::LogVerbatim("ClusterCrawler")
    <<"Cl1 hit "<<wire<<":"<<(int)allhits[hit].Time
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
    TmpStore(allhits, tcl);
    unsigned short itnew = tcl.size()-1;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"DoMerge "<<cl1.ID<<" "<<cl2.ID<<" -> "<<tcl[itnew].ID;
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
  void ClusterCrawlerAlg::PrintClusters(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
  {

    float aveRMS = 0;
    // prints clusters to the screen for code development
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"  ID CTP nht Stop  Proc  beg_W:T  bTheta"
      <<" begChg end_W:T  eTheta eChg  bVx  eVx aveRMS\n";
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      if(fDebugPlane >= 0 && fDebugPlane != (int) tcl[ii].CTP) continue;
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
        aveRMS += allhits[iht].RMS;
      }
      aveRMS /= (float)tcl[ii].tclhits.size();
      myprt<<std::right<<std::setw(5)<<std::fixed<<std::setprecision(1)<<aveRMS;
      myprt<<"\n";
    } // ii
    // print out vertices
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      mf::LogVerbatim("ClusterCrawler")
        <<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time
        <<" wght "<<(int)vtx[iv].Wght<<" topo "<<vtx[iv].Topo;
    }    
    // Check for incompatible hit->cluster cluster->hit associations
    for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if(tcl[icl].ID < 0) continue;
      for(unsigned short ii = 0; ii< tcl[icl].tclhits.size(); ++ii) {
        unsigned short iht = tcl[icl].tclhits[ii];
        if(allhits[iht].InClus != tcl[icl].ID) {
          mf::LogVerbatim("PrintClusters")
            <<"Association error: cluster "<<tcl[icl].ID
            <<" Hit "<<allhits[iht].WireNum<<":"<<(int)allhits[iht].Time
            <<" InClus is incorrect "<<allhits[iht].InClus
            <<" Hit index "<<iht;
          return;
        }
      } // ii
    } // icl

  } // PrintClusters

/////////////////////////////////////////
    void ClusterCrawlerAlg::TmpGet(std::vector<ClusterStore>& tcl, 
        unsigned short it1)
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
  void ClusterCrawlerAlg::TmpStore(std::vector<CCHitFinderAlg::CCHit>& allhits, 
     std::vector<ClusterStore>& tcl)
  {

    if(fcl2hits.size() < 3) return;
    
    if(fcl2hits.size() == USHRT_MAX) return;
    
    ++NClusters;

    // flag all the hits as used
    for(unsigned short it = 0; it < fcl2hits.size(); ++it) {
      unsigned short hit = fcl2hits[it];
      if(allhits[hit].InClus < 0) {
        mf::LogError("ClusterCrawler")<<"Trying to use obsolete hit "<<hit
        <<" on wire "<<allhits[hit].WireNum<<" on cluster "<<NClusters
        <<" in plane "<<plane<<" ProcCode "<<clProcCode;
      }
      allhits[hit].InClus = NClusters;
    }

    // ensure that the cluster begin/end info is correct

    // define the begin/end charge if it wasn't done already
    if(clEndChg <= 0) {
      // use the next to the last two hits. The End hit may have low charge
      unsigned int ih0 = fcl2hits.size() - 2;
      unsigned int hit = fcl2hits[ih0];
      clEndChg = allhits[hit].Charge;
      hit = fcl2hits[ih0 - 1];
      clEndChg += allhits[hit].Charge;
      clEndChg = clEndChg / 2.;
    }
    if(clBeginChg <= 0) {
      // use the 2nd and third hit. The Begin hit may have low charge
      unsigned int hit = fcl2hits[1];
      clBeginChg = allhits[hit].Charge;
      hit = fcl2hits[2];
      clBeginChg += allhits[hit].Charge;
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
    clstr.BeginAng    = atan(fScaleF * clBeginSlp);
    clstr.BeginWir    = allhits[hitb].WireNum;
    clstr.BeginTim    = allhits[hitb].Time;
    clstr.BeginChg    = clBeginChg;
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndAng      = atan(fScaleF * clEndSlp);
    clstr.EndWir      = allhits[hite].WireNum;
    clstr.EndTim      = allhits[hite].Time;
    clstr.EndChg      = clEndChg;
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.BeginVtx    = -99;
    clstr.EndVtx      = -99;
    clstr.CTP         = clCTP;
    clstr.tclhits     = fcl2hits;
    tcl.push_back(clstr);
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::LACrawlUS(std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<VtxStore>& vtx)
  {
    // Crawl a large angle cluster upstream. Similar to CrawlUS but require
    // that a hit be added on each wire


    unsigned short dhit = fcl2hits[0];
    short dwir = allhits[dhit].WireNum;
    prt = false;
  if(fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
    prt = abs(allhits[dhit].Time - fDebugHit) < 20;

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"******************* LACrawlUS PASS "<<pass<<" Hits ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned short iht = fcl2hits[fcl2hits.size() - 1 - ii];
      myprt<<allhits[iht].WireNum<<":"<<(int)allhits[iht].Time<<" ";
    }
  }

    // Merge all the hits if this is a very large angle cluster
    if(fabs(clpar[1]) > 5) {
      unsigned short ii, nmult = 0;
      for(ii = 0; ii < fcl2hits.size(); ++ii) 
        if(allhits[fcl2hits[ii]].numHits > 1) ++nmult;
      if(nmult == fcl2hits.size()) {
        for(ii = 0; ii < fcl2hits.size(); ++ii) {
          MergeHits(allhits, fcl2hits[ii]);
        } // ii
        FitCluster(allhits);
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
    unsigned short lastwire = allhits[lasthit].WireNum;
    bool ChkCharge = false;
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"LACrawlUS: next wire "<<nextwire;
      // stop crawling if there is a nearby vertex
      if(CrawlVtxChk(allhits, vtx, nextwire)) {
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"LACrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // AddLAHit will merge the hit on nextwire if necessary
      AddLAHit(allhits, vtx, nextwire, ChkCharge, HitOK, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"LACrawlUS: HitOK "<<HitOK<<" SigOK "<<SigOK;
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
          MergeHits(allhits, hit);
        }
        // update the fit
        FitCluster(allhits);
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
        mf::LogError("ClusterCrawler")
          <<"LACrawlUS: chifits size error "<<chifits.size()<<" "<<fcl2hits.size();
        return;
      }
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
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
        float dt02 = allhits[hit2].Time - allhits[hit0].Time;
        float dw02 = allhits[hit2].WireNum - allhits[hit0].WireNum;
        float th02 = atan( fScaleF * dt02 / dw02);
        // and the 3rd and 5th hit
        unsigned short ih3 = ih2 - 1;
        unsigned short hit3 = fcl2hits[ih3];
        unsigned short ih5 = ih3 - 2;
        unsigned short hit5 = fcl2hits[ih5];
        float dt35 = allhits[hit5].Time - allhits[hit3].Time;
        float dw35 = allhits[hit5].WireNum - allhits[hit3].WireNum;
        float th35 = atan(fScaleF * dt35 / dw35);
        float dth = fabs(th02 - th35);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" Kink angle "<<std::setprecision(3)<<dth<<" cut "<<fKinkAngCut[pass];
        if(dth > fKinkAngCut[pass]) {
          // hit a kink. Lop of the first 3 hits, refit and keep crawling?
          for(short jj = 0; jj < 3; ++jj) {
            fcl2hits.pop_back();
            chifits.pop_back();
            hitnear.pop_back();
          }
          FitCluster(allhits);
          // See if this is a second kink and it is close to the first
          // kink (which had hits removed).
          if(kinkOnWire > 0) {
            if(kinkOnWire - nextwire < 4) {
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Hit a second kink. kinkOnWire = "<<kinkOnWire<<" Stopping";
              // set the kink stop code
              clStopCode = 3;
              break;
            }
          }
          kinkOnWire = nextwire;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Removed kink hits";
        } // kinkang check
      } // chifits test
    } // nextwire

    CheckClusterHitFrac(allhits, prt);

    clProcCode += 300;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"LACrawlUS done. Nhits = "<<fcl2hits.size();
    prt = false;
  } // LACrawlUS

/////////////////////////////////////////
  void ClusterCrawlerAlg::CrawlUS(std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<VtxStore>& vtx)
  {
    // Crawl along a trail of hits moving upstream

    if(fcl2hits.size() < 2) return;

    unsigned short dhit = fcl2hits[0];
    short dwir = allhits[dhit].WireNum;
    prt = false;
  if(fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
    prt = abs(allhits[dhit].Time - fDebugHit) < 20;

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"******************* CrawlUS PASS "<<pass<<" Hits: ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned short iht = fcl2hits[fcl2hits.size() - 1 - ii];
      myprt<<allhits[iht].WireNum<<":"<<(int)allhits[iht].Time<<" ";
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
    if(lasthit > allhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"CrawlUS bad lasthit "<<lasthit;
    }
    unsigned short lastwire = allhits[lasthit].WireNum;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"CrawlUS: last wire "<<lastwire<<" hit "<<lasthit;
    
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"CrawlUS: next wire "<<nextwire;
      // stop crawling if there is a nearby vertex
      if(CrawlVtxChk(allhits, vtx, nextwire)) {
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"CrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // add hits and check for PH and width consistency
      AddHit(allhits, nextwire, HitOK, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"CrawlUS: HitOK "<<HitOK<<" SigOK "<<SigOK;
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" PostSkip && nmissed = "<<nmissed;
            clStopCode = 2;
            for(short jj = 0; jj < nHitAfterSkip; ++jj) {
              fcl2hits.pop_back();
              chifits.pop_back();
              hitnear.pop_back();
            } // pop_back
            FitCluster(allhits);
            if(clChisq > 90.) {
              fcl2hits.clear();
              return;
            }
            FitCluster(allhits);
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
          if(prt) mf::LogVerbatim("ClusterCrawler")<<"No hit or signal on wire "<<nextwire;
          break;
        } // else SigOK false
      } // !HitOK
      else {
        if(clChisq > 99.) {
          if(fcl2hits.size() < 3) return;
          // a fit error occurred. Lop off the leading hit and refit
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Fit failed ";
          fcl2hits.pop_back();
          chifits.pop_back();
          hitnear.pop_back();
          FitCluster(allhits);
          if(clChisq > 99.) {
            // something really bad happened. Bail out
            fcl2hits.clear();
            return;
          }
          FitClusterChg(allhits);
          continue;
        } // clChisq > 99
        // monitor the onset of a kink. Look for a progressive increase
        // in chisq for the previous 0 - 2 hits.
        if(chifits.size() > 5 && fKinkChiRat[pass] > 0) {
          if(chifits.size() != fcl2hits.size()) {
            mf::LogError("ClusterCrawler")
              <<"CrawlUS: chifits size error "<<chifits.size()<<" "<<fcl2hits.size();
            return;
          }
          unsigned short chsiz = chifits.size()-1;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
          if( chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
              chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
            if(fcl2hits.size() != chifits.size()) {
              mf::LogError("ClusterCrawler")
              <<"bad kink check size "<<chifits.size()<<" "<<fcl2hits.size()
              <<" plane "<<plane<<" cluster "<<dwir<<":"<<dhit;
              continue;
            }
            // find the kink angle (crudely) from the 0th and 2nd hit
            unsigned short ih0 = fcl2hits.size() - 1;
            unsigned short hit0 = fcl2hits[ih0];
            unsigned short ih2 = ih0 - 2;
            unsigned short hit2 = fcl2hits[ih2];
            float dt02 = allhits[hit2].Time - allhits[hit0].Time;
            float dw02 = allhits[hit2].WireNum - allhits[hit0].WireNum;
            float th02 = atan( fScaleF * dt02 / dw02);
            // and the 3rd and 5th hit
            unsigned short ih3 = ih2 - 1;
            unsigned short hit3 = fcl2hits[ih3];
            unsigned short ih5 = ih3 - 2;
            unsigned short hit5 = fcl2hits[ih5];
            float dt35 = allhits[hit5].Time - allhits[hit3].Time;
            float dw35 = allhits[hit5].WireNum - allhits[hit3].WireNum;
            float th35 = atan(fScaleF * dt35 / dw35);
            float dth = fabs(th02 - th35);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" Kink angle "<<std::setprecision(3)<<dth<<" cut "<<fKinkAngCut[pass];
            // cut on the allowed kink angle
            if(dth > fKinkAngCut[pass]) {
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" Stopped tracking ";
              // kill the last 3 hits, refit and return
              for(short jj = 0; jj < 3; ++jj) {
                fcl2hits.pop_back();
                chifits.pop_back();
                hitnear.pop_back();
              }
              FitCluster(allhits);
              FitClusterChg(allhits);
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" Set clBeginChg "<<clBeginChg;
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"chirat "<<chirat
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Bad fit chisq - short cluster. Break";
            break;
          }
          FitCluster(allhits);
          FitClusterChg(allhits);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Bad fit chisq - removed hits. Continue...";
        } // clChisq > fChiCut[pass]
      } // !HitOK check
    } // nextwire

    // count the number of hits on adjacent wires at the leading edge and
    // ensure that the count is consistent with fMinWirAfterSkip
    bool reFit = false;
    if((unsigned short)fcl2hits.size() > fMinWirAfterSkip[pass] + 3) {
      unsigned short ih0 = fcl2hits.size() - 1;
      unsigned short hit0 = fcl2hits[ih0];
      unsigned short uswir = allhits[hit0].WireNum;
      unsigned short nAdjHit = 0;
      for(unsigned short ii = ih0 - 1; ii > 0; --ii) {
        unsigned short nxtwir = allhits[ fcl2hits[ii] ].WireNum;
        if(nxtwir != uswir + 1) break;
        ++nAdjHit;
        // break if there are enough hits
        if( nAdjHit == fMinWirAfterSkip[pass] ) break;
        uswir = nxtwir;
      } // ii
      // lop off hits?
      if(nAdjHit < fMinWirAfterSkip[pass]) {
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Lop off "<<nAdjHit<<" hits ";
        for(unsigned short ii = 0; ii < nAdjHit + 1; ++ii) {
          fcl2hits.pop_back();
          chifits.pop_back();
          hitnear.pop_back();
          reFit = true;
        }
      }
    } // fcl2hits.size() > fMinWirAfterSkip[pass] + 3
    if(reFit) {
      FitCluster(allhits);
      FitClusterChg(allhits);
    }
    
    CheckClusterHitFrac(allhits, prt);

    prt = false;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::CheckClusterHitFrac(
      std::vector<CCHitFinderAlg::CCHit>& allhits, bool prt)
  {


    // Find the fraction of the wires on the cluster that have
    // hits.
    unsigned short iht = fcl2hits[fcl2hits.size() - 1];
    clEndWir = allhits[iht].WireNum;
    float hitFrac = (float)fcl2hits.size() / (float)(clBeginWir - clEndWir + 1);

    if(hitFrac < 0.7) {
      fcl2hits.clear();
      if(prt) mf::LogVerbatim("ClusterCrawler")
        <<"CheckClusterHitFrac: Poor hit fraction "<<hitFrac;
      return;
    } // hitFrac
    
    // lop off the last hit if it is part of a hit multiplet
    if(allhits[iht].numHits > 1) {
      fcl2hits.resize(fcl2hits.size() - 1);
    }

    // check for short track ghosts
    if(fcl2hits.size() > 4) return;
    unsigned short nsing = 0;
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      iht = fcl2hits[ii];
      if(allhits[iht].numHits == 1) ++nsing;
    } // ii
    hitFrac = (float)nsing / (float)fcl2hits.size();
    
    if(hitFrac < 0.7) {
      fcl2hits.clear();
      if(prt) mf::LogVerbatim("ClusterCrawler")
        <<"CheckClusterHitFrac: Poor short track hit fraction "<<hitFrac;
      return;
    } // hitFrac
        
  } // CheckClusterHitFrac

/////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterMid(
    std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short ihtin,
    short nhit)
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
        if(ihit > allhits.size()-1) {
          mf::LogError("ClusterCrawler")<<"FitClusterMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = allhits[ihit].WireNum;
        }
        // use hits after finding the first desired hit
        if(UseEm) {
          unsigned short wire = allhits[ihit].WireNum;
          xwir.push_back(wire - wire0);
          ytim.push_back(allhits[ihit].Time);
          // pass the error^2 to the fitter
          float terr = fHitErrFac * allhits[ihit].RMS;
	  ytimerr2.push_back(terr * terr);
          fAveChg += allhits[ihit].Charge;
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
        if(ihit > allhits.size()-1) {
          mf::LogVerbatim("ClusterCrawler")<<"FitClusterMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = allhits[ihit].WireNum;
        }
        // use hits after finding the first desired hit

        if(UseEm) {
          unsigned short wire = allhits[ihit].WireNum;
          xwir.push_back(wire - wire0);
          ytim.push_back(allhits[ihit].Time);
          float terr = fHitErrFac * allhits[ihit].RMS;
	  ytimerr2.push_back(terr * terr);
          fAveChg += allhits[ihit].Charge;
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
  void ClusterCrawlerAlg::FitCluster(std::vector<CCHitFinderAlg::CCHit>& allhits)
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
      unsigned short wire = allhits[ihit].WireNum;
      if(first) {
        wire0 = wire;
        first = false;
      }
      xwir.push_back(wire - wire0);
      ytim.push_back(allhits[ihit].Time);
      // Scale error by hit multiplicity to account for bias in hit
      // multiplet fitting
      float terr = fHitErrFac * allhits[ihit].RMS * allhits[ihit].numHits;
      ytimerr2.push_back(angfactor * terr * terr);
      if(iht == nht) break;
      ++iht;
    }

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"FitCluster W:T ";
    unsigned short cnt = 0;
    for(std::vector<unsigned short>::reverse_iterator it = fcl2hits.rbegin();
       it != fcl2hits.rend(); ++it) {
      unsigned short ihit = *it;
      unsigned short wire = allhits[ihit].WireNum;
      myprt<<wire<<":"<<(short)allhits[ihit].Time<<" ";
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

  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"nht "<<nht<<" fit par "<<(int)clpar[0]<<"+/-"<<(int)intcpterr
    <<" "<<clpar[1]<<"+/-"<<slopeerr
    <<" clChisq "<<clChisq;
  
  }
/////////////////////////////////////////
  float ClusterCrawlerAlg::AngleFactor(float slope)
  {
    // returns an angle dependent cluster projection error factor for fitting
    // and hit finding
    
    float slp = fabs(slope);
    if(slp > 15) slp = 15;
    // return a value between 1 and 4
    float angfac = 1 + 0.03 * slp * slp;
    return angfac;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterChg(
    std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // Fits the charge of hits on the fcl2hits vector to a line, or simply
    // uses the average of 1 or 2 hits as determined by NHitsAve

    unsigned short ih0 = fcl2hits.size() - 1;
    
    if(pass >= fNumPass) {
      mf::LogError("ClusterCrawler")<<"FitClusterChg bad pass "<<pass;
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
      fAveChg = allhits[fcl2hits[ih0]].Charge;
      fChgSlp = 0.;
    } else if(fNHitsAve[pass] == 2) {
      // average the last two points if requested
      fAveChg = (allhits[fcl2hits[ih0]].Charge + 
                 allhits[fcl2hits[ih0 - 1]].Charge) / 2.;
      fChgSlp = 0.;
    } else if((unsigned short)fcl2hits.size() > fitLen){
      // do a real fit
      std::vector<float> xwir;
      std::vector<float> ychg;
      std::vector<float> ychgerr2;
      // origin of the fit
      unsigned short wire0 = allhits[fcl2hits[fcl2hits.size()-1]].WireNum;
      // find the mean and rms of the charge
      unsigned short npt = 0;
      unsigned short imlast = 0;
      float ave = 0.;
      float rms = 0.;
      // this loop intentionally ignores the Begin hit
      for(unsigned int ii = fcl2hits.size() - 1; ii > 0; --ii) {
        ++npt;
        float chg = allhits[fcl2hits[ii]].Charge;
        ave += chg;
        rms += chg * chg;
        if(npt == fitLen) {
          imlast = ii;
          break;
        }
      }
      float fnpt = npt;
      ave /= fnpt;
      rms = sqrt((rms - fnpt * ave * ave) / (fnpt - 1));
      float chgcut = ave + rms;
      for(unsigned int ii = fcl2hits.size() - 1; ii > imlast; --ii) {
        unsigned short wire = allhits[fcl2hits[ii]].WireNum;
        float chg = allhits[fcl2hits[ii]].Charge;
        if(chg > chgcut) continue;
        xwir.push_back((float)(wire - wire0));
        ychg.push_back(chg);
        ychgerr2.push_back(chg);
      }
      if(ychg.size() < 3) return;
      float intcpt; float slope; float intcpterr;
      float slopeerr; float chidof;
      fLinFitAlg.LinFit(xwir, ychg, ychgerr2, intcpt, slope, intcpterr, slopeerr, chidof);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"FitClusterChg wire "<<wire0
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
  void ClusterCrawlerAlg::AddLAHit(std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<VtxStore>& vtx,
    unsigned short kwire, bool& ChkCharge, bool& HitOK, bool& SigOK)
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
    unsigned short wire0 = allhits[lastClHit].WireNum;
    // the projected time of the cluster on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    float chgrat;
    
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"AddLAHit: wire "<<kwire<<" prtime "<<(int)prtime
    <<" max time diff "<<timeDiff;
    unsigned short imbest = 0;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(allhits[khit].InClus < 0) continue;
      chgrat = allhits[khit].Charge / allhits[lastClHit].Charge;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Chk W:T "<<kwire<<":"<<(short)allhits[khit].Time
    <<" Charge "<<(short)allhits[khit].Charge
    <<" chgrat "<<std::setw(8)<<std::setprecision(2)<<chgrat
    <<" InClus "<<allhits[khit].InClus
    <<" mult "<<allhits[khit].numHits
    <<" RMS "<<std::setprecision(2)<<allhits[khit].RMS
    <<" Chi2 "<<std::setw(8)<<std::setprecision(2)<<allhits[khit].ChiDOF
    <<" LoT "<<(int)allhits[khit].LoTime
    <<" HiT "<<(int)allhits[khit].HiTime;
      // projected time outside the Signal time window?
      if(prtime < allhits[khit].LoTime - 20) continue;
      if(prtime > allhits[khit].HiTime + 20) continue;
      SigOK = true;
      // hit used?
      if(allhits[khit].InClus > 0) continue;
      // ignore very low charge hits
      if(chgrat < 0.1) continue;
      // ignore very high charge hits for intermediate angle clusters
      if(fabs(clpar[1]) < 5 && chgrat > 2) continue;
      dtime = fabs(prtime - allhits[khit].Time);
      // Use a tighter requirement for crude hits that may span a very large
      // time range. Crude hits have a large chisq and/or large rms
      if((allhits[khit].ChiDOF > 500 || allhits[khit].RMS > 20) &&
          dtime > 20) continue;
      if(dtime < timeDiff) {
        HitOK = true;
        imbest = khit;
        timeDiff = dtime;
      }
    } // khit
    
  if(prt) {
    if(!HitOK) mf::LogVerbatim("ClusterCrawler")<<" no hit found ";
  }
    if(!HitOK) return;

  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Pick hit time "<<(int)allhits[imbest].Time
    <<" hit index "<<imbest;
    
    // merge the hits in a multiplet?
    bool doMerge = false;
    // assume there is no nearby hit
    short hnear = 0;
    if(allhits[imbest].numHits > 1) {
      doMerge = true;
      // don't merge if we are close to a vertex
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(vtx[ivx].CTP != clCTP) continue;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" close vtx chk W:T "<<vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time;
        if(abs(kwire - vtx[ivx].Wire) < 5 &&
           abs(allhits[imbest].Time - vtx[ivx].Time) < 20 ) {
  if(prt) mf::LogVerbatim("ClusterCrawler")
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
        float chicut = 6 * (1 - 1/(1 + fabs(clpar[1])));
        // look for a big separation between adjacent hits
        for(unsigned short jj = 0; jj < allhits[imbest].numHits; ++jj) {
          unsigned short jht = allhits[imbest].LoHitID + jj;
          // ignore obsolete hits
          if(allhits[jht].InClus < 0) continue;
          // count used hits
          if(allhits[jht].InClus > 0) ++nused;
          if(allhits[jht].InClus == 0) multipletChg += allhits[jht].Charge;
          // check the neighbor hit separation
          if(jj > 0) {
            // pick the larger RMS of the two hits
            float hitRMS = allhits[jht].RMS;
            if(allhits[jht - 1].RMS > hitRMS) hitRMS = allhits[jht-1].RMS;
            float tdiff = fabs(allhits[jht].Time - allhits[jht-1].Time) / hitRMS;
            if(tdiff > chicut) doMerge = false;
          } // jj > 0
        } // jj
  if(prt) {
    if(!doMerge) mf::LogVerbatim("ClusterCrawler")
      <<" Hits are well separated. Don't merge them";
  }
        if(doMerge && nused == 0) {
          // compare the charge with the last hit added?
          if(ChkCharge) {
            // there is a nearby hit
            hnear = 1;
            float chgrat = multipletChg / allhits[lastClHit].Charge;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" merge hits charge check "
    <<(int)multipletChg<<" Previous hit charge "<<(int)allhits[lastClHit].Charge;
            if(chgrat > 2) doMerge = false;
          }
        } // doMerge && nused == 0
      } // doMerge true
      if(doMerge) {
        // there is a nearby hit and it will be merged
        hnear = -1;
        MergeHits(allhits, imbest);
      }
    } // allhits[imbest].numHits > 1
    
    // attach to the cluster and fit
    fcl2hits.push_back(imbest);
    FitCluster(allhits);
    chifits.push_back(clChisq);
    hitnear.push_back(hnear);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" >>ADD W:T "<<kwire<<":"<<(int)allhits[imbest].Time
    <<std::setprecision(3)<<" clChisq "<<clChisq
    <<" charge "<<(int)allhits[imbest].Charge;
    // decide what to do with a bad fit
    if(clChisq > fChiCut[pass]) {
      fcl2hits.pop_back();
      chifits.pop_back();
      hitnear.pop_back();
      FitCluster(allhits);
      HitOK = false;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Bad fit. Removed hit. New clChisq "
    <<std::setprecision(3)<<clChisq
    <<" nhits "<<fcl2hits.size();
    }
      // stop tracking if previous chisq values were close to the cut.
      // this is an indicator that the track is wandering too much for this pass
      if(chifits.size() > 2 && 
        chifits[chifits.size()-2] > 0.8 * fChiCut[pass]) SigOK = false;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"  Set SigOK = "<<SigOK;

  } // AddLAHit


/////////////////////////////////////////
  void ClusterCrawlerAlg::AddHit(std::vector<CCHitFinderAlg::CCHit>& allhits,
    unsigned short kwire, bool& HitOK, bool& SigOK)
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
    unsigned short wire0 = allhits[lastClHit].WireNum;

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
    float prtimerr2 = fabs(kwire-wire0)*clparerr[1]*clparerr[1];
    // apply an angle dependent scale factor to the hit error
    float hiterr = AngleFactor(clpar[1]) * fHitErrFac * allhits[lastClHit].RMS;
    float err = sqrt(prtimerr2 + hiterr * hiterr);
    // Time window for accepting a hit.
    float prtimeLo = prtime - 4 * err;
    float prtimeHi = prtime + 4 * err;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"AddHit: wire "<<kwire
    <<" prtime Lo "<<(int)prtimeLo<<" Hi "<<(int)prtimeHi
    <<" fAveChg "<<(int)fAveChg;

    // loop through the hits
    int imbest = -1;
    float best = 9999.;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(allhits[khit].InClus < 0) continue;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Chk W:T "<<kwire<<":"<<(short)allhits[khit].Time
    <<" InClus "<<allhits[khit].InClus
    <<" mult "<<allhits[khit].numHits
    <<" RMS "<<std::setprecision(2)<<allhits[khit].RMS
    <<" Chi2 "<<std::setprecision(2)<<allhits[khit].ChiDOF
    <<" Charge "<<(int)allhits[khit].Charge
    <<" LoT "<<(int)allhits[khit].LoTime
    <<" HiT "<<(int)allhits[khit].HiTime;
      // check for signal
      if(prtime < allhits[khit].LoTime) continue;
      if(prtime > allhits[khit].HiTime) continue;
      SigOK = true;
      // check for good hit
      if(allhits[khit].Time < prtimeLo) continue;
      if(allhits[khit].Time > prtimeHi) continue;
      // hit used?
      if(allhits[khit].InClus > 0) continue;
      float dtime = fabs(allhits[khit].Time - prtime);
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

    if(imbest < 0) return;

    if(prt) mf::LogVerbatim("ClusterCrawler")
      <<" Best hit time "<<(int)allhits[imbest].Time;

    // merge hits in a doublet?
    // assume there is no nearby hit
    short hnear = 0;
    bool doMerge = false;
    if(fHitMergeChiCut > 0 && allhits[imbest].numHits == 2) {
      doMerge = true;
      // don't merge hits if near a vertex
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(abs(kwire - vtx[ivx].Wire) < 10 &&
           abs(allhits[imbest].Time - vtx[ivx].Time) < 20 ) doMerge = false;
      } // ivx
      if(doMerge) {
        // find the neighbor hit
        unsigned short imbestn = 0;
        if(allhits[imbest].LoHitID == (unsigned int)imbest) {
          imbestn = imbest + 1;
        } else {
          imbestn = imbest - 1;
        }
        // is the neighbor close and unused?
        float hitSep = fabs(allhits[imbest].Time - allhits[imbestn].Time);
        hitSep = hitSep / allhits[imbest].RMS;
        if(hitSep < fHitMergeChiCut && allhits[imbestn].InClus == 0) {
          // there is a nearby hit
          hnear = 1;
          // Is the charge of the doublet more similar to the charge of the
          // previously added hits than the single hit
          float totChg = allhits[imbest].Charge + allhits[imbestn].Charge;
          float lastHitChg = fAveChg;
          if(lastHitChg < 0) lastHitChg = allhits[lastClHit].Charge;
          // decide whether to merge
      if(prt) mf::LogVerbatim("ClusterCrawler")
        <<" merge hits charge check: totChg "<<totChg<<" lastHitChg "<<lastHitChg
        <<" hit chg "<<allhits[imbest].Charge;
          if(abs(totChg - lastHitChg) < abs(allhits[imbest].Charge - lastHitChg)) {
            // the total charge of both hits is a better match than the 
            // charge of the hit selected
            MergeHits(allhits, imbest);
            // the nearby hit was merged
            hnear = -1;
      if(prt) mf::LogVerbatim("ClusterCrawler")
        <<" Merging hits "<<imbest<<" and "<<imbestn
        <<" New Time "<<allhits[imbest].Time
        <<" New Chg "<<allhits[imbest].Charge;
          } // merge hit
        } // bestn < fHitMergeChiCut
      } // doMerge true
    } // fHitMergeChiCut > 0 && allhits[imbest].numHits == 2

    // Make a charge similarity cut if the average charge is defined
    bool fitChg = true;
    if(fAveChg > 0.) {

      float chgrat = (allhits[imbest].Charge - fAveChg) / fAveChg;
    if(prt) mf::LogVerbatim("ClusterCrawler")
      <<" Chgrat "<<std::setprecision(2)<<chgrat;

      // charge is way too high?
      if(chgrat > 2 * fChgCut[pass]) {
        if(prt) mf::LogVerbatim("ClusterCrawler")<<" fails high charge cut";
        return;
      }

      // Determine if the last hit added was a large (low) charge hit
      // This will be used to prevent adding large (low) charge hits on two
      // consecutive fits. This cut is only applied to hits on adjacent wires
      float bigchgcut = 1.5 * fChgCut[pass];
      bool lasthitbig = false;
      bool lasthitlow = false;
      if(abs(wire0 - kwire) == 1) {
        float lastchgrat = (allhits[lastClHit].Charge - fAveChg) / fAveChg;
        lasthitbig = ( lastchgrat > bigchgcut);
        lasthitlow = ( lastchgrat < -fChgCut[pass]);
      }
      
      // the last hit added was low charge and this one is as well
      if(lasthitlow && chgrat < -fChgCut[pass]) {
        if(prt) mf::LogVerbatim("ClusterCrawler")<<" fails low charge cut. Stop crawling.";
        SigOK = false;
        return;
      } // lasthitlow
    
      // the last hit was high charge and this one is also
      if(lasthitbig && chgrat > fChgCut[pass]) {
        if(prt) mf::LogVerbatim("ClusterCrawler")<<" fails 2nd high charge cut";
        return;
      } // lasthitbig

    
      // require that large charge hits have a very good projection error
      if(chgrat > fChgCut[pass]) {
        if(best > 2 * err) {
          if(prt) mf::LogVerbatim("ClusterCrawler")<<" high charge && bad dT= "
            <<best<<" err= "<<err;
          return;
        }
      } // chgrat > fChgCut[pass]

      // decide whether to fit the charge
      fitChg = (chgrat < fabs(fChgCut[pass]) );
    } // fAveChg > 0
    
    // we now have a hit that meets all the criteria. Fit it
    fcl2hits.push_back(imbest);
    FitCluster(allhits);
    chifits.push_back(clChisq);
    hitnear.push_back(hnear);
    // nearby hit check
    ChkClusterNearbyHits(allhits, prt);
    HitOK = true;

  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
      <<" >>ADD W:T "<<kwire<<":"<<(short)allhits[imbest].Time<<" dT "<<best
      <<std::setprecision(2)<<" Chisq "<<clChisq
      <<" Chg "<<(int)allhits[imbest].Charge<<" AveChg "<<(int)fAveChg
      <<" hitnear "<<allhits[imbest].numHits
      <<" fcl2hits size "<<fcl2hits.size();
  }
    if(!fitChg) return;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" Fit charge ";
    FitClusterChg(allhits);
  } // AddHit


//////////////////////////////////////
    void ClusterCrawlerAlg::ChkClusterNearbyHits(
        std::vector<CCHitFinderAlg::CCHit>& allhits, bool prt)
    {
      // analyze the hitnear vector
      //  0 = no nearby hit exists
      //  1 = a nearby hit exists but was not merged
      // -1 = a nearby hit was merged
      
      if(fHitMergeChiCut <= 0) return;
      
      if(hitnear.size() != fcl2hits.size()) {
        mf::LogWarning("ClusterCrawler")<<"Coding error: hitnear size != fcl2hits";
        return;
      }
      
      // Analyze the last 6 hits added but don't consider the first few hits
      if(hitnear.size() < 12) return;
      
      unsigned short ii, indx;
      unsigned short merged = 0;
      unsigned short notmerged = 0;
      for(ii = 0; ii < 6; ++ii) {
        indx = hitnear.size() - 1 - ii;
        if(hitnear[indx] > 0) ++notmerged;
        if(hitnear[indx] < 0) ++merged;
      }
      
//  if(prt) mf::LogVerbatim("ClusterCrawler")
//    <<"ChkClusterNearbyHits: nearby hits merged "<<merged
//    <<" not merged "<<notmerged;

      if(notmerged < 2) return;
      
      // a number of nearby hits were not merged while crawling, so the 
      // average charge is probably wrong. Look at the last 6 hits added
      // and merge them if they are close
      unsigned short iht, oht;
      float hitSep, chgRat;
      for(ii = 0; ii < 6; ++ii) {
        indx = fcl2hits.size() - 1 - ii;
        iht = fcl2hits[indx];
        if(allhits[iht].numHits == 2) {
          // hit doublet. Get the index of the other hit
          if(allhits[iht].LoHitID == iht) {
            oht = iht + 1;
          } else {
            oht = iht - 1;
          } // allhits[iht].LoHitID == iht
          hitSep = fabs(allhits[iht].Time - allhits[oht].Time);
          hitSep /= allhits[iht].RMS;
          if(hitSep < fHitMergeChiCut && allhits[oht].InClus == 0) {
            // check charge assuming the hits are merged
            chgRat = allhits[iht].Charge + allhits[oht].Charge;
            chgRat /= allhits[iht].Charge;
            // merge em
            if(chgRat < 4) {
      if(prt) mf::LogVerbatim("ClusterCrawler")<<"Merging hit doublet "<<iht;
              MergeHits(allhits, iht);
              hitnear[indx] = -1;
            }
          } // hitSep OK and not in a cluster
        } // hit doublet
      } // ii
      
      // now re-fit
      FitCluster(allhits);
      FitClusterChg(allhits);

      if(prt) mf::LogVerbatim("ClusterCrawler")
        <<"ChkClusterNearbyHits refit cluster. fAveChg= "<<fAveChg;
      
    } // ChkClusterHitNear

//////////////////////////////////////
    void ClusterCrawlerAlg::FitVtx(std::vector<ClusterStore>& tcl,
        std::vector<VtxStore>& vtx, unsigned short iv, float& ChiDOF)
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
      if(fabs(vwire - vtx[iv].Wire) > 2) return;
      if(fabs(vtime - vtx[iv].Time) > 10) return;
      vtx[iv].Wire = vwire;
      vtx[iv].Time = vtime;
    } // FitVtx

//////////////////////////////////////
    void ClusterCrawlerAlg::Vtx3ClusterMatch(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
       std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
       std::vector<Vtx3Store>& vtx3, unsigned int cstat, unsigned int tpc)
      {
        // Look for clusters that end/begin near the expected wire/time
        // for incomplete 3D vertices
        if(vtx3.size() == 0) return;
        
        unsigned short thePlane, theWire;
        float dw, theTime, dt;
        short dwb, dwe;

        for(unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
          // A complete 3D vertex with matching 2D vertices in all planes?
          if(vtx3[ivx].Wire < 0) continue;
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
            dwb = abs(theWire - tcl[icl].BeginWir);
            dwe = abs(theWire - tcl[icl].EndWir);
            if(dwb < dwe && dwb < 5) {
              // cluster begin is closer
              if(theWire < tcl[icl].BeginWir) continue;
              dw = theWire - tcl[icl].BeginWir;
              dt = tcl[icl].BeginTim + tcl[icl].BeginSlp * dw - theTime;
              if(fabs(dt) > 10) continue;
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
              if(fabs(dt) > 10) continue;
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
    void ClusterCrawlerAlg::Vtx3ClusterSplit(std::vector<CCHitFinderAlg::CCHit>& allhits,
       std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
       std::vector<Vtx3Store>& vtx3, unsigned int cstat, unsigned int tpc)
      {
        // Try to split clusters in a view in which there is no 2D vertex
        // assigned to a 3D vertex
        if(vtx3.size() == 0) return;
        
      vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);
        
        unsigned short lastplane = 5, kcl, kclID;
        float dth, theTime;
        unsigned short thePlane, theWire, plane;
        unsigned short loWire, hiWire;

        for(unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
          // Complete 3D vertex with matching 2D vertices in all planes?
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Vtx3ClusterSplit ivx "<<ivx
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
            (int)thePlane, (int)tpc, (int)cstat);
          // get the hit range if necessary
          if(thePlane != lastplane) {
            clCTP = EncodeCTP(cstat, tpc, thePlane);
            GetHitRange(allhits, clCTP, WireHitRange, fFirstWire, fLastWire);
            lastplane = thePlane;
          }
          // make a list of clusters that have hits near this point on nearby wires xxx
          std::vector<short> clIDs;
          if(theWire > fFirstWire + 5) { loWire = theWire - 5; } else { loWire = fFirstWire; }
          if(theWire < fLastWire  - 5) { hiWire = theWire + 5; } else { hiWire = fLastWire; }
    if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"3DVtx "<<ivx
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
              if(allhits[khit].Charge < 0) continue;
              // ignore un-assigned hits
              if(allhits[khit].InClus <= 0) continue;
              if((unsigned short)allhits[khit].InClus > tcl.size() + 1) {
                mf::LogError("ClusterCrawler")<<"Invalid hit InClus. "<<khit
                  <<" "<<allhits[khit].InClus;
                continue;
              }
              // check an expanded time range
              if(theTime < allhits[khit].LoTime - 10) continue;
              if(theTime > allhits[khit].HiTime + 10) continue;
              kclID = allhits[khit].InClus;
              kcl = kclID - 1;
              // ignore obsolete clusters
              if(tcl[kcl].ID < 0) continue;
              
              // put the cluster in the list if it's not there already
    if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Bingo "<<ivx<<" plane "<<thePlane
      <<" wire "<<wire<<" hit "<<allhits[khit].WireNum<<":"<<(int)allhits[khit].Time
        <<" inClus "<<allhits[khit].InClus
        <<" P:W:T "<<thePlane<<":"<<tcl[kcl].BeginWir<<":"<<(int)tcl[kcl].BeginTim;
              if(std::find(clIDs.begin(), clIDs.end(), kclID) == clIDs.end()) {
                // ignore long straight clusters
                if(tcl[kcl].tclhits.size() > 100 ) {
                  dth = tcl[kcl].BeginAng - tcl[kcl].EndAng;
    if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Long straight check: nhits "
      <<tcl[kcl].tclhits.size()<<" dth "<<dth;
                  if(fabs(dth) < 0.05) continue;
                } // tcl[kcl].tclhits.size() > 100
                clIDs.push_back(kclID);
              } // std::find
            } // khit
          } // wire
          if(clIDs.size() == 0) continue;
    if(vtxprt) {
      for(unsigned int ii = 0; ii < clIDs.size(); ++ii) {
        mf::LogVerbatim("ClusterCrawler")<<" clIDs "<<clIDs[ii];
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
              if(allhits[iht].WireNum < theWire) {
                nhitfit = 3;
                if(jj > 3) nhitfit = -3;
                FitClusterMid(allhits, tcl, icl, iht, nhitfit);
                dth = clpar[0] + (theWire - allhits[iht].WireNum) * clpar[1] - theTime;
    if(vtxprt) {
      for(unsigned short kk = 0; kk < clIDs.size(); ++kk) {
        mf::LogVerbatim("ClusterCrawler")<<" cls "<<icl<<" dt "<<dth;
      }
    }
                if(fabs(dth) > 5) clIDs[ii] = -1;
                didit = true;
                break;
              } // allhits[iht].WireNum < theWire
              if(didit) break;
            } // jj
            if(didit) break;
          } // ii
    if(vtxprt) {
      mf::LogVerbatim("ClusterCrawler")<<"clIDs size after fit "<<clIDs.size();
      for(ii = 0; ii < clIDs.size(); ++ii) {
        mf::LogVerbatim("ClusterCrawler")<<" clIDs "<<clIDs[ii];
      }
    }

          // see if any candidates remain
          unsigned short nok = 0;
          for(ii = 0; ii < clIDs.size(); ++ii) if(clIDs[ii] >= 0) ++nok;
          if(nok == 0) continue;

          // make a new 2D vertex
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Make new 2D vtx in plane "<<thePlane
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
              if(allhits[iht].WireNum < theWire) {
                pos = jj;
                break;
              }
            } // jj
            if(pos == 0) {
              // vertex is DS of the cluster Begin
              tcl[icl].BeginVtx = ivnew;
              ++nvcl;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Attach to Begin "<<icl;
            } else if(pos > tcl[icl].tclhits.size()) {
              // vertex is US of the cluster Eend
              tcl[icl].EndVtx = ivnew;
              ++nvcl;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Attach to End "<<icl;
            } else {
              // vertex is in the middle of the cluster
              SplitCluster(allhits, tcl, icl, pos, ivnew);
              tcl[icl].ProcCode += 10000;
              tcl[tcl.size()-1].ProcCode += 10000;
              nvcl += 2;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
    <<"Split cluster "<<clIDs[ii]<<" at pos "<<pos;
            } // pos check
          } // ii
          // try to attach other clusters
//          VertexCluster(tcl, vtx, ivnew);
          // Fit the vertex position
          float chisq = 0;
          FitVtx(tcl, vtx, ivnew, chisq);
          // mark the 3D vertex as complete
          vtx3[ivx].Wire = -1;
        } // ivx
        
      } // Vtx3ClusterSplit


//////////////////////////////////////
    void ClusterCrawlerAlg::VtxMatch(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      std::vector<Vtx3Store>& vtx3, unsigned int cstat, unsigned int tpc)
    {
      // Create 3D vertices from 2D vertices. 3D vertices that are matched
      // in all three planes have Ptr2D >= 0 for all planes
      
      unsigned int nPln = geom->Cryostat(cstat).TPC(tpc).Nplanes();
      if(nPln != 3) return;
      
      vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);
      
      // wire spacing in cm
      float wirePitch = geom->WirePitch(0, 1, 0, tpc, cstat);
      
      // Y,Z limits of the detector
      double YHi = geom->DetHalfHeight(tpc, cstat);
      double YLo = -YHi;
      double ZLo = 0.;
      double ZHi = geom->DetLength(tpc, cstat);
      
      // create a vector of vertex indices in each plane
      std::vector<std::vector<unsigned short>> vIndex;
      std::vector<unsigned short> temp;
      for(unsigned short ipl = 0; ipl < 3; ++ipl) {
        temp.clear();
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].Wght < 0) continue;
          geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
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
/*
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
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
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
    <<"jpl "<<jpl<<" jvx "<<jvx<<" W:T "<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time
    <<" jX "<<jX;
              if(fabs(jX - iX) > fVertex3DCut) continue;
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"2DMatchX "<<jX - iX;
              jWire = vtx[jvx].Wire;
              geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
              if(y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
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

  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
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
                float dr = 0.5 * sqrt(dW * dW + dX * dX);
                if(dr < best) {
                  best = dr;
                  xbest = (kvxX + 2 * kX) / 3;
                  ybest = y;
                  zbest = z;
                  t2dIndex[ipl] = ivx;
                  t2dIndex[jpl] = jvx;
                  t2dIndex[kpl] = kvx;
                }

  if(vtxprt) mf::LogVerbatim("ClusterCrawler")
    <<" kvx "<<kvx<<" kpl "<<kpl
    <<" wire "<<(int)vtx[kvx].Wire<<" kTime "<<(int)vtx[kvx].Time
    <<" dr "<<dr;

              } // kk
            } // jj
          } // jpl
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"3DMatch best "<<best;
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

 if(vtxprt) mf::LogVerbatim("ClusterCrawler")
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
      mf::LogVerbatim("ClusterCrawler")
        <<"vtx3 "<<it<<" Ptr2D "<<vtx3[it].Ptr2D[0]<<
        " "<<vtx3[it].Ptr2D[1]<<" "<<vtx3[it].Ptr2D[2]
        <<" wire "<<vtx3[it].Wire;
    }
  }

    } // VtxMatch

//////////////////////////////////
    void ClusterCrawlerAlg::GetHitRange(std::vector<CCHitFinderAlg::CCHit>& allhits,
      CTP_t CTP, 
      std::vector< std::pair<short, short> >& WireHitRange,
      unsigned short& firstwire, unsigned short& lastwire)
    {
      art::ServiceHandle<geo::Geometry> geom;
      // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
      bool first = true;
      lastwire = 0;
      unsigned short firsthit = 0;
      geo::PlaneID planeID = DecodeCTP(CTP);
      unsigned short lasthit = 0;
      // find the first and last wire with a hit
      for(unsigned short hit = 0; hit < allhits.size(); ++hit) {
        if(allhits[hit].WirID.Plane != planeID.Plane) continue;
        if(allhits[hit].WirID.TPC != planeID.TPC) continue;
        if(allhits[hit].WirID.Cryostat != planeID.Cryostat) continue;
        unsigned short theWireNum = allhits[hit].WireNum;
        if(theWireNum != allhits[hit].WireNum)
          mf::LogError("ClusterCrawler")<<"GetHitRange WireNum screwup ";
        if(first) {
          firsthit = hit;
          firstwire = theWireNum;
          first = false;
        }
        lastwire = theWireNum;
        lasthit = hit;
      } //hit

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
        uint32_t chan = geom->PlaneWireToChannel
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
        CCHitFinderAlg::CCHit& theHit = allhits[hit];
        unsigned short thiswire = theHit.WireNum;
        if(thiswire > lastwire) {
          unsigned short index = lastwire - firstwire;
          short itmp1 = lastfirsthit;
          short itmp2 = thishit;
          WireHitRange[index] = std::make_pair(itmp1,itmp2);
          lastwire = thiswire;
          lastfirsthit = thishit;
        } else if(thiswire < lastwire) {
          mf::LogError("ClusterCrawler")<<"ERROR: Hits not sorted!!";
          return;
        }
        ++thishit;
      } //hit
      // define for the last wire
      unsigned short index = lastwire - firstwire;
      short itmp1 = lastfirsthit;
      short itmp2 = thishit;
      WireHitRange[index] = std::make_pair(itmp1,itmp2);
    } // GetHitRange

/////////////////////////////////////////
    void ClusterCrawlerAlg::SortByLength(std::vector<ClusterStore>& tcl,
        CTP_t inCTP, std::map<unsigned short, unsigned short>& sortindex)
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

} // namespace cluster
