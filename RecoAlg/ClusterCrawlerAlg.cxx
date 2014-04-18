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
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "RecoAlg/CCHitFinderAlg.h"

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
    fWidCut             = pset.get< std::vector<float> >("WidCut");
    fChiCut             = pset.get< std::vector<float> >("ChiCut");
    fMaxWirSkip         = pset.get< std::vector<unsigned short> >("MaxWirSkip");
    fMinWirAfterSkip    = pset.get< std::vector<unsigned short> >("MinWirAfterSkip");
    fKinkChiRat         = pset.get< std::vector<float> >("KinkChiRat");
    fKinkAngCut         = pset.get< std::vector<float> >("KinkAngCut");
    fDoMerge            = pset.get< std::vector<bool>  >("DoMerge");
    fTimeDelta          = pset.get< std::vector<float> >("TimeDelta");
    fFindVertices       = pset.get< std::vector<bool>  >("FindVertices");

    fHitErrFac          = pset.get<             float  >("HitErrFac");
    fLAClusAngleCut     = pset.get<             float  >("LAClusAngleCut");
    fClHitMergeChiCut   = pset.get<             float  >("ClHitMergeChiCut");
    fClGhostHitFrac     = pset.get<             float  >("ClGhostHitFrac");
    fAllowNoHitWire     = pset.get<    unsigned short  >("AllowNoHitWire");
    fDebugPlane         = pset.get<             short  >("DebugPlane");
    fDebugWire          = pset.get<             short  >("DebugWire");
    fDebugHit           = pset.get<             short  >("DebugHit");

    // some error checking
    bool badinput = false;
    if(fNumPass > fMaxHitsFit.size()) badinput = true;
    if(fNumPass > fMinHits.size()) badinput = true;
    if(fNumPass > fNHitsAve.size()) badinput = true;
    if(fNumPass > fChgCut.size()) badinput = true;
    if(fNumPass > fWidCut.size()) badinput = true;
    if(fNumPass > fChiCut.size()) badinput = true;
    if(fNumPass > fMaxWirSkip.size()) badinput = true;
    if(fNumPass > fMinWirAfterSkip.size()) badinput = true;
    if(fNumPass > fKinkChiRat.size()) badinput = true;
    if(fNumPass > fKinkAngCut.size()) badinput = true;
    if(fNumPass > fDoMerge.size()) badinput = true;
    if(fNumPass > fTimeDelta.size()) badinput = true;

    if(badinput) throw cet::exception("ClusterCrawler")<<"Bad input from fcl file ";

  }

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
    clStopCode = 0; clProcCode = 0; clAssn = 0;        fFirstWire = 0;
    fLastWire = 0; fAveChg = 0.; fChgSlp = 0.; fAveWid = 0.; pass = 0;
    fScaleF = 0; tcl.clear(); vtx.clear(); WireHitRange.clear();
  }


  void ClusterCrawlerAlg::RunCrawler(std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // Run the ClusterCrawler algorithm - creating seed clusters and crawling upstream.

    CrawlInit();
    
    if(allhits.size() < 3) return;

    
    for(cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      for(tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
        for(plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
          WireHitRange.clear();
          // define a code to ensure clusters are compared within the same plane
          clCTP = EncodeCTP(cstat, tpc, plane);
          // fill the WireHitRange vector with first/last hit on each wire
          // dead wires and wires with no hits are flagged < 0
          GetHitRange(allhits, clCTP, WireHitRange, fFirstWire, fLastWire);
          fFirstHit = WireHitRange[0].first;
          // get the scale factor to convert dTick/dWire to dX/dU. This is used
          // to make the kink and merging cuts
          art::Ptr<recob::Wire> theWire = allhits[fFirstHit].Wire;
          uint32_t channel = theWire->RawDigit()->Channel();
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
          ClusterLoop(allhits);
        } // plane
      } // tpc
    } // cstat
    
    WireHitRange.clear(); 
    fcl2hits.clear();
    chifits.clear();
    
  } // RunCrawler
    
////////////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterLoop(std::vector<CCHitFinderAlg::CCHit>& allhits) {
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
            // skip multiple hits except on the last pass
            if(pass < fNumPass - 1 && allhits[ihit].numHits > 1) continue;
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
              fAveWid = -1.;
              fAveChg = -1.;
              clBeginChg = -1.;
              clEndChg = -1.;
              clStopCode = 0;
              clProcCode = pass;
              clAssn = -1; 
              fcl2hits.push_back(ihit);
              chifits.push_back(0.);
              fcl2hits.push_back(jhit);
              chifits.push_back(0.);
              clpar[0] = allhits[jhit].Time;
              clpar[1] = (allhits[ihit].Time - allhits[jhit].Time) / (iwire - jwire);
              clChisq = 0;
              // now look for hits to add on the intervening wires
              bool SigOK = false;
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
              // Hits will be added in the proper order by cl2Crawl
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
              // save the cluster begin info
              clBeginWir = iwire;
              clBeginTim = allhits[ihit].Time;
              clBeginSlp = clpar[1];
              clBeginSlpErr = clparerr[1];
              // decide whether to crawl a large angle cluster. Requirements are:
              // 1) the user has set the LACluster angle cut > 0, AND
              // 2) the cluster slope exceeds the cut, AND
              // 3) this is the last pass
              // Note that if condition 1 is met, normal cluster crawling is done
              // only if the slope is less than the cut
              fCrawlLACluster = false;
              if(fLAClusSlopeCut > 0) {
                // LA cluster crawling requested
                if(fabs(clBeginSlp) > fLAClusSlopeCut) {
                  // skip if this is not the last pass
                  if(pass != fNumPass - 1) continue;
                  // Crawl with LA cluster code. Set the flag to ignore the charge ratio cut
                  fCrawlLACluster = true;
                  fAveWid = -1.;
                  LACrawlUS(allhits, vtx);
                } else {
                  CrawlUS(allhits, vtx);
                } // fabs(clBeginSlp) > fLAClusAngleCut
              } else {
                // allow clusters of any angle
                CrawlUS(allhits, vtx);
              } // fLAClusSlopeCut > 0
              // do a quality check. fcl2hits size set 0 if bad cluster
              QACheck(allhits);
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
              if(pass < fNumPass - 2) {
                // Is it long enough for the next pass?
                if(fcl2hits.size() >= fMinHits[pass+1]) {
                  clEndSlp = clpar[1]; // save the slope at the end
                  clEndSlpErr = clparerr[1];
                  // set a special code 
                  clProcCode += 2000;
                  TmpStore(allhits, tcl);
                  ClusterAdded = true;
                  nHitsUsed += fcl2hits.size();
                  AllDone = (nHitsUsed == allhits.size());
                  break;
                } // long enough
              } // pass < fNumPass
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
      
      // merge unused hits on cluster hit multiplets
      MergeClusterHits(allhits, tcl);
      
      // merge clusters that share a significant fraction of hits in 
      // the same hit multiplet
      MergeGhostClusters(allhits, tcl);

      // split clusters using vertices
      VtxClusterSplit(allhits, tcl, vtx);

  if(fDebugPlane == (short)plane) {
    mf::LogVerbatim("ClusterCrawler")<<"Clustering done in plane "<<plane;
    PrintClusters(allhits, tcl, vtx);
  }
    
  } // ClusterLoop

//////////////////////////////////////////
  void ClusterCrawlerAlg::QACheck(
      std::vector<CCHitFinderAlg::CCHit>& allhits)
  {

    if(fcl2hits.size() == 0) return;

    // ignore short large angle clusters that have multiplicity 1 hits
    if(fcl2hits.size() < 5 && fabs(clpar[1]) > 3.) {
      for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
        if(allhits[fcl2hits[ii]].numHits == 1) {
          fcl2hits.clear();
          return;
        }
      } // ii
    } // fcl2hits.size() < 5
    // Check the fraction of wires that have hits
    unsigned short iht = fcl2hits[fcl2hits.size() - 1];
    clEndWir = allhits[iht].WireNum;
    float hitFrac = (float)fcl2hits.size() / (float)(clBeginWir - clEndWir + 1);
    if(hitFrac < 0.7) {
      fcl2hits.clear();
      return;
    }
  } // cl2QACheck


//////////////////////////////////////////
  bool ClusterCrawlerAlg::CrawlVtxChk(
    std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<VtxStore>& vtx, unsigned short kwire)
  {
    
    // returns true if the cluster is near a vertex on wire kwire
    if(vtx.size() == 0) return false;
    float wire0 = allhits[fcl2hits[fcl2hits.size() - 1]].WireNum;
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      if(vtx[iv].Wire != kwire) continue;
      float prtime = fabs(clpar[0]+((float)kwire-wire0)*clpar[1] - vtx[iv].Time);
      if(prtime < fTimeDelta[pass]) return true;
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

      if(vtx.size() == 0) return;
      unsigned short tclsize = tcl.size();
      if(tclsize < 2) return;

      for(unsigned short icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore short clusters
        if(tcl[icl].tclhits.size() < 5) continue;
        // check vertices
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].CTP != clCTP) continue;
          // abandoned vertex?
          if(vtx[ivx].Wght <= 0) continue;
          // already assigned to this vertex?
          if(tcl[icl].BeginVtx == ivx) continue;
          if(tcl[icl].EndVtx == ivx) continue;
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
//  if(vtxprt) std::cout<<"Chk vertex "<<ivx<<std::endl;
          short ihvx = -1;
          // nSplit is the index of the hit in the cluster where we will
          // split it if all requirements are met
          unsigned short nSplit = 0;
          unsigned short nLop = 0;
          for(unsigned short ii = tcl[icl].tclhits.size()-1; ii > 0 ; --ii) {
            unsigned short iht = tcl[icl].tclhits[ii];
            ++nLop;
            if(allhits[iht].WireNum >= vtx[ivx].Wire) {
              nSplit = ii;
//  if(vtxprt) std::cout<<"Split at wire "<<allhits[iht].WireNum<<" nSplit "<<nSplit<<std::endl;
              ihvx = iht;
              break;
            }
          } // ii
          // found the wire. Now make a rough time cut
          if(ihvx > 0) {
            short dtime = abs(allhits[ihvx].Time - vtx[ivx].Time);
//  if(vtxprt) std::cout<<"dtime "<<dtime<<std::endl;
            if(dtime > 5) continue;
          }
          // check the angle between the crossing cluster icl and the
          // clusters that comprise the vertex. 
          // First decide which end of cluster icl to use to define the angle
          float iclAng = 0.;
          if(nSplit > tcl[icl].tclhits.size() / 2) {
            iclAng = atan(fScaleF * tcl[icl].EndSlp);
          } else {
            iclAng = atan(fScaleF * tcl[icl].BeginSlp);
          }
//  if(vtxprt) std::cout<<"iclAng "<<iclAng<<std::endl;
          // check angle wrt the the vertex clusters
          bool angOK = false;
          for(unsigned short jcl = 0; jcl < tclsize; ++jcl) {
            if(tcl[jcl].ID < 0) continue;
            if(tcl[jcl].CTP != clCTP) continue;
            if(tcl[jcl].BeginVtx == ivx) {
              float clAng =  atan(fScaleF * tcl[jcl].BeginSlp);
              if(fabs(clAng - iclAng) > 0.4) {
                // large angle difference. Set the flag
                angOK = true;
                break;
              }
            } // tcl[jcl].BeginVtx == ivx
            if(tcl[jcl].EndVtx == ivx) {
              float clAng =  atan(fScaleF * tcl[jcl].EndSlp);
              if(fabs(clAng - iclAng) > 0.4) {
                // large angle difference. Set the flag
                angOK = true;
                break;
              }
            } // tcl[jcl].EndVtx == ivx
          } // jcl
          // time to split or chop
          if(angOK) {
//  if(vtxprt) std::cout<<"Split/Chop "<<nLop<<std::endl;
            if(nLop < 3) {
              // lop off hits at the US end
              // Put the cluster in the local arrays
              TmpGet(allhits, tcl, icl);
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
              tcl[icl].ProcCode += 1000;
              // correct the split position
              ++nSplit;
              SplitCluster(allhits, tcl, icl, nSplit, ivx);
            }
            break;
          } // angOK
        } // ivx
      } // icl

    } // cl2VtxClusterSplit


//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeGhostClusters(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl)
    {
      // Merge clusters if they share a fraction of hits in the same hit multiplet

      if(fClGhostHitFrac < 0.) return;
      
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // vector of cluster ID's with "shared" hits and the number of shared hits
        std::vector<std::pair< unsigned short, unsigned short> > oClus;
        for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
          unsigned short iht = tcl[icl].tclhits[ii];
  std::cout<<" W:T "<<allhits[iht].WireNum<<":"<<(int)allhits[iht].Time
    <<" mult "<<allhits[iht].numHits<<std::endl;
          if(allhits[iht].numHits == 1) continue;
          for(unsigned short jj = 0; jj < allhits[iht].numHits; ++jj) {
            unsigned short jht = allhits[iht].LoHitID + jj;
            // unused or obsolete hit
            if(allhits[jht].InClus <= 0) continue;
            // hit used in cluster icl
            if(allhits[jht].InClus == tcl[icl].ID) continue;
            // get the cluster ID
            unsigned short jclID = allhits[jht].InClus;
  std::cout<<" jclID "<<jclID<<std::endl;
            // find this cluster ID in oClus
            bool found = false;
            for(unsigned short kk = 0; kk < oClus.size(); ++kk) {
              if(oClus[kk].first == jclID) {
                // 
                oClus[kk].second += 1;
                found = true;
                break;
              } // oClus[kk].first == tcl[jcl].ID
            } // kk
            // not found? add a new one
            if(!found) {
              oClus.push_back(std::make_pair(jclID, 1));
            } // !found
          } // jj
        } // ii
        // find the cluster that shares the most hits in the multiplet with
        // cluster icl
        unsigned short big = 0;
        unsigned short imbig = 0;
        for(unsigned short ii = 0; ii < oClus.size(); ++ii) {
  std::cout<<"pair "<<ii<<" "<<oClus[ii].first<<" "<<oClus[ii].second<<std::endl;
          if(oClus[ii].second > big) {
            big = oClus[ii].second;
            imbig = oClus[ii].first;
          }
        } // ii
  std::cout<<"Chk "<<tcl[icl].ID<<" size "<<tcl[icl].tclhits.size()
    <<" big "<<big<<" imbig "<<imbig<<std::endl;
      } // icl

    } // MergeGhostClusters

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeClusterHits(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl)
    {
      // merge unused cluster hits in a multiplet into one hit. This
      // is called after clustering is completed in a plane
      
      if(fClHitMergeChiCut < 0.) return;
      
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        unsigned short nMerged = 0;
        for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
          unsigned short theHit = tcl[icl].tclhits[ii];
  if(allhits[theHit].InClus != tcl[icl].ID) {
    std::cout<<"MergeClusterHits bad hit assignment "<<allhits[theHit].InClus<<" "
      <<tcl[icl].ID<<std::endl;
  }
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"MergeClusterHits: cluster "<<tcl[icl].ID<<" check hit "
    <<allhits[theHit].WireNum<<":"<<(int)allhits[theHit].Time;
          // check for hit multiplet
          if(allhits[theHit].numHits > 1) {
            // count the number of unused hits
            unsigned short nUnused = 0;
            unsigned short nClose = 0;
            for(unsigned short jj = 0; jj < allhits[theHit].numHits; ++jj) {
              unsigned short jht = allhits[theHit].LoHitID + jj;
              // hit used in a cluster?
              if(allhits[jht].InClus != 0) continue;
              ++nUnused;
              // count the number of close unused hits
              float hitSep = fabs(allhits[jht].Time - allhits[theHit].Time);
              hitSep = hitSep / allhits[theHit].RMS;
              if(hitSep < fClHitMergeChiCut) ++nClose;
            } // jj
            // all of the hits in the multiplet used?
            if(nUnused == 0) continue;
            // require all of the hits be close (Is this necessary????)
            if(nUnused != nClose) continue;
            MergeHits(allhits, theHit);
            ++nMerged;
          } // allhits[iht].numHits > 1
        } // ii (cluster hit)
        if(nMerged == 0) continue;
        // refit the cluster?
        // set the processor code
        tcl[icl].ProcCode += 3000;
      } // icl
    } // MergeClusterHits

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeHits(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short theHit)
    {
      // Merge all of the unused separate hits in the multiplet of which 
      // theHit is a member into one hit (= theHit). Set the charge of the 
      // merged hits < 0 to indicate they are obsolete. Hits in the multiplet
      // that are associated with an existing cluster are not affected

      if(theHit > allhits.size() - 1) {
        mf::LogError("ClusterCrawler")<<"Bad theHit";
        return;
      }
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Inside MergeHits "
    <<allhits[theHit].WireNum<<":"<<(int)allhits[theHit].Time;

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
      unsigned short nHitMult = allhits[theHit].numHits;
      // find the time range for the hit multiplet
      for(unsigned short jj = 0; jj < nHitMult; ++jj) {
        unsigned short jht = allhits[theHit].LoHitID + jj;

  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<"MergeHits P:W:T "<<plane<<":"<<allhits[jht].WireNum<<":"<<(int)allhits[jht].Time
    <<" Amp "<<allhits[jht].Amplitude
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
      } // jj
      // all hits in the multiplet used?
      if(nGaus == 0) return;
      if(loTime < 0) loTime = 0;
      ++hiTime;
      // define a signal shape
      std::vector<double> signal;
      // fill the signal with zeros
      for(unsigned short time = loTime; time < hiTime; ++time) signal.push_back(0.);
      // now add the Gaussians for each hit
      double chgsum = 0.;
      for(unsigned short jj = 0; jj < nHitMult; ++jj) {
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
    <<"theHit "<<allhits[theHit].WireNum<<":"<<(int)aveTime<<" RMS "<<allhits[theHit].RMS
    <<" chg "<<(int)chgsum<<" Amp "<<allhits[theHit].Amplitude;
  }

    } // MergeHits

/////////////////////////////////////////
    void ClusterCrawlerAlg::FindVertices(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {
      // try to make 2D vertices
      
      if(tcl.size() < 2) return;

      // form vertices starting with the longest
      std::map<unsigned short, unsigned short> sortindex;
      cl2SortByLength(tcl,sortindex);
      
      float nwires = fNumWires;
      float maxtime = fMaxTime;

  vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")<<"FindVertices plane "<<plane<<" pass "<<pass;
    PrintClusters(allhits,tcl, vtx);
  }

      for(unsigned short ii1 = 0; ii1 < sortindex.size() - 1; ++ii1) {
        unsigned short it1 = sortindex[ii1];
        // ignore obsolete clusters
        if(tcl[it1].ID < 0) continue;
        // ignore already attached clusters
        if(tcl[it1].BeginVtx >= 0 && tcl[it1].EndVtx >= 0) continue;
        float es1 = tcl[it1].EndSlp;
        float eth1 = atan(fScaleF * tcl[it1].EndSlp);
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float bs1 = tcl[it1].BeginSlp;
        float bth1 = atan(fScaleF * tcl[it1].BeginSlp);
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        for(unsigned short ii2 = ii1 + 1; ii2 < sortindex.size(); ++ii2) {
          unsigned short it2 = sortindex[ii2];
          if(tcl[it2].ID < 0) continue;
          // ignore already attached clusters
          if(tcl[it2].BeginVtx >= 0 && tcl[it2].EndVtx >= 0) continue;
          // try to attach cluster to existing vertices at either end
          ClusterVertex(allhits, tcl, vtx, it2);
          // ignore if both clusters are short
          if(tcl[it1].tclhits.size() < 10 &&
             tcl[it2].tclhits.size() < 10) continue;
          float es2 = tcl[it2].EndSlp;
          float eth2 = atan(fScaleF * tcl[it2].EndSlp);
          unsigned short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float bs2 = tcl[it2].BeginSlp;
          float bth2 = atan(fScaleF * tcl[it2].BeginSlp);
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
//  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
  // topo 1: check for vtx US of the ends of both clusters
          float dth = fabs(eth1 - eth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0 && dth > 0.3) {
            float dsl = es2 - es1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            float fvw = 0.5 + (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
            if(fvw > 0. && fvw < nwires) {
              // vertex wire in the detector
              unsigned short vw = fvw;
              // require vtx in the range of wires with hits AND
              // vtx US of both clusters AND
              // vtx not too far US of both clusters
              if(vw >= fFirstWire && 
                 vw <= ew1      && vw <= ew2 &&
                 vw  > ew1 - 10 && vw  > ew2 - 10) {
                float fvt = et1 + (vw - ew1) * es1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo1 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check this against existing vertices and update
                  ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 1);
                } // fvt in detector
              } // vw topo 1 check
            } // fvw in detector
          } // topo 1
  // topo 2: check for vtx US of it1 and DS of it2
          dth = fabs(eth1 - bth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0 && dth > 0.3) {
            float dsl = bs2 - es1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            float fvw = 0.5 + (et1 - ew1 * es1 - bt2 + bw2 * bs2) / dsl;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              unsigned short vw = fvw;
              // require vtx US of cluster 1 end (within 10 wires) AND
              // vtx DS of cluster 2 begin (within 10 wires)
              if(vw <= ew1      && vw >= bw2 &&
                 vw  > ew1 - 10 && vw  < ew2 + 10) {
                float fvt = et1 + (vw - ew1) * es1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo2 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                if(fvt > 0. && fvt < maxtime) {
                  ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 2);
                } // fvt in detector
              } // vw topo 2 check
            } // fvw in detector
          } // topo 2
  // topo 3: check for vtx DS of it1 and US of it2
          dth = fabs(bth1 - eth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0 && dth > 0.3) {
            float dsl = bs1 - es2;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            float fvw = 0.5 + (et2 - ew2 * es2 - bt1 + bw1 * bs1) / dsl;
            if(fvw > 0 && fvw < nwires) {
              unsigned short vw = fvw;
              // require vtx US of cluster 2 (within 10 wires) AND
              // vtx DS of cluster 1 (within 10 wires)
              if(vw <= ew2      && vw >= bw1 &&
                 vw  > ew2 - 10 && vw  < bw1 + 10) {
                float fvt = et2 + (vw - ew2) * es2;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo3 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                if(fvt > 0. && fvt < maxtime) {
                  ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 3);
                } // fvt in detector
              } // vw topo 3 check
            } // fvw in detector
          } // topo 3
  // topo 4: check for vtx DS of it1 and DS of it2
          dth = fabs(bth1 - bth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0 && dth > 0.3) {
            float dsl = bs2 - bs1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            // convert to integer if within the detector (vw, vt)
            float fvw = 0.5 + (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              unsigned short vw = fvw;
              // require vtx in the range of wires with hits AND
              // vtx DS of both clusters AND
              // vtx not too far DS of both clusters
              if(vw >= bw1 && 
                 vw >= bw2 && vw <= fLastWire &&
                 vw <  bw2 + 10 && vw <  bw1 + 10) {
                float fvt = bt1 + (vw - bw1) * bs1;
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID
      <<" topo4 vtx wire "<<vw<<" time "<<(int)fvt;
  }
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check this against existing vertices and update
                  ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 4);
                } // fvt in detector
              } // vw topo 4 check
            } // fvw in detector
          } // it2
        } // topo 4
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
      
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != clCTP) continue;
        if(vtx[iv].Wght < 0) continue;
        // fit the vertices
        float ChiDOF = 0.;
        FitVtx(tcl, vtx, iv, ChiDOF);
      }

      // set the vertex weights
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != clCTP) continue;
        // cluster weight = number of hits, saturated at 10
        float cw = tcl[it].tclhits.size();
        if(cw > 10) cw = 10;
        if(tcl[it].BeginVtx >=0) vtx[tcl[it].BeginVtx].Wght += cw;
        if(tcl[it].EndVtx >=0) vtx[tcl[it].EndVtx].Wght += cw;
      }
      

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
    } // FindVertices

/////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterVertex(std::vector<CCHitFinderAlg::CCHit>& allhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned short it)
    {
      // try to attach cluster it to an existing vertex
      
      if(vtx.size() == 0) return;
      
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        // ignore vertices in the wrong cryostat/TPC/Plane
        if(vtx[iv].CTP != clCTP) continue;
        // ignore deleted vertices
        if(vtx[iv].Wght < 0) continue;
        // determine which end to match - begin or end
        if(tcl[it].EndVtx < 0 &&  vtx[iv].Wire <= tcl[it].EndWir + 2) {
          // project cluster to US vertex
          float tdiff = tcl[it].EndTim + (vtx[iv].Wire - tcl[it].EndWir) * tcl[it].EndSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            bool SigOK = false;
            ChkSignal(allhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].EndVtx = iv;
              // re-fit it
              float ChiDOF = 99.;
              FitVtx(tcl, vtx, iv, ChiDOF);
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Add End "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<ChiDOF;
              return;
            } // SigOK
          } // tdiff
        } else if(tcl[it].BeginVtx < 0 && vtx[iv].Wire >= tcl[it].BeginWir - 2) {
          // project cluster to DS vertex
          float tdiff = tcl[it].BeginTim + (vtx[iv].Wire - tcl[it].BeginWir) * tcl[it].BeginSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            bool SigOK = false;
            ChkSignal(allhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].BeginVtx = iv;
              // re-fit it
              float ChiDOF = 99.;
              FitVtx(tcl, vtx, iv, ChiDOF);
  if(vtxprt) mf::LogVerbatim("ClusterCrawler")<<"Add Begin "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<ChiDOF;
              return;
            } // SigOK
          } // tdiff
        } // vtx[iv].Wire
      } // iv
    }



/////////////////////////////////////////
    void ClusterCrawlerAlg::ChkVertex(std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        short vw, float fvt, unsigned short it1, unsigned short it2, short topo)
      {
        // Checks the vertex (vw, vt) against the existing set of vertices.
        // If there a match, clusters it1 and/or it2 are associated with it
        // if there is signal between the existing vertex and the start of
        // the cluster. The topo flag indicates the type of vertex that was
        // found: 1 = US of it1 && US of it2, 2 = US of it1 && DS of it2,
        // 3 = DS of it1 && US of it2, 4 = DS of it1 and DS of it2.
        // didit is set true if a cluster is attached to a (new or old) vertex
        
  if(vtxprt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"ChkVertex "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" topo "<<topo
      <<" vw "<<vw<<" vt "<<(int)fvt;
  }
        // check vertex and clusters for proximity to existing vertices
        bool SigOK = false;
        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if(vtx[iv].CTP != clCTP) continue;
          if( abs( vw - vtx[iv].Wire) < 4 &&
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
            } // cluster it1
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
          newvx.Wght = 0;
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
      // (wire1, time1) and (wire2, time2).
      SigOK = false;
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
      if(wireb == wiree) {
        clpar[1] = 0.;
      } else {
        clpar[1] = (timeb - timee) / (wireb - wiree);
      }
      for(short wire = wiree; wire < wireb + 1; ++wire) {
        // assume there is no signal on this wire
        bool WireSigOK = false;
        float prtime = timee + (wire - wire0) * clpar[1];
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
          // skip obsolete hits
          if(allhits[khit].Charge < 0) continue;
          // outside the hit RAT?
          if(prtime > allhits[khit].HiTime) continue;
          if(prtime < allhits[khit].LoTime) continue;
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
      clBeginWir = tcl[icl].BeginWir;
      clBeginTim = tcl[icl].BeginTim;
      clBeginChg = tcl[icl].BeginChg;
      clStopCode = 5;
      clProcCode = tcl[icl].ProcCode;
      fcl2hits.clear();
      chifits.clear();
      for(unsigned short ii = 0; ii < pos; ++ii) {
        unsigned short iht = tcl[icl].tclhits[ii];
  if(allhits[iht].InClus != tcl[icl].ID) {
    std::cout<<"SplitCluster bad hit "<<iht<<" "<<allhits[iht].InClus
      <<" "<<tcl[icl].ID<<" ProcCode "<<clProcCode<<std::endl;
    return;
  }
        fcl2hits.push_back(iht);
      }
      // determine the pass in which this cluster was created
      pass = tcl[icl].ProcCode - 10 * (tcl[icl].ProcCode / 10);
      if(pass > fNumPass-1) pass = fNumPass-1;
      // fit the end hits
      FitCluster(allhits);
      clEndSlp = clpar[1];
      clEndSlpErr = clparerr[1];
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
      clEndWir = tcl[icl].EndWir;
      clEndTim = tcl[icl].EndTim;
      clEndChg = tcl[icl].EndChg;
      clStopCode = 5;
      clProcCode = tcl[icl].ProcCode;
      fcl2hits.clear();
      chifits.clear();
      for(unsigned short ii = pos; ii < tcl[icl].tclhits.size(); ++ii) {
        unsigned short iht = tcl[icl].tclhits[ii];
  if(allhits[iht].InClus != tcl[icl].ID) {
    std::cout<<"SplitCluster bad hit "<<iht<<" "<<allhits[iht].InClus
      <<" "<<tcl[icl].ID<<" ProcCode "<<clProcCode<<std::endl;
  }
        fcl2hits.push_back(iht);
        // define the Begin parameters
        if(fcl2hits.size() == fMaxHitsFit[pass] ||
           fcl2hits.size() == fMinHits[pass]) {
          FitCluster(allhits);
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
        }
        if((unsigned short)fcl2hits.size() == fNHitsAve[pass] + 1) {
          FitClusterChg(allhits);
          clBeginChg = fAveChg;
        }
      } // ii
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
      // is done.

      prt = (fDebugPlane == (short)plane && fDebugWire < 0);
      
      unsigned short tclsize = tcl.size();
      
      for(unsigned short it1 = 0; it1 < tclsize - 1; ++it1) {
        // ignore already merged clusters
        if(tcl[it1].ID < 0) continue;
        if(tcl[it1].CTP != clCTP) continue;
        float bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        float arg = fScaleF * bs1;
        float bth1 = atan(arg);
        // more compact notation for begin/end, wire/time/chg/slp/theta, 1/2
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        float bc1 = tcl[it1].BeginChg;
        float es1 = tcl[it1].EndSlp;
        // convert slope to angle
        arg = fScaleF * es1;
        float eth1 = atan(arg);
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float ec1 = tcl[it1].EndChg;
        unsigned short nh1 = tcl[it1].tclhits.size();
        unsigned short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
        if(pass1 > fNumPass) pass1 = fNumPass;
        for(unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
          // ignore already merged clusters
          if(tcl[it1].ID < 0) continue;
          if(tcl[it2].ID < 0) continue;
          // only merge if they are in the right cryostat/TPC/plane
          if(tcl[it2].CTP != clCTP) continue;
          float bs2 = tcl[it2].BeginSlp;
          float bs2e = tcl[it2].BeginSlpErr;
          // convert slope to angle
          arg = fScaleF * bs2;
          float bth2 = atan(arg);
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
          float bc2 = tcl[it2].BeginChg;
          float es2 = tcl[it2].EndSlp;
          arg = fScaleF * es2;
          float eth2 = atan(arg);
          unsigned short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float ec2 = tcl[it2].EndChg;
          unsigned short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
          if(pass2 > fNumPass) pass2 = fNumPass;
          // use the more promiscuous pass for cuts
          float angcut = fKinkAngCut[pass1];
          if(fKinkAngCut[pass2] > angcut) angcut = fKinkAngCut[pass2];
          unsigned short skipcut = fMaxWirSkip[pass1];
          if(fMaxWirSkip[pass2] > skipcut) skipcut = fMaxWirSkip[pass2];
          float chgcut = fChgCut[pass1];
          if(fChgCut[pass2] > chgcut) chgcut = fChgCut[pass2];
          float timecut = fTimeDelta[pass];
          if(fTimeDelta[pass2] > timecut) timecut = fTimeDelta[pass2];
          // increase the time cut for large angle clusters
          timecut *= (2 - 1/(1 + fabs(clpar[1])));
          
          bool bothLong = (nh1 > 5 && tcl[it2].tclhits.size() > 5);
          float dth = fabs(bth2 - eth1);

  if(prt && bw2 < ew1 ) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(ew1 - bw2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }

          if(bw2 < ew1 && (ew1 - bw2)  < skipcut && dth < angcut) {
            // look for US and DS broken clusters w similar angle.
            // US cluster 2 merge with DS cluster 1?
            // This is the most likely occurrence given the order in which
            // clusters are created so put it first.
            float chgrat = 2 * fabs((bc2 - ec1) / (bc2 + ec1));
            // ignore the charge cut for large angle clusters
            if(fabs(es1) > fLAClusSlopeCut) chgrat = 0.;
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.1) chgrat = 0.;
            // project bw2,bt2 to ew1
            float dtim = fabs(bt2 + (ew1-bw2)*bs2 - et1);
            // error on the projection
            float dtime = fabs(bt2 + (ew1-bw2)*bs2e - et1);
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<" dtim "<<dtim<<" err "<<dtime<<" timecut "<<(int)timecut
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              // ensure there is a signal between cluster ends
              bool SigOK = true;
              if(fAllowNoHitWire == 0) ChkSignal(allhits,ew1,et1,bw2,bt2,SigOK);
              if(SigOK) {
                DoMerge(allhits, tcl, vtx, it2, it1, 10);
                tclsize = tcl.size();
                break;
              }
            }
          } // US cluster 2 merge with DS cluster 1?
          
          dth = fabs(bth1 - eth2);
  if(prt && bw1 < ew2 && (ew2 - bw1)  < skipcut) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Chk2 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(bw1 - ew2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          if( bw1 < ew2 && (ew2 - bw1)  < skipcut && dth < angcut ) {
            // look for US and DS broken clusters w similar angle
            // US cluster 1 merge with DS cluster 2?
            float chgrat = 2 * fabs((bc1 - ec2) / (bc1 + ec2));
            // ignore the charge cut for large angle clusters
            if(fabs(es2) > fLAClusSlopeCut) chgrat = 0.;
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.1) chgrat = 0.;
            // project sw1,st1 to ew2
            float dtim = fabs(bt1 + (ew2-bw1)*bs1 - et2);
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<" dtim "<<dtim<<" err "<<dtim<<" timecut "<<(int)timecut
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              bool SigOK = true;
              if(fAllowNoHitWire == 0) ChkSignal(allhits,bw1,bt1,ew2,et2,SigOK);
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
            float dth = fabs(eth2 - eth1);
            float dtim = fabs(et1 +(ew2 - ew1 - 1)*es1 - et2);
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
            float dth = fabs(eth2 - eth1);
            float dtim = fabs(et2 +(ew1 - ew2 - 1)*es2 - et1);
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"US chgrat "<<chgrat<<" cut "<<fChgCut[pass];
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * fabs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"DS chgrat "<<chgrat<<" cut "<<fChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    SigOK = false;
    ChkSignal(allhits, wiron1, timon1, bw2, bt2, SigOK);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"DS SigOK? "<<SigOK;
    if( !SigOK ) return;

  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Merge em";
    // success. Merge them
    DoMerge(allhits, tcl, vtx, it1, it2, 100);
    didit = true;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::DoMerge(
      std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      unsigned short it1, unsigned short it2, short inProcCode)
  {
    // Merge clusters. Cluster 1 has precedence for assignment of hits
    
    ClusterStore& cl1 = tcl[it1];
    ClusterStore& cl2 = tcl[it2];
    // mark cl1 and cl2 obsolete
    cl1.ID = -cl1.ID;
    cl2.ID = -cl2.ID;
    
    // find the low and high wire for both clusters
    unsigned short hiwire = cl1.BeginWir;
    if(cl2.BeginWir > hiwire) hiwire = cl2.BeginWir;
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
    } // iht
    // now cluster 1
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned short hit = cl1.tclhits[iht];
      // un-assign the hit from this cluster
      allhits[hit].InClus = 0;
      unsigned short wire = allhits[hit].WireNum;
      unsigned short index = wire - lowire;
      wirehit[index] = hit;
    } // iht
    // make the new cluster
    fcl2hits.clear();
    chifits.clear();
    for(unsigned short ii = wirehit.size()-1; ii > 0; --ii) {
      if(wirehit[ii] >= 0) {
        unsigned short hit = wirehit[ii];
        fcl2hits.push_back(hit);
        if(fcl2hits.size() == 4) {
          // re-fit the Begin end of the cluster
          FitCluster(allhits);
          if(clChisq > 99.) {
            mf::LogError("ClusterCrawler")<<"cl2DoMerge bad Begin fit "<<clChisq;
            return;
          }
          // define the first wire/time
          unsigned short jj = fcl2hits[0];
          clBeginWir = allhits[jj].WireNum;
          clBeginTim = allhits[jj].Time;
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
          FitClusterChg(allhits);
          clBeginChg = fAveChg;
          // reset fAveChg in case the charge fit fails at the end
          fAveChg = -1;
        } // fcl2hits.size() == 4
      } // hit >= 0
    } // ii
    // re-fit the End of the cluster with the current pass params
    FitCluster(allhits);
    if(clChisq > 99.) {
      mf::LogError("ClusterCrawler")<<"DoMerge bad End fit "<<clChisq;
      return;
    } // clChisq > 99
    // define the last wire/time
    unsigned short jj = fcl2hits[fcl2hits.size()-1];
    clEndWir = allhits[jj].WireNum;
    clEndTim = allhits[jj].Time;
    clEndSlp = clpar[1];
    clEndSlpErr = clparerr[1];
    // define the end charge
    FitClusterChg(allhits);
    clEndChg = fAveChg;

    clStopCode = cl1.StopCode;
    // append it to the tcl vector
    TmpStore(allhits, tcl);
    unsigned short itnew = tcl.size()-1;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"DoMerge "<<cl1.ID<<" "<<cl2.ID<<" -> "<<tcl[itnew].ID;
    // stuff the processor code with the current pass
    tcl[itnew].ProcCode = inProcCode + pass;
    // delete a vertex between these two?
    if(tcl[it1].BeginVtx >= 0) {
      if(tcl[it1].BeginVtx == tcl[it2].EndVtx) {
        vtx[tcl[it2].EndVtx].Wght = -1;
        tcl[it1].BeginVtx = -99;
        tcl[it2].EndVtx = -99;
      }
    } else if(tcl[it1].EndVtx >= 0) {
      if(tcl[it1].EndVtx == tcl[it2].BeginVtx) {
        vtx[tcl[it2].BeginVtx].Wght = -1;
        tcl[it1].EndVtx = -99;
        tcl[it2].BeginVtx = -99;
      }
    } // tcl[it1].BeginVtx
    // try to preserve the vertex assignments
    tcl[itnew].EndVtx = tcl[it1].EndVtx;
    if(tcl[itnew].EndVtx < 0 && tcl[it2].EndVtx > 0) 
       tcl[itnew].EndVtx = tcl[it2].EndVtx;

    tcl[itnew].BeginVtx = tcl[it1].BeginVtx;
    if(tcl[itnew].BeginVtx < 0 && tcl[it2].BeginVtx > 0) 
       tcl[itnew].BeginVtx = tcl[it2].BeginVtx;

  // check for errors
  for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
    if(tcl[icl].ID < 0) continue;
    for(unsigned short ii = 0; ii< tcl[icl].tclhits.size(); ++ii) {
      unsigned short iht = tcl[icl].tclhits[ii];
      if(allhits[iht].InClus != tcl[icl].ID) {
        std::cout<<"DoMerge Bad "<<tcl[icl].ID<<" "<<allhits[iht].InClus
          <<" "<<tcl[icl].ProcCode<<" "<<iht<<std::endl;
        std::cout<<"Merging "<<tcl[it1].ID<<" "<<tcl[it2].ID
          <<" --> "<<tcl[itnew].ID<<std::endl;
        return;
      }
    } // ii
  } // icl


  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::PrintClusters(
    std::vector<CCHitFinderAlg::CCHit>& /*allhits*/, 
     std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
  {
    // prints clusters to the screen for code development
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"  ID CTP nht Stop  Proc  beg_W:T  bTheta Therr"
      <<" begChg  end_W:T  eTheta Therr endChg bVx eVx\n";
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
      float arg = fScaleF * tcl[ii].BeginSlp;
      float theta = atan(arg);
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<theta;
      float thetaerr =  fScaleF * tcl[ii].BeginSlpErr / (1 + arg * arg);
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<thetaerr;
      myprt<<std::right<<std::setw(5)<<(short)tcl[ii].BeginChg;
      iTime = tcl[ii].EndTim;
      myprt<<std::right<<std::setw(6)<<tcl[ii].EndWir<<":"<<iTime;
      if(iTime < 10) {
        myprt<<"   ";
      } else if(iTime < 100) {
        myprt<<"  ";
      } else if(iTime < 1000) myprt<<" ";
      arg = fScaleF * tcl[ii].EndSlp;
      theta = atan(arg);
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<theta;
      thetaerr =  fScaleF * tcl[ii].EndSlpErr / (1 + arg * arg);
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<thetaerr;
      myprt<<std::right<<std::setw(5)<<(short)tcl[ii].EndChg;
      myprt<<std::right<<std::setw(5)<<tcl[ii].BeginVtx;
      myprt<<std::right<<std::setw(5)<<tcl[ii].EndVtx;
      myprt<<"\n";
    } // ii
    // print out vertices
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      mf::LogVerbatim("ClusterCrawler vertices")
        <<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time
        <<" wght "<<(int)vtx[iv].Wght<<" topo "<<vtx[iv].Topo;
    }    
  } // cl2Print

/////////////////////////////////////////
    void ClusterCrawlerAlg::TmpGet(std::vector<CCHitFinderAlg::CCHit>& /*allhits*/,
        std::vector<ClusterStore>& tcl, unsigned short it1)
    {
      // copies temp cluster it1 into the fcl2hits vector, etc. This is 
      // effectively the inverse of cl2TmpStore
      
      if(it1 > tcl.size()) return;


      clBeginSlp = tcl[it1].BeginSlp;
      clBeginSlpErr = tcl[it1].BeginSlpErr;
      clBeginWir = tcl[it1].BeginWir;
      clBeginTim = tcl[it1].BeginTim;
      clBeginChg = tcl[it1].BeginChg;
      clEndSlp = tcl[it1].EndSlp;
      clEndSlpErr = tcl[it1].EndSlpErr;
      clEndWir = tcl[it1].EndWir;
      clEndTim = tcl[it1].EndTim;
      clEndChg = tcl[it1].EndChg;
      clStopCode = tcl[it1].StopCode;
      clProcCode = tcl[it1].ProcCode;
      clAssn = tcl[it1].Assn;
      clCTP = tcl[it1].CTP;
      fcl2hits = tcl[it1].tclhits;
    }


/////////////////////////////////////////
  void ClusterCrawlerAlg::TmpStore(std::vector<CCHitFinderAlg::CCHit>& allhits, 
     std::vector<ClusterStore>& tcl)
  {

    if(fcl2hits.size() < 3) return;
    
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
    if(clEndChg < 0.) {
      // use the next to the last two hits. The End hit may have low charge
      unsigned int ih0 = fcl2hits.size() - 2;
      unsigned int hit = fcl2hits[ih0];
      clEndChg = allhits[hit].Charge;
      hit = fcl2hits[ih0 - 1];
      clEndChg += allhits[hit].Charge;
      clEndChg = clEndChg / 2.;
    }
    if(clBeginChg < 0.) {
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
    clstr.BeginWir    = allhits[hitb].WireNum;
    clstr.BeginTim    = allhits[hitb].Time;
    clstr.BeginChg    = clBeginChg;
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndWir      = allhits[hite].WireNum;
    clstr.EndTim      = allhits[hite].Time;
    clstr.EndChg      = clEndChg;
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.Assn        = clAssn;
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
  if(dwir == fDebugWire && fDebugHit > 0)
    prt = abs(allhits[dhit].Time - fDebugHit) < 20;

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"******************* LACrawlUS PASS "<<pass<<" Hits ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned short iht = fcl2hits[fcl2hits.size() - 1 - ii];
      myprt<<allhits[iht].WireNum<<":"<<(int)allhits[iht].Time<<" ";
    }
  }

    bool SigOK = true;
    bool HitOK = true;
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
      if(!HitOK) continue;
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
      if(chsiz < 8) continue;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
      if( chifits[chsiz-2] > fKinkChiRat[pass] * chifits[chsiz-3] && 
          chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
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
          // hit a kink. Lop of the first 3 hits, refit and stop crawling
          for(short jj = 0; jj < 3; ++jj) {
            fcl2hits.pop_back();
            chifits.pop_back();
          }
          FitCluster(allhits);
          // set the kink stop code
          clStopCode = 3;
          break;
        } // kinkang check
      } // chifits test
      // chisq check
      if(clChisq > 10.) {
  if(prt)mf::LogVerbatim("ClusterCrawler")<<" Bad chisq "<<clChisq;
        for(unsigned short nlop = 0; nlop < 4; ++nlop) {
          unsigned short cfsize = chifits.size() - 1;
          // stop lopping if chisq is good
          if(chifits[cfsize] < 1.) break;
          float chirat = chifits[cfsize] / chifits[cfsize - 1];
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"chirat "<<chirat
    <<" last hit "<<fcl2hits[fcl2hits.size()-1];
          if(chirat < 1.2) break;
          fcl2hits.pop_back();
          chifits.pop_back();
          if(fcl2hits.size() < 4) break;
          if(chifits.size() < 4) break;
        } // nlop
        FitCluster(allhits);
        clStopCode = 4;
      } // lChisq > fChiCut[pass]
    } // nextwire
/*
    // fit the charge at the end
    FitClusterChg(allhits);
    clEndChg = fAveChg;
*/
    clProcCode += 300;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"LACrawlUS done ";
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
  if(dwir == fDebugWire && fDebugHit > 0)
    prt = abs(allhits[dhit].Time - fDebugHit) < 10;

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
          if(PostSkip && nmissed > 1) {
            if((short)(fcl2hits.size() - nHitAfterSkip) < 4) {
              fcl2hits.clear();
              return;
            }
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" PostSkip && nmissed = "<<nmissed;
            for(short jj = 0; jj < nHitAfterSkip; ++jj) {
              fcl2hits.pop_back();
              chifits.pop_back();
            } // pop_back
            FitCluster(allhits);
            if(clChisq > 90.) {
              fcl2hits.clear();
              return;
            }
            FitCluster(allhits);
            clStopCode = 2;
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
        if(chifits.size() > 7 && fKinkChiRat[pass] > 0) {
  if(chifits.size() != fcl2hits.size()) {
    mf::LogVerbatim("ClusterCrawler")<<"chifits problem "<<chifits.size()<<" "<<fcl2hits.size();
    unsigned short hh = fcl2hits[0];
    mf::LogVerbatim("ClusterCrawler")<<" "<<allhits[hh].WireNum<<":"<<hh;
  }
          unsigned short chsiz = chifits.size()-1;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
          if( chifits[chsiz-2] > fKinkChiRat[pass] * chifits[chsiz-3] && 
              chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
              chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
            if(fcl2hits.size() < 8) {
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
              }
              FitCluster(allhits);
              FitClusterChg(allhits);
              // set the kink stop code but keep looking
              clStopCode = 3;
//              break;
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
        // set the Begin charge after fAveChg is defined
        if(clBeginChg < 0 && fcl2hits.size() >= fNHitsAve[pass]) {
          FitClusterChg(allhits);
          // project the charge to the Begin end of the cluster
          clBeginChg = fAveChg + (clBeginWir - nextwire) * fChgSlp;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Set clBeginChg "<<clBeginChg;
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
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Check nAdjHit "<<nAdjHit;
      // lop off hits?
      if(nAdjHit < fMinWirAfterSkip[pass]) {
        for(unsigned short ii = 0; ii < nAdjHit + 1; ++ii) {
          fcl2hits.pop_back();
          chifits.pop_back();
          reFit = true;
        }
      }
    } // fcl2hits.size() > fMinWirAfterSkip[pass] + 3
    if(reFit) {
      FitCluster(allhits);
      FitClusterChg(allhits);
    }
    clEndChg = fAveChg;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"CrawlUS done ";
    prt = false;
  }

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
    // fAveChg, fAveAmp and fAveRMS is re-calculated as well.
    
    
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
      fAveRMS = 0.;
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
          fAveAmp += allhits[ihit].Amplitude;
          fAveRMS += allhits[ihit].RMS;
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
      fAveRMS = 0.;
      unsigned short hitcnt = 0;
      bool UseEm = false;
//      for(unsigned short it = cls.tclhits.size() - 1; it >= 0; --it) {
//        unsigned short ihit = cls.tclhits[it];
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
          fAveAmp += allhits[ihit].Amplitude;
          fAveRMS += allhits[ihit].RMS;
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
    LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
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
    // apply an angle dependent scale factor. The error should be
    // wire pitch / sqrt(12) for a cluster at 90 degrees. This formula 
    // simply doubles the error, which I think is reasonable for uBooNE and
    // ArgoNeuT 
    float angfactor = 2 - 1/(1 + fabs(clpar[1]));

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
      float terr = fHitErrFac * allhits[ihit].RMS;
      ytimerr2.push_back(angfactor * terr * terr);
      if(iht == nht) break;
      ++iht;
    }

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"ClusterFit W:T ";
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
    LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
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
    } else if((unsigned short)fcl2hits.size() >= fNHitsAve[pass]){
      // do a real fit
      std::vector<float> xwir;
      std::vector<float> ychg;
      std::vector<float> ychgerr2;
      unsigned short wire0 = 0;
      bool first = true;
      unsigned int npt = 0;
      // this loop intentionally ignores the Begin hit
      for(unsigned int ii = fcl2hits.size() - 1; ii > 0; --ii) {
        unsigned short wire = allhits[fcl2hits[ii]].WireNum;
        if(first) {
          wire0 = wire;
          first = false;
        }
        xwir.push_back((float)(wire - wire0));
        float chg = allhits[fcl2hits[ii]].Charge;
        ychg.push_back(chg);
        ychgerr2.push_back(0.5 * chg);
        if(npt == fNHitsAve[pass]) break;
        ++npt;
      }
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"FitChg Wire: "<<wire0<<" nht "<<ychg.size();
      if(ychg.size() < 3) return;
      float intcpt; float slope; float intcpterr;
      float slopeerr; float chidof;
      LinFit(xwir, ychg, ychgerr2, intcpt, slope, intcpterr, slopeerr, chidof);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"FitChg chidof "<<chidof;
      if(chidof > 20.) return;
      fAveChg = intcpt;
      fChgSlp = slope;
    }
  } // fitchg

/////////////////////////////////////////
  void ClusterCrawlerAlg::AddLAHit(std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<VtxStore>& vtx,
    unsigned short kwire, bool& ChkCharge, bool& HitOK, bool& SigOK)
  {
    // A variant of AddHit for large angle clusters. The main differences are
    // that 1) no charge similarity cut is applied and 2) the projected cluster
    // position must lie within the time window of a hit multiplet (if one
    // is found)
    
    SigOK = false;
    HitOK = false;
    
    // not in the range of wires with hits
    if(kwire < fFirstWire || kwire > fLastWire) return;

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
    
    // the last hit added to the cluster
    unsigned short lastClHit = fcl2hits[fcl2hits.size()-1];
    unsigned short wire0 = allhits[lastClHit].WireNum;
    // the projected time of the cluster on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"AddLAHit: prtime= "<<(short)prtime;

    unsigned short imbest = 0;
    float best = 999.;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(allhits[khit].InClus < 0) continue;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Chk W:T "<<kwire<<":"<<(short)allhits[khit].Time
    <<" InClus "<<allhits[khit].InClus
    <<" mult "<<allhits[khit].numHits
    <<" RMS "<<allhits[khit].RMS
    <<" LoTime "<<(int)allhits[khit].LoTime
    <<" HiTime "<<(int)allhits[khit].HiTime;
      // projected time within the Signal time window?
      if(prtime < allhits[khit].LoTime) continue;
      if(prtime > allhits[khit].HiTime) continue;
      SigOK = true;
      // hit used?
      if(allhits[khit].InClus > 0) continue;
      // projected time within the Hit time window?
      HitOK = true;
      float dtime = fabs(prtime - allhits[khit].Time);
      if(dtime < best) {
        imbest = khit;
        best = dtime;
      }
    } // khit
    
  if(prt) {
    if(!HitOK) mf::LogVerbatim("ClusterCrawler")<<" no hit found ";
  }
    if(!HitOK) return;

  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Pick hit time "<<(int)allhits[imbest].Time;
    
    // merge the hits in a multiplet?
    if(allhits[imbest].numHits > 1) {
      bool doMerge = true;
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
        // look for a big separation between adjacent hits
        for(unsigned short jj = 0; jj < allhits[imbest].numHits; ++jj) {
          unsigned short jht = allhits[imbest].LoHitID + jj;
          if(allhits[jht].InClus > 0) ++nused;
          if(allhits[jht].InClus == 0) multipletChg += allhits[jht].Charge;
          // check the neighbor hit separation
          if(jj > 0) {
            // pick the larger RMS of the two hits
            float hitRMS = allhits[jht].RMS;
            if(allhits[jht - 1].RMS > hitRMS) hitRMS = allhits[jht-1].RMS;
            float tdiff = fabs(allhits[jht].Time - allhits[jht-1].Time) / hitRMS;
            if(tdiff > 2.5) doMerge = false;
          } // jj > 0
        } // jj
  if(prt) {
    if(!doMerge) mf::LogVerbatim("ClusterCrawler")
      <<" Hits are well separated. Don't merge them";
  }
        if(doMerge && nused == 0) {
          // compare the charge with the last hit added?
          if(ChkCharge) {
            float chgrat = multipletChg / allhits[lastClHit].Charge;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" merge hits charge check "
    <<(int)multipletChg<<" Previous hit charge "<<(int)allhits[lastClHit].Charge;
            if(chgrat > 1.7) doMerge = false;
          }
        } // doMerge && nused == 0
      } // doMerge true
      if(doMerge) MergeHits(allhits, imbest);
    } // allhits[imbest].numHits > 1
    
    // attach to the cluster and fit
    fcl2hits.push_back(imbest);
    FitCluster(allhits);
    chifits.push_back(clChisq);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" >>ADD W:T "<<kwire<<":"<<(short)allhits[imbest].Time
    <<" clChisq "<<clChisq;
    

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
        (wire0 - kwire) > fAllowNoHitWire) return;
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
    float angfactor = 2 - 1/(1 + fabs(clpar[1]));
    float hiterr = angfactor * fHitErrFac * allhits[lastClHit].RMS;
    float err = sqrt(prtimerr2 + hiterr * hiterr);
    // Time window for accepting a hit.
    float prtimeLo = prtime - 3 * err;
    float prtimeHi = prtime + 3 * err;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"AddHit: prtime Lo "<<(int)prtimeLo<<" Hi "<<(int)prtimeHi
    <<" fAveChg "<<(int)fAveChg;

    // loop through the hits
    short imbest = -1;
    float best = 9999.;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(allhits[khit].InClus < 0) continue;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<" Chk W:T "<<kwire<<":"<<(short)allhits[khit].Time
    <<" InClus "<<allhits[khit].InClus
    <<" mult "<<allhits[khit].numHits
    <<" RMS "<<allhits[khit].RMS
    <<" Charge "<<(int)allhits[khit].Charge
    <<" LoTime "<<(int)allhits[khit].LoTime
    <<" HiTime "<<(int)allhits[khit].HiTime;
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

    // Make a charge similarity cut if the average charge is defined
    bool fitChg = true;
    if(fAveChg > 0) {

      float chgrat = (allhits[imbest].Charge - fAveChg) / fAveChg;

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
        lasthitlow = ( lastchgrat < fChgCut[pass]);
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
        if(prt) mf::LogVerbatim("ClusterCrawler")<<" fails high charge && bad delta T";
        if(best > 1.5 * err) return;
      }

      // decide whether to fit the charge
      fitChg = ( !lasthitlow && !lasthitbig && chgrat < fabs(fChgCut[pass]) );
    } // fAveChg > 0
    
    // we now have a hit that meets all the criteria. Fit it
    fcl2hits.push_back(imbest);
    FitCluster(allhits);
    chifits.push_back(clChisq);
    HitOK = true;
  if(prt) mf::LogVerbatim("ClusterCrawler")
      <<" >>ADD W:T "<<kwire<<":"<<(short)allhits[imbest].Time<<" best "<<best
      <<" clChisq "<<clChisq;
    if(!fitChg) return;
  if(prt) mf::LogVerbatim("ClusterCrawler")<<" Fit charge ";
    FitClusterChg(allhits);
  } // AddHit

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
      LinFit(x, y, ey2, tv, wv, tverr, wverr, ChiDOF);
      if(ChiDOF < 5) {
        vtx[iv].Wire = (int)(wv + 0.5);
        vtx[iv].Time = -tv;
        if(vtx[iv].Time < 0 || vtx[iv].Time > 3200) {
          //mf::LogError("ClusterCrawler")<<"FitVtx: Bad fit time "<<vtx[iv].Time
	  //<<" on vtx "<<iv; // commenting out as it gives incorrect message in 1-big-window regime
        }
      } // ChiDOF < 5
    } // FitVtx


/////////////////////////////////////////
    void ClusterCrawlerAlg::LinFit(std::vector<float>& x, std::vector<float>& y, 
      std::vector<float>& ey2, float& Intercept, float& Slope, 
      float& InterceptError, float& SlopeError, float& ChiDOF) 
    {
      // fit a line ala Bevington linfit.F. The number of points fit is defined by
      // the size of the y vector. 

      ChiDOF = 999.;

      if(y.size() < 2) return;
      if(x.size() < y.size() || ey2.size() < y.size()) return;
      
      float sum = 0.;
      float sumx = 0.;
      float sumy = 0.;
      float sumxy = 0.;
      float sumx2 = 0.;
      float sumy2 = 0.;

      for(unsigned short ii = 0; ii < y.size(); ++ii) {
        float weight = 1. / ey2[ii];
        sum += weight;
        sumx += weight * x[ii];
        sumy += weight * y[ii];
        sumx2 += weight * x[ii] * x[ii];
        sumxy += weight * x[ii] * y[ii];
        sumy2 += weight * y[ii] * y[ii];
      }
      // calculate coefficients and std dev
      float delta = sum * sumx2 - sumx * sumx;
      if(delta == 0.) return;
      float A = (sumx2 * sumy - sumx * sumxy) / delta;
      float B = (sumxy * sum  - sumx * sumy) / delta;
      Intercept = A;
      Slope = B;
      if(x.size() == 2) {
        ChiDOF = 0.;
        return;
      }
      float ndof = x.size() - 2;
      float varnce = (sumy2 + A*A*sum + B*B*sumx2 - 
                      2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
      if(varnce > 0.) {
        InterceptError = sqrt(varnce * sumx2 / delta);
        SlopeError = sqrt(varnce * sum / delta);
      } else {
        InterceptError = 0.;
        SlopeError = 0.;
      }
      sum = 0.;
      // calculate chisq
      for(unsigned short ii = 0; ii < y.size(); ++ii) {
        float arg = y[ii] - A - B * x[ii];
        sum += arg * arg / ey2[ii];
      }
      ChiDOF = sum / ndof;
    }

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
        art::Ptr<recob::Wire> theWire = allhits[hit].Wire;
        uint32_t channel = theWire->Channel();
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        if(wids[0].Plane != planeID.Plane) continue;
        if(wids[0].TPC != planeID.TPC) continue;
        if(wids[0].Cryostat != planeID.Cryostat) continue;
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
    void ClusterCrawlerAlg::cl2SortByLength(std::vector<ClusterStore>& tcl,
        std::map<unsigned short, unsigned short>& sortindex)
    {
      // sorts the temporary cluster vector by decreasing number of hits,
      // while ignoring abandoned clusters. Returns index map with the
      // sort order
      
      // form a vector of pairs of the number of hits and the index
      std::vector< std::pair<unsigned short, unsigned short> > index;
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        if(tcl[ii].ID > 0 && tcl[ii].CTP == clCTP) 
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
