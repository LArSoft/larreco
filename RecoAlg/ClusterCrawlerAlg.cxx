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

struct CluLen{
  int index;
  int length;
};

bool greaterThan (CluLen c1, CluLen c2) { return (c1.length > c2.length);}
bool lessThan (CluLen c1, CluLen c2) { return (c1.length < c2.length);}

namespace cluster {

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
		
		fMinAmp 						= pset.get< float >("MinAmp");
		fChgNearWindow			= pset.get< float >("ChgNearWindow");
		fChgNearCut 				= pset.get< float >("ChgNearCut");

    fChkClusterDS       = pset.get< bool   >("ChkClusterDS");
    fVtxClusterSplit    = pset.get< bool   >("VtxClusterSplit");
    fFindStarVertices   = pset.get< bool   >("FindStarVertices");
    fHitErrFac          = pset.get< float  >("HitErrFac");
    fMinHitFrac         = pset.get< float  >("MinHitFrac");
    
    fLAClusAngleCut     = pset.get< float  >("LAClusAngleCut");
		fLAClusMaxHitsFit 	= pset.get< unsigned short>("LAClusMaxHitsFit");
    fHitMergeChiCut     = pset.get< float  >("HitMergeChiCut");
    fMergeOverlapAngCut = pset.get< float  >("MergeOverlapAngCut");
    fAllowNoHitWire     = pset.get< unsigned short  >("AllowNoHitWire");
		fVertex2DCut				= pset.get< float  >("Vertex2DCut");
    fVertex3DCut        = pset.get< float  >("Vertex3DCut");
    
    fDebugPlane         = pset.get< short  >("DebugPlane");
    fDebugWire          = pset.get< short  >("DebugWire");
    fDebugHit           = pset.get< short  >("DebugHit");

    fuBCode             = pset.get< bool >("uBCode", false);

    
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
  bool SortByLowHit(unsigned short i, unsigned short j) {return ((i > j));}

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
    inClus.clear();
  } // ClusterCrawlerAlg::ClearResults()
  
  
//------------------------------------------------------------------------------
  void ClusterCrawlerAlg::CrawlInit() {
    prt = false; vtxprt = false;
    NClusters = 0;  clBeginSlp = 0; clBeginSlpErr = 0; clBeginTim = 0;
    clBeginWir = 0; clBeginChg = 0; clBeginChgNear = 0; clEndSlp = 0;      clEndSlpErr = 0;
    clEndTim = 0;   clEndWir = 0;   clEndChg = 0; clEndChgNear = 0; clChisq = 0;
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
    
    inClus.resize(fHits.size());
    for(unsigned int iht = 0; iht < inClus.size(); ++iht) inClus[iht] = 0;
    
    // don't do anything...
    if(fNumPass == 0) return;
    
    if(fHits.size() < 3) return;
    if(fHits.size() > USHRT_MAX) { // not really, but let's keep things sane
      mf::LogWarning("CC")<<"Too many hits for ClusterCrawler "<<fHits.size();
      return;
    }

    const dataprov::LArProperties* larprop = art::ServiceHandle<util::LArPropertiesService>()->getLArProperties();
    const dataprov::DetectorProperties* detprop = art::ServiceHandle<util::DetectorPropertiesService>()->getDetectorProperties();

    for (geo::TPCID const& tpcid: geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      for(plane = 0; plane < TPC.Nplanes(); ++plane){
        WireHitRange.clear();
        // define a code to ensure clusters are compared within the same plane
        clCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
        cstat = tpcid.Cryostat;
        tpc = tpcid.TPC;
        // fill the WireHitRange vector with first/last hit on each wire
        // dead wires and wires with no hits are flagged < 0
        GetHitRange(clCTP, WireHitRange, fFirstWire, fLastWire);

// sanity check
/*
  std::cout<<"Plane "<<plane<<" wire range "<<fFirstWire<<" "<<fLastWire;
  unsigned int nhts = 0;
  for(unsigned short wire = fFirstWire; wire < fLastWire; ++wire) {
    unsigned short index = wire - fFirstWire;
    if(index > fHits.size() - 1) continue;
    if(WireHitRange[index].first < 0) continue;
    unsigned int fhit = WireHitRange[index].first;
    unsigned int lhit = WireHitRange[index].second;
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
  std::cout<<" is OK. nhits "<<nhts<<"\n";
*/
	if (WireHitRange.empty()||(fFirstWire == fLastWire)){
	  mf::LogWarning("CC")<<"No hits in "<<tpcid<<" plane "<<plane;
	  continue;
	}
	else {
	  LOG_DEBUG("CC")
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
        float tickToDist = detprop->DriftVelocity(detprop->Efield(),larprop->Temperature());
        tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
        fScaleF = tickToDist / wirePitch;
        // convert Large Angle Cluster crawling cut to a slope cut
        if(fLAClusAngleCut > 0) 
          fLAClusSlopeCut = std::tan(3.142 * fLAClusAngleCut / 180.) / fScaleF;
        fMaxTime = detprop->NumberTimeSamples();
        fNumWires = geom->Nwires(plane, tpc, cstat);
        
        // look for clusters
        ClusterLoop();
        // analyze hits
//        AnalyzeHits();
      } // plane
      if(fVertex3DCut > 0) {
        // Match vertices in 3 planes
        VtxMatch(tpcid);
        Vtx3ClusterMatch(tpcid);
				if(fHammerCluster) FindHammerClusters();
        // split clusters using 3D vertices
        Vtx3ClusterSplit(tpcid);
      }
      if(fDebugPlane >= 0) {
        mf::LogVerbatim("CC")<<"Clustering done in TPC ";
        PrintClusters();
      }
    } // for all tpcs
    
    // remove the hits that have become obsolete
    RemoveObsoleteHits();
    
    WireHitRange.clear(); 
    fcl2hits.clear();
    chifits.clear();
    hitNear.clear();
    chgNear.clear();
    
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
        unsigned int span = 3;
        if(fMinHits[pass] < span) span = fMinHits[pass];
        for(auto iwire = fLastWire; iwire > fFirstWire + span; iwire--) {
          auto index = iwire - fFirstWire;
          // skip bad wires or no hits on the wire
          if(WireHitRange[index].first < 0) continue;
          auto ifirsthit = (unsigned int)WireHitRange[index].first;
          auto ilasthit = (unsigned int)WireHitRange[index].second;
          for(auto ihit = ifirsthit; ihit < ilasthit; ++ihit) {
            bool ClusterAdded = false;
            recob::Hit const& hit = fHits[ihit];
            // skip used hits
            if(ihit > fHits.size()-1) {
              mf::LogError("CC")<<"ClusterLoop bad ihit "<<ihit;
              return;
            }
            // skip used and obsolete hits
            if(inClus[ihit] != 0) continue;
            if(fDebugPlane == (short)plane && (short)iwire == fDebugWire && fDebugHit > 0)
              prt = std::abs(hit.PeakTime() - fDebugHit) < 20;
            // Check for a hit signal on the next DS wire
            bool SigOK = ChkSignal(iwire + 1, hit.PeakTime() - 10,iwire + 1, hit.PeakTime() + 10);
            // Don't start a seed cluster if there is a hit signal DS. 
            // This is an indicator that we might be trying
            // to start a cluster just US of shower blob.
            if(prt) {
              if(SigOK) mf::LogVerbatim()<<"Seed hit is near DS hits";
              if(hit.Multiplicity() > 1) mf::LogVerbatim()<<"Seed hit has multiplicity "<<hit.Multiplicity();
            }
            if(SigOK && hit.Multiplicity() > 1) continue;
            if((iwire - span + 1) < fFirstWire) continue;
            unsigned short jwire = iwire - span + 1;
            if(jwire > fLastWire) continue;
            unsigned short jindx = jwire - fFirstWire;
            if(WireHitRange[jindx].first < 0) continue;
            // Find the hit on wire jwire that best matches a line between
            // a nearby vertex and hit ihit. No constraint if useHit < 0
            unsigned int useHit = 0;
            bool doConstrain = false;
            VtxConstraint(iwire, ihit, jwire, useHit, doConstrain);
            unsigned int jfirsthit = (unsigned int)WireHitRange[jindx].first;
            unsigned int jlasthit = (unsigned int)WireHitRange[jindx].second;
            if(jfirsthit > fHits.size()-1 || jfirsthit > fHits.size()-1) {
              mf::LogError("CC")<<"ClusterLoop jwire "<<jwire<<" bad firsthit "<<jfirsthit<<" lasthit "<<jlasthit<<" fhits size "<<fHits.size();
              exit(1);
            }
            for(unsigned int jhit = jfirsthit; jhit < jlasthit; ++jhit) {
              if(jhit > fHits.size()-1) {
                mf::LogError("CC")<<"ClusterLoop bad jhit "<<jhit<<" firsthit "<<jfirsthit<<" lasthit "<<jlasthit<<" fhits size"<<fHits.size();
                exit(1);
              }
              recob::Hit const& other_hit = fHits[jhit];
              // skip used and obsolete hits
              if(inClus[jhit] != 0) continue;
              // Vertex constraint
              if(doConstrain && jhit != useHit) continue;
              // start a cluster with these two hits
              fcl2hits.clear();
              chifits.clear();
              hitNear.clear();
              chgNear.clear();
              fAveChg = -1.;
              clEndChg = -1.;
              clStopCode = 0;
              clProcCode = pass;
              fcl2hits.push_back(ihit);
              chifits.push_back(0.);
              hitNear.push_back(0);
              chgNear.push_back(0); // This will be defined if the cluster survives the cuts
              fcl2hits.push_back(jhit);
              chifits.push_back(0.);
              hitNear.push_back(0);
              chgNear.push_back(0); // This will be defined if the cluster survives the cuts
              clLA = false;
              clpar[0] = other_hit.PeakTime();
              clpar[1] = (hit.PeakTime() - other_hit.PeakTime()) / (iwire - jwire);
              clpar[2] = fHits[jhit].WireID().Wire;
/*
              // check for hit width consistency. Large angle clusters should have
              // large rms hits or smaller rms hits if they part of a multiplet
              if(std::abs(clpar[1]) > 10*geom->WirePitch(hit.View())/0.3) {
                if(pass == 0) continue;
                if(prt) mf::LogVerbatim("CC")<<" other_hit "<<other_hit.WireID().Wire<<":"<<(int)other_hit.PeakTime()<<" slope "<<clpar[1]<<" hit.RMS "<<hit.RMS()<<" other_hit multiplicity "<<other_hit.Multiplicity();
                if(hit.RMS() < 4 && hit.Multiplicity() < 3) continue;
                if(other_hit.RMS() < 4 && other_hit.Multiplicity() < 3) continue;
              } // std::abs(clpar[1]) > 5
*/
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
                // store the cluster
                if(!TmpStore()) continue;
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

      // Merge overlapping clusters
      if(fMergeOverlapAngCut > 0) MergeOverlap();
      // Check the DS end of clusters
      if(fChkClusterDS) ChkClusterDS();
      // split clusters using vertices
      if(fVtxClusterSplit) VtxClusterSplit();
      // Look for 2D vertices with star topology - short, back-to-back clusters
      if(fFindStarVertices) FindStarVertices();

      if(fDebugPlane == (short)plane) {
        mf::LogVerbatim("CC")<<"Clustering done in plane "<<plane;
        PrintClusters();
      }
      
      CheckHitClusterAssociations();
    
  } // ClusterLoop()

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeOverlap()
    {
      // Tries to merge overlapping clusters schematically shown below. The minimal condition is that both
      // clusters have a length of at least minLen wires with minOvrLap of the minLen wires in the overlap region
      //      End  Begin
      // icl  ------
      // jcl     ------
      //         End  Begin
      // This can occur when tracking cosmic rays through delta ray showers.
      // If successfull the hits in clusters icl and jcl are merged into one long cluster
      // and a short cluster if there are sufficient remaining hits.
      // This routine changes the "pass" variable to define cuts and should NOT be used inside any pass loops
      
      unsigned short icl, jcl;
      
      prt = (fDebugWire == 666);
      if(prt) mf::LogVerbatim("CC")<<"MergeOverlap check. clCTP "<<clCTP;
      
      unsigned short minLen = 10;
      unsigned short minOvrLap = 3;

      unsigned short tclsize = tcl.size();
      unsigned short overlapSize, ii, indx, bWire, eWire;
      unsigned int iht;
      for(icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].BeginVtx >= 0) continue;
        if(tcl[icl].tclhits.size() < minLen) continue;
        bool gotone = false;
        for(jcl = 0; jcl < tclsize; ++jcl) {
          if(icl == jcl) continue;
          if(tcl[jcl].ID < 0) continue;
          if(tcl[jcl].CTP != clCTP) continue;
          if(tcl[jcl].EndVtx >= 0) continue;
          if(tcl[jcl].tclhits.size() < minLen) continue;
          // icl Begin is not far enough DS from the end of jcl
          if(tcl[icl].BeginWir < tcl[jcl].EndWir + minOvrLap) continue;
          // and it doesn't end within the wire boundaries of jcl
          if(tcl[icl].BeginWir > tcl[jcl].BeginWir - minOvrLap) continue;
          // jcl End isn't far enough US from the end of icl
          if(tcl[jcl].EndWir < tcl[icl].EndWir + minOvrLap) continue;
          // require roughly similar angles at the overlap region
          if(std::abs(tcl[icl].BeginAng - tcl[jcl].EndAng) > 0.2) continue;
            // and similar time
          if(std::abs(tcl[icl].BeginTim - tcl[jcl].EndTim) > 50) continue;
          overlapSize = tcl[icl].BeginWir - tcl[jcl].EndWir + 1;
          bWire = tcl[jcl].EndWir;
          eWire = tcl[icl].BeginWir;
          if(prt) mf::LogVerbatim("CC")<<" Candidate icl ID "<<tcl[icl].ID<<" "<<tcl[icl].EndWir<<"-"<<tcl[icl].BeginWir<<" jcl ID "<<tcl[jcl].ID<<" "<<tcl[jcl].EndWir<<"-"<<tcl[jcl].BeginWir<<" overlapSize "<<overlapSize<<" bWire "<<bWire<<" eWire "<<eWire;
          std::vector<int> iclHit(overlapSize, -1);
          std::vector<int> jclHit(overlapSize, -1);
          // enter the hit times in the overlap region into the two vectors
          for(ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            iht = tcl[icl].tclhits[ii];
            if(fHits[iht].WireID().Wire < bWire) break;
            indx = fHits[iht].WireID().Wire - bWire;
            if(indx > overlapSize - 1) {
              mf::LogError("CC")<<"MergeOverlap: icl indx error "<<indx<<" overlapSize "<<overlapSize;
              return;
            } // bad indx
            iclHit[indx] = iht;
          } // ii
          for(ii = 0; ii < tcl[jcl].tclhits.size(); ++ii) {
            iht = tcl[jcl].tclhits[tcl[jcl].tclhits.size() - ii - 1];
            if(fHits[iht].WireID().Wire > eWire) break;
            indx = fHits[iht].WireID().Wire - bWire;
            if(indx > overlapSize - 1) {
              mf::LogError("CC")<<"MergeOverlap: jcl indx error "<<indx<<" overlapSize "<<overlapSize;
              return;
            } // bad indx
            jclHit[indx] = iht;
          } // ii
          if(prt) {
            mf::LogVerbatim("CC")<<"overlap vectors";
            for(ii = 0; ii < iclHit.size(); ++ii) mf::LogVerbatim("CC")<<bWire+ii<<" "<<iclHit[ii]<<" "<<jclHit[ii];
          }
          // draw a line between a hit on icl just US of the overlap region and a hit on jcl just DS of the overlap region
          // use iclpar to hold this information
          float iclpar[3];
          iclpar[0] = -1;
          for(ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            iht = tcl[icl].tclhits[ii];
            if(fHits[iht].WireID().Wire < bWire) {
              iclpar[0] = fHits[iht].PeakTime();
              iclpar[2] = fHits[iht].WireID().Wire;
              if(prt) mf::LogVerbatim("CC")<<"Use hit "<<iht<<" W T "<<fHits[iht].WireID().Wire<<" "<<fHits[iht].PeakTime();
              break;
            }
          } // ii
          if(iclpar[0] < 0) continue;
          // now find a hit on the jcl cluster just DS of the overlap region
          unsigned short hitFit = 0;
          short nHitFit = -1;
          for(ii = 0; ii < tcl[jcl].tclhits.size(); ++ii) {
            iht = tcl[jcl].tclhits[tcl[jcl].tclhits.size() - ii - 1];
            if(fHits[iht].WireID().Wire >= eWire) {
              iclpar[1] = (fHits[iht].PeakTime() - iclpar[0]) / (fHits[iht].WireID().Wire - iclpar[2]);
              if(prt) mf::LogVerbatim("CC")<<" and hit "<<iht<<" W T "<<fHits[iht].WireID().Wire<<" "<<fHits[iht].PeakTime()<<" slope "<<iclpar[1];
              hitFit = iht;
              nHitFit = tcl[jcl].tclhits.size() - ii;
              if(nHitFit > 4) nHitFit = 4;
              break;
            }
          } // ii
          if(nHitFit < 0) continue;
          // fit the hits just DS of the overlap region on jcl to make a tighter angle cut
          FitClusterMid(jcl, hitFit, -nHitFit);
          if(clChisq > 5) continue;
          // fit parameters are stored in clpar
          // make a tight angle cut
          float dth = std::abs(atan(fScaleF * iclpar[1]) - atan(fScaleF * clpar[1]));
          if(prt) mf::LogVerbatim("CC")<<"iclpar "<<iclpar[1]<<" clpar "<<clpar[1]<<" dth "<<dth;
          if(dth > fMergeOverlapAngCut) continue;
          // prepare to make a new cluster
          TmpGet(jcl);
          // resize it
          unsigned short jclNewSize;
          for(jclNewSize = 0; jclNewSize < fcl2hits.size(); ++jclNewSize) {
            iht = fcl2hits[jclNewSize];
            if(fHits[iht].WireID().Wire <= eWire) break;
          } // jclNewSize
          if(prt) {
            mf::LogVerbatim("CC")<<"jcl old size "<<fcl2hits.size()<<" newSize "<<jclNewSize;
            iht = fcl2hits[fcl2hits.size()-1];
            unsigned short iiht = fcl2hits[jclNewSize-1];
            mf::LogVerbatim("CC")<<"jcl old last wire "<<fHits[iht].WireID().Wire<<" After resize last wire "<<fHits[iiht].WireID().Wire;
          }
          fcl2hits.resize(jclNewSize);
          // now add hits in the overlap region. Window for assigning hits (ala AddHit) inflated somewhat
          float hiterr = 6 * AngleFactor(clpar[1]) * fHitErrFac * fHits[hitFit].RMS();
          if(prt) mf::LogVerbatim("CC")<<"hiterr "<<hiterr;
          // clobber icl and jcl
          MakeClusterObsolete(icl);
          MakeClusterObsolete(jcl);
          unsigned short wire;
          float prtime, best;
          unsigned int iHit, jHit, goodHit;
          unsigned short nadd = 0;
          for(ii = 0; ii < iclHit.size(); ++ii) {
            indx = iclHit.size() - ii - 1;
            wire = bWire + indx;
            prtime = iclpar[0] + iclpar[1] * (wire - iclpar[2]);
            // no hit on this wire TODO Look for the missing hit
            if(iclHit[indx] < 0 && jclHit[indx] < 0) continue;
            // Check for two hits in a multiplet and merge them
            if(iclHit[indx] >= 0 && jclHit[indx] >= 0) {
              iHit = iclHit[indx];
              jHit = jclHit[indx];
              if(areInSameMultiplet(fHits[iHit], fHits[jHit])) {
                bool didMerge = false;
                MergeHits(iHit, didMerge);
                fcl2hits.push_back(iHit);
                iclHit[indx] = -1;
                jclHit[indx] = -1;
                ++nadd;
              }
            } // iclHit[indx] >= 0 && jclHit[indx] >= 0
            else {
              // one hit on icl or jcl
              // check for a good hit on icl
              goodHit = 0;
              best = hiterr;
              if(iclHit[indx] >= 0 && std::abs(fHits[iclHit[indx]].PeakTime() - prtime) < best) {
                goodHit = iclHit[indx];
                best = fabs(fHits[iclHit[indx]].PeakTime() - prtime);
              }
              // check for a better one on jcl
              if(jclHit[indx] >= 0 && std::abs(fHits[jclHit[indx]].PeakTime() - prtime) < best) {
                goodHit = jclHit[indx];
                best = fabs(fHits[jclHit[indx]].PeakTime() - prtime);
              }
              if(best < hiterr) {
                fcl2hits.push_back(goodHit);
                // flag hits used
                if(iclHit[indx] >= 0 && (unsigned int)iclHit[indx] == goodHit) iclHit[indx] = -1;
                if(jclHit[indx] >= 0 && (unsigned int)jclHit[indx] == goodHit) jclHit[indx] = -1;
                ++nadd;
              } // best < hiterr
            } // one hit on icl or jcl
          } // ii
          // now paste in the icl hits that are US of the overlap region
          for(ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            iht = tcl[icl].tclhits[ii];
            if((unsigned short)fHits[iht].WireID().Wire >= bWire) continue;
            fcl2hits.push_back(iht);
          }
          clBeginSlp = tcl[jcl].BeginSlp;
          clBeginSlpErr = tcl[jcl].BeginSlpErr;
          clBeginAng = tcl[jcl].BeginAng;
          clBeginWir = tcl[jcl].BeginWir;
          clBeginTim = tcl[jcl].BeginTim;
          clBeginChg = tcl[jcl].BeginChg;
          clBeginChgNear = tcl[jcl].BeginChgNear;
          // End info from icl
          clEndSlp = tcl[icl].EndSlp;
          clEndSlpErr = tcl[icl].EndSlpErr;
          clEndAng = tcl[icl].EndAng;
          clEndWir = tcl[icl].EndWir;
          clEndTim = tcl[icl].EndTim;
          clEndChg = tcl[icl].EndChg;
          clEndChgNear = tcl[icl].EndChgNear;
          clStopCode = tcl[icl].StopCode;
          clProcCode = tcl[icl].ProcCode + 500;
          if(!TmpStore()) {
            // Merged cluster is fubar. Try to recover
            RestoreObsoleteCluster(icl);
            RestoreObsoleteCluster(jcl);
            continue;
          }
          tcl[tcl.size() - 1].BeginVtx = tcl[jcl].BeginVtx;
          tcl[tcl.size() - 1].EndVtx = tcl[icl].BeginVtx;
          gotone = true;
          // after this point any failure should result in a jcl loop break
          if(prt) mf::LogVerbatim("CC")<<"MergeOverlap new long cluster ID "<<tcl[tcl.size()-1].ID<<" in clCTP "<<clCTP;
          // try to make a new cluster from the remnants in iclHit and jclHit
          // put the valid hit indices into the iclHit vector
          for(ii = 0; ii < iclHit.size(); ++ii) if(iclHit[ii] < 0) iclHit[ii] = jclHit[ii];
          if(prt) {
            mf::LogVerbatim("CC")<<"updated iclHit";
            for(ii = 0; ii < iclHit.size(); ++ii) mf::LogVerbatim("CC")<<bWire+ii<<" "<<iclHit[ii];
          }
          fcl2hits.clear();
          unsigned short nmiss = 0;
          for(ii = 0; ii < iclHit.size(); ++ii) {
            indx = iclHit.size() - ii - 1;
            if(iclHit[indx] < 0) {
              ++nmiss;
              nadd = 0;
              continue;
            }
            if(inClus[iclHit[indx]] != 0) {
              std::cout<<"Hit not free "<<iclHit[indx]<<" on wire "<<fHits[iclHit[indx]].WireID().Wire<<"\n";
              continue;
            }
            fcl2hits.push_back((unsigned int)iclHit[indx]);
            ++nadd;
            if(nmiss > 0) continue;
            nmiss = 0;
          } // ii
          if(fcl2hits.size() < 3) break;
          // check the hitfraction
          iht = fcl2hits.size() - 1;
          float dw = fHits[fcl2hits[0]].WireID().Wire - fHits[fcl2hits[iht]].WireID().Wire + 1;
          float hitFrac = (float)fcl2hits.size() / dw;
          if(hitFrac < fMinHitFrac) break;
          // do the End fit
          pass = fNumPass - 1;
          FitCluster();
          if(clChisq > fChiCut[pass]) break;
          clEndSlp = clpar[1];
          clEndSlpErr = clparerr[1];
          clEndAng = std::atan(fScaleF * clEndSlp);
          clEndWir = fHits[fcl2hits[iht]].WireID().Wire;
          clEndTim = fHits[fcl2hits[iht]].PeakTime();
          // find the charge at the end
          FitClusterChg();
          clEndChg = fAveChg;
          clEndChgNear = 1; // This is a fake but not a bad fake
          if(fcl2hits.size() > fMaxHitsFit[pass]) {
            // do some temporary damage to fit the Begin end
            std::vector<unsigned int> tfcl2hits = fcl2hits;
            fcl2hits.resize(fMaxHitsFit[pass]);
            FitCluster();
            if(clChisq > fChiCut[pass]) break;
            clBeginSlp = clpar[1];
            clBeginSlpErr = clparerr[1];
            clBeginAng = std::atan(fScaleF * clBeginSlp);
            clBeginWir = fHits[fcl2hits[0]].WireID().Wire;
            clBeginTim = fHits[fcl2hits[0]].PeakTime();
            // find the charge at the end
            FitClusterChg();
            clBeginChg = fAveChg;
            clBeginChgNear = 1; // This is a fake but not a bad fake
            fcl2hits = tfcl2hits;
          }
          clProcCode = 666;
          if(!TmpStore()) continue;
          if(prt) mf::LogVerbatim("CC")<<"MergeOverlap new runt cluster ID "<<tcl[tcl.size()-1].ID<<" in clCTP "<<clCTP;
          if(gotone) break;
        } // jcl
      } // icl
      
      
    } // MergeOverlap()

//////////////////////////////////////////
  void ClusterCrawlerAlg::MakeClusterObsolete(unsigned short icl) {
    short& ID = tcl[icl].ID;
    if (ID <= 0) {
      mf::LogError("CC")<<"Trying to make already-obsolete cluster obsolete ID = "<<ID;
      return; // already obsolete
    }
    ID = -ID;                     // mark the cluster as obsolete
    
    // release the hits
    for(unsigned int iht = 0; iht < tcl[icl].tclhits.size(); ++iht) inClus[tcl[icl].tclhits[iht]] = 0;

  } // ClusterCrawlerAlg::MakeClusterObsolete()

//////////////////////////////////////////
  void ClusterCrawlerAlg::RestoreObsoleteCluster(unsigned short icl) {
    short& ID = tcl[icl].ID;
    if(ID > 0) {
      mf::LogError("CC")<<"Trying to restore non-obsolete cluster ID = "<<ID;
      return;
    }
    ID = -ID;
    
    for(unsigned short iht = 0; iht < tcl[icl].tclhits.size(); ++iht) inClus[tcl[icl].tclhits[iht]] = ID;
    
  } // RestoreObsoleteCluster()

//////////////////////////////////////////
  void ClusterCrawlerAlg::CheckHitClusterAssociations()
  {
    // check hit - cluster associations
    
    if(fHits.size() != inClus.size()) {
      mf::LogError("CC")<<"CHCA: Sizes wrong "<<fHits.size()<<" "<<inClus.size();
      return;
    }
    
    unsigned int iht, nErr = 0;
    short clID;
    
    // check cluster -> hit association
    for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if(tcl[icl].ID < 0) continue;
      clID = tcl[icl].ID;
      for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
        iht = tcl[icl].tclhits[ii];
        if(iht > fHits.size() - 1) {
          mf::LogError("CC")<<"CHCA: Bad tclhits index "<<iht<<" fHits size "<<fHits.size();
          return;
        } // iht > fHits.size() - 1
        if(inClus[iht] != clID) {
          mf::LogError("CC")<<"CHCA: Bad cluster -> hit association. clID "<<clID<<" hit inClus "<<inClus[iht]<<" ProcCode "<<tcl[icl].ProcCode<<" CTP "<<tcl[icl].CTP;
          ++nErr;
          if(nErr > 10) return;
        }
      } // ii
    } // icl
    
    // check hit -> cluster association
    unsigned short icl;
    for(iht = 0; iht < fHits.size(); ++iht) {
      if(inClus[iht] <= 0) continue;
      icl = inClus[iht] - 1;
      // see if the cluster is obsolete
      if(tcl[icl].ID < 0) {
        mf::LogError("CC")<<"CHCA: Hit associated with an obsolete cluster. hit W:T "<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()
        <<" tcl[icl].ID "<<tcl[icl].ID;
        ++nErr;
        if(nErr > 10) return;
      }
      if (std::find(tcl[icl].tclhits.begin(), tcl[icl].tclhits.end(), iht) == tcl[icl].tclhits.end()) {
        mf::LogError("CC")<<"CHCA: Bad hit -> cluster association. hit index "<<iht
          <<" W:T "<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" inClus "<<inClus[iht];
        ++nErr;
        if(nErr > 10) return;
      }
    } // iht
    
  } // CheckHitClusterAssociations()
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::RemoveObsoleteHits() {
    
    unsigned int destHit = 0;
    
    if(fHits.size() != inClus.size()) {
      mf::LogError("CC")<<"RemoveObsoleteHits size mis-match "<<fHits.size()<<" "<<inClus.size();
      return;
    }
    
    unsigned short icl;
    for(unsigned int srcHit = 0; srcHit < fHits.size(); ++srcHit) {
      if(inClus[srcHit] < 0) continue;
      if(srcHit != destHit) {
        fHits[destHit] = std::move(fHits[srcHit]);
        inClus[destHit] = inClus[srcHit];
        if(inClus[destHit] > 0) {
          // hit is in a cluster. Find it and change the index
          icl = inClus[destHit] - 1;
          auto& hits = tcl[icl].tclhits;
          auto iHitIndex = std::find(hits.begin(), hits.end(), srcHit);
          if (iHitIndex == hits.end()) {
            mf::LogError("CC")<< "RemoveObsoleteHits: Hit #" << srcHit << " not found in cluster ID "<< inClus[destHit];
          } else {
            *iHitIndex = destHit; // update the index
          }
        } // inClus[destHit] > 0
      }
      ++destHit;
    } // srcHit
    
    fHits.resize(destHit);
    inClus.resize(destHit);
    
  } // RemoveObsoleteHits()
  
//////////////////////////////////////////
    void ClusterCrawlerAlg::ChkClusterDS() {
      // Try to extend clusters DS by a few wires. 
      // DS hits may not have been  included in a cluster if they have high 
      // multiplicity or high charge. 
      // Ref ClusterLoop cuts for starting a seed cluster.

      prt = (fDebugPlane == 3);

      if(prt) mf::LogVerbatim("CC")<<"ChkClusterDS clCTP "<<clCTP;

      const unsigned short tclsize = tcl.size();
      bool didMerge;
      unsigned short icl, ii, nhm;
      unsigned int iht;
      
      // first merge any hits on the DS end of clusters
      for(icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore clusters that have a Begin vertex
        if(tcl[icl].BeginVtx >= 0) continue;
        // check the first few hits
        nhm = 0;
        for(ii = 0; ii < 3; ++ii) {
          iht = fcl2hits[ii];
          if(fHits[iht].Multiplicity() > 1) {
            MergeHits(iht, didMerge);
            if(didMerge) ++nhm;
          }
        } // ii
        if(nhm > 0) {
          // update the Begin parameters in-place
          FitClusterMid(icl, 0, 3);
          tcl[icl].BeginTim = clpar[0];
          tcl[icl].BeginSlp = clpar[1];
          tcl[icl].BeginAng = atan(fScaleF * clpar[1]);
          tcl[icl].BeginSlpErr = clparerr[1];
          tcl[icl].BeginChg = fAveChg;
          tcl[icl].ProcCode += 5000;
          if(prt) mf::LogVerbatim("CC")<<"ChkClusterDS: Merge hits on cluster "<<tcl[icl].ID;
        } // nhm > 0
      } // icl

      float thhits, prevth, hitrms, rmsrat;
      bool ratOK;
      std::vector<unsigned int> dshits;
      for(unsigned short icl = 0; icl < tclsize; ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // ignore clusters that have a Begin vertex
        if(tcl[icl].BeginVtx >= 0) continue;
        // ignore clusters with lots of nearby charge
        if(tcl[icl].BeginChgNear > fChgNearCut) continue;
        // find the angle using the first 2 hits
        unsigned int ih0 = tcl[icl].tclhits[1];
        unsigned int ih1 = tcl[icl].tclhits[0];
        const float slp = (fHits[ih1].PeakTime()    - fHits[ih0].PeakTime()) / 
                          (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
        prevth = std::atan(fScaleF * slp);
        // move the "origin" to the first hit
        ih0 = ih1;
        unsigned short wire = fHits[ih0].WireID().Wire;
        hitrms = fHits[ih0].RMS();
        float time0 = fHits[ih0].PeakTime();
        float prtime;
        dshits.clear();
        // follow DS a few wires. Stop if any encountered
        // hit is associated with a cluster
        for(unsigned short ii = 0; ii < 2; ++ii) {
          ++wire;
          if(wire > fLastWire) break;
          prtime = time0 + slp;
          if(prt) mf::LogVerbatim("CC")<<"ChkClusterDS: Try to extend "
            <<tcl[icl].ID<<" to W:T "<<wire<<" hitrms "<<hitrms<<" prevth "<<prevth<<" prtime "<<(int)prtime;
          unsigned short indx = wire - fFirstWire;
          // stop if no hits on this wire
          if(WireHitRange[indx].first == -2) break;
        //  if(WireHitRange[index].first == -1) break; // FIXME: trouble; why??
          unsigned int firsthit = WireHitRange[indx].first;
          unsigned int lasthit = WireHitRange[indx].second;
          bool hitAdded = false;
          for(ih1 = firsthit; ih1 < lasthit; ++ih1) {
            if(inClus[ih1] != 0) continue;
            if(prtime < fHits[ih1].PeakTimeMinusRMS(5)) continue;
            if(prtime > fHits[ih1].PeakTimePlusRMS(5)) continue;
            const float slp = (fHits[ih1].PeakTime() - fHits[ih0].PeakTime()) /
                              (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
            thhits = std::atan(fScaleF * slp);
            if(prt) mf::LogVerbatim("CC")<<" theta "<<thhits;
            if(fabs(thhits - prevth) > 0.4) continue;
            rmsrat = fHits[ih1].RMS() / hitrms;
            ratOK = rmsrat > 0.3 && rmsrat < 3;
            // require small angle and not wildly different width compared
            // to the first hit in the cluster
            // TODO handle hit multiplets here
            if(ratOK && std::abs(thhits - prevth) < 0.3) {
              dshits.push_back(ih1);
              hitAdded = true;
              prevth = thhits;
              ih0 = ih1;
  if(prt) mf::LogVerbatim("CC")<<" Add hit "<<fHits[ih1].WireID().Wire
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
            unsigned int iht = tcl[icl].tclhits[ii];
            inClus[iht] = 0;
            fcl2hits.push_back(iht);
          }
          clProcCode += 5900;
          pass = fNumPass - 1;
          FitClusterChg();
          clBeginChg = fAveChg;
          // declare the old one obsolete
          MakeClusterObsolete(icl);
          // add the new one
          if(!TmpStore()) {
            mf::LogError("CC")<<"ChkClusterDS TmpStore failed while extending cluster ID "<<tcl[icl].ID;
            continue;
          }
          const size_t newcl = tcl.size() -1;
  if(prt) mf::LogVerbatim("CC")<<" Store "<<newcl;
          tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
          tcl[newcl].EndVtx = tcl[icl].EndVtx;
        } // dshits.size() > 0
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
    
    unsigned int iht = fcl2hits.size() - 1;
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
    unsigned int iht = fcl2hits.size() - 1;
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
        unsigned short iwire, unsigned int ihit, unsigned short jwire,
        unsigned int& useHit, bool& doConstrain)
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
      unsigned int jfirsthit = WireHitRange[jindx].first;
      unsigned int jlasthit = WireHitRange[jindx].second;
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != clCTP) continue;
        // vertex must be US of the cluster
        if(vtx[iv].Wire > jwire) continue;
        // but not too far US
        if(vtx[iv].Wire < jwire - 10) continue;
        recob::Hit const& hit = fHits[ihit];
        clpar[0] = hit.PeakTime();
        clpar[1] = (vtx[iv].Time - hit.PeakTime()) / (vtx[iv].Wire - iwire);
        clpar[2] = hit.WireID().Wire;
        float prtime = clpar[0] + clpar[1] * (jwire - iwire);
        for(unsigned int jhit = jfirsthit; jhit < jlasthit; ++jhit) {
          if(inClus[jhit] != 0) continue;
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
      
      vtxprt = (fDebugPlane >= 0 && fDebugWire == 5555);

  if(vtxprt) mf::LogVerbatim("CC")<<"VtxClusterSplit ";

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
          if(vtx[ivx].NClusters == 0) continue;
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
          if(vtxprt) mf::LogVerbatim("CC")<<" Chk cluster ID "<<tcl[icl].ID<<" with vertex "<<ivx;
          short ihvx = -99;
          // nSplit is the index of the hit in the cluster where we will
          // split it if all requirements are met
          unsigned short nSplit = 0;
          unsigned short nLop = 0;
          for(unsigned short ii = tcl[icl].tclhits.size()-1; ii > 0 ; --ii) {
            unsigned int iht = tcl[icl].tclhits[ii];
            ++nLop;
            if(fHits[iht].WireID().Wire >= vtx[ivx].Wire) {
              nSplit = ii;
  if(vtxprt) mf::LogVerbatim("CC")<<"Split cluster "<<tcl[icl].ID<<" at wire "<<fHits[iht].WireID().Wire
		<<" nSplit "<<nSplit;
              ihvx = iht;
              break;
            }
          } // ii
          // found the wire. Now make a rough time cut
          if(ihvx < 0) continue;
//          short dtime = std::abs(short(fHits[ihvx].PeakTime() - vtx[ivx].Time));
//  if(vtxprt) mf::LogVerbatim("CC")<<" Delta time "<<dtime;
          if(fabs(fHits[ihvx].PeakTime() - vtx[ivx].Time) > 10) continue;
          // check the angle between the crossing cluster icl and the
          // clusters that comprise the vertex. 
          // First decide which end of cluster icl to use to define the angle
          float iclAng = 0.;
          if(nSplit > tcl[icl].tclhits.size() / 2) {
            iclAng = tcl[icl].EndAng;
          } else {
            iclAng = tcl[icl].BeginAng;
          }
  				if(vtxprt) mf::LogVerbatim("CC")<<" iclAng "<<iclAng;
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
  if(vtxprt) mf::LogVerbatim("CC")<<"Split/Chop at pos "<<nLop;
            if(nLop < 3) {
              // lop off hits at the US end
              // Put the cluster in the local arrays
              TmpGet(icl);
              for(unsigned short ii = 0; ii < nLop; ++ii) {
                unsigned int iht = fcl2hits[fcl2hits.size()-1];
                inClus[iht] = 0;
                fcl2hits.pop_back();
              }
              // store it
              clProcCode += 1000;
              // declare this cluster obsolete
              MakeClusterObsolete(icl);
              // store the new one
              if(!TmpStore()) continue;
              unsigned short newcl = tcl.size() - 1;
              tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
              tcl[newcl].EndVtx = ivx;
            } else {
              // split the cluster into two
              // correct the split position
              ++nSplit;
              if(SplitCluster(icl, nSplit, ivx)) {
                tcl[tcl.size()-1].ProcCode += 1000;
                tcl[tcl.size()-2].ProcCode += 1000;
              }
            }
            didit = true;
          } // angOK
          if(didit) break;
        } // ivx
      } // icl
      
    } // VtxClusterSplit()

//////////////////////////////////////////
    void ClusterCrawlerAlg::MergeHits(const unsigned int theHit, bool& didMerge) {
      // Merge all unused separate hits in the multiplet of which 
      // theHit is a member into one hit (= theHit).
      // Mark the merged hits other than theHit obsolete.
      // Hits in the multiplet that are associated with an existing cluster are
      // not affected.
      // Hit multiplicity is reworked (including all the hits in the multiplet).
      // Used hits have the multiplicity and index corrected too; the local
      // index reflects the peak time.
      // Note that theHit may or may not be marked free (usually, it is not)

      didMerge = false;
      
      if(theHit > fHits.size() - 1) {
        mf::LogError("CC")<<"Bad theHit";
        return;
      }
      
      // don't bother trying to merge an already merged hit
      if(fHits[theHit].GoodnessOfFit() == 6666) return;
      
      recob::Hit const& hit = fHits[theHit];

  if(prt) mf::LogVerbatim("CC")<<"MergeHits "<<hit.WireID().Wire<<":"<<(int)hit.PeakTime()
    <<" numHits "<<hit.Multiplicity();

      // number of hits in this hit multiplet
      std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(theHit);
      
      // ensure that this is a high multiplicity hit:
      if (MultipletRange.second <= MultipletRange.first) return;
      
      // do a quick check to see how many hits are available to be merged
      unsigned short nAvailable = 0;
      for(size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
        if(jht == theHit) continue;
        if(inClus[jht] != 0) continue;
        ++nAvailable;
      } // jht
      if(nAvailable == 0) return;
      
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
        if(inClus[jht] < 0) continue;
        recob::Hit const& other_hit = fHits[jht];

  if(prt) {
    mf::LogVerbatim("CC")
    <<" P:W:T "<<plane<<":"<<other_hit.WireID().Wire<<":"<<(int)other_hit.PeakTime()
    <<" Amp "<<(int)other_hit.PeakAmplitude()
    <<" RMS "<<other_hit.RMS()
    <<" Charge "<<(int)other_hit.Integral()
    <<" inClus "<<inClus[jht];
  }
        // error checking
        if((other_hit.StartTick() != hit.StartTick())
          || (other_hit.WireID() != hit.WireID()))
        {
          mf::LogError("CC")<<"Hit multiplet ID error "
            << other_hit.WireID() << " @" << other_hit.StartTick()
              << " " << other_hit.LocalIndex() << "/" << other_hit.Multiplicity()
            << " vs. " << hit.WireID() << " @" << hit.StartTick()
              << " " << hit.LocalIndex() << "/" << hit.Multiplicity()
            ;
          return;
        }
        if (other_hit.Multiplicity() != hit.Multiplicity()) {
          mf::LogError("CC")
            << " hit #" << jht << " in the same multiplet as #" << theHit
            << " has different multiplicity!"
            << "\n hit #" << theHit << ": " << hit
            << "\n hit #" << jht << ": " << other_hit;
          return;
        }
        // hit is not used by another cluster
        if(inClus[jht] != 0) continue;
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
          if(inClus[jht] != 0) continue;
          // declare this hit obsolete
          inClus[jht] = -1;
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
        mf::LogError("CC")<<"MergeHits: bad sum";
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
        6666,                   // GoodnessOfFit (flag for merged hit)
        hit.DegreesOfFreedom(),
        hit.View(),
        hit.SignalType(),
        hit.WireID()
        );
  if(prt) {
    mf::LogVerbatim("CC")
    <<" theHit "<<fHits[theHit].WireID().Wire<<":"<<(int)aveTime
    <<" RMS "<<std::setprecision(1)<<fHits[theHit].RMS()
    <<" chg "<<(int)chgsum<<" Amp "<<(int)fHits[theHit].PeakAmplitude();
  }
      
      FixMultipletLocalIndices(MultipletRange.first, MultipletRange.second);
      didMerge = true;
      
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
          if(inClus[iHit] < 0) continue;
//          if (!isHitPresent(iHit)) continue;
          ++multiplicity;
        } // for
      } // if no valid multiplicity is given
      
      // second pass: assign the correct multiplicity
      short int local_index = 0;
      for (size_t iHit = begin; iHit < end; ++iHit) {
        if(inClus[iHit] < 0) continue;
//        if (!isHitPresent(iHit)) continue;
        
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
    mf::LogVerbatim("CC")<<"FindStarVertices";
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
          if(vtxprt) mf::LogVerbatim("CC")<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" topo5 vtx wire "<<fvw;
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
  if(vtxprt) mf::LogVerbatim("CC")<<" vertex time "<<fvt
    <<" lotime "<<lotime<<" hitime "<<hitime;
          unsigned short vw = (0.5 + fvw);
          // ensure that the vertex is near a hit on both clusters
          unsigned int pos1 = 0;
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
          newvx.WireErr = 1;
          newvx.Time = fvt;
          newvx.TimeErr = 1;
          newvx.Topo = 5;
          newvx.CTP = clCTP;
          vtx.push_back(newvx);
          unsigned short ivx = vtx.size() - 1;
          if(vtxprt) mf::LogVerbatim("CC")<<" new vertex "<<ivx<<" cluster "<<tcl[it1].ID<<" split pos "<<pos1;
          if(!SplitCluster(it1, pos1, ivx)) continue;
          tcl[tcl.size()-1].ProcCode += 1000;
          tcl[tcl.size()-2].ProcCode += 1000;
          if(vtxprt) mf::LogVerbatim("CC")<<" new vertex "<<ivx<<" cluster "<<tcl[it2].ID<<" split pos "<<pos2;
          if(!SplitCluster(it2, pos2, ivx)) continue;
          tcl[tcl.size()-1].ProcCode += 1000;
          tcl[tcl.size()-2].ProcCode += 1000;
          FitVtx(ivx);
          break;
        } // it2
      } // it1
            
      if(vtx.size() > vtxSizeIn) {
        // try to split other clusters
        VtxClusterSplit();
        // try to attach other clusters to it
        VertexCluster(vtx.size() - 1);
        FitAllVtx(clCTP);
      } // new vertex

  if(vtxprt) {
    mf::LogVerbatim("CC")<<"Vertices "<<vtx.size();
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      mf::LogVerbatim("CC")
        <<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time
        <<" NClusters "<<vtx[iv].NClusters<<" topo "<<vtx[iv].Topo;
    }
    PrintClusters();
  }
      
    } // FindStarVertices()

//////////////////////////////////////////
    void ClusterCrawlerAlg::VertexCluster(unsigned short iv)
    {
      // try to attach clusters to the specified vertex
      if(vtx[iv].NClusters == 0) return;

      short dwb, dwe, dtb, dte;
      bool sigOK;
      
      vtxprt = (fDebugPlane >= 0 && fDebugWire == 3333);
      
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != vtx[iv].CTP) continue;
        
        dwb = vtx[iv].Wire - tcl[icl].BeginWir;
        dtb = vtx[iv].Time - tcl[icl].BeginTim;
        dwe = vtx[iv].Wire - tcl[icl].EndWir;
        dte = vtx[iv].Time - tcl[icl].EndTim;
        
        float drb = dwb * dwb + dtb * dtb;
        float dre = dwe * dwe + dte * dte;
        
        bool bCloser = (drb < dre);
        
        // ignore clusters in showers
        if(bCloser) {
          if(tcl[icl].BeginChgNear > fChgNearCut) continue;
        } else {
          if(tcl[icl].EndChgNear > fChgNearCut) continue;
        }
        
        if(vtxprt) mf::LogVerbatim("CC")<<"VertexCluster: Try icl ID "<<tcl[icl].ID<<" w vtx "<<iv<<" dwb "<<dwb<<" dwe "<<dwe<<" drb "<<drb<<" dre "<<dre<<" Begin closer? "<<bCloser;
        
        if(tcl[icl].BeginVtx < 0 && bCloser && dwb > -3 && dwb < 3 && tcl[icl].EndVtx != iv) {
          sigOK = ChkSignal(tcl[icl].BeginWir, tcl[icl].BeginTim, vtx[iv].Wire, vtx[iv].Time);
          if(vtxprt) mf::LogVerbatim("CC")<<" Attach cluster Begin to vtx? "<<iv<<" sigOK "<<sigOK;
          if(sigOK) {
            if(vtxprt) mf::LogVerbatim("CC")<<" check ClusterVertexChi "<<ClusterVertexChi(icl, 0, iv);
            if(ClusterVertexChi(icl, 0, iv) < 3) {
              // do a fit and check the vertex error
              tcl[icl].BeginVtx = iv;
              FitVtx(iv);
              if(vtx[iv].ChiDOF > 5 || vtx[iv].WireErr > 2) {
                tcl[icl].BeginVtx = -99;
                FitVtx(iv);
              }
            } // good DoCA
          } // sigOK
        } // check BEGIN
        
        if(tcl[icl].EndVtx < 0 && !bCloser && dwe > -3 && dwe < 3 && tcl[icl].BeginVtx != iv) {
          sigOK = ChkSignal(tcl[icl].EndWir, tcl[icl].EndTim, vtx[iv].Wire, vtx[iv].Time);
          if(vtxprt) mf::LogVerbatim("CC")<<" Attach cluster End to vtx? "<<iv<<" sigOK "<<sigOK;
          if(sigOK) {
            if(vtxprt) mf::LogVerbatim("CC")<<" check ClusterVertexChi "<<ClusterVertexChi(icl, 1, iv);
            if(ClusterVertexChi(icl, 1, iv) < 3) {
              // do a fit and check the vertex error
              tcl[icl].EndVtx = iv;
              FitVtx(iv);
              if(vtx[iv].ChiDOF > 5 || vtx[iv].WireErr > 2) {
                tcl[icl].EndVtx = -99;
                FitVtx(iv);
              }
            } // good DoCA
          } // sigOK
        } // check END
      } // icl
    } // VertexCluster

//////////////////////////////////////////
    void ClusterCrawlerAlg::FitAllVtx(CTP_t inCTP)
      {

        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if(vtx[iv].CTP != inCTP) continue;
          FitVtx(iv);
      }

    } // FitAllVtx
  
/////////////////////////////////////////
    void ClusterCrawlerAlg::FindVertices()
    {
      // try to make 2D vertices
      
      if(tcl.size() < 2) return;

      // form vertices starting with the longest
      std::vector<CluLen> clulens;
      CluLen clulen;
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].BeginVtx >= 0 && tcl[icl].EndVtx >= 0) continue;
        clulen.index = icl;
        clulen.length = tcl[icl].tclhits.size();
        clulens.push_back(clulen);
      }
      if(clulens.size() == 0) return;
      std::sort (clulens.begin(),clulens.end(), greaterThan);
      
      float nwires = fNumWires;
      float maxtime = fMaxTime;
      
      unsigned short vtxSizeIn = vtx.size();

  vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
  if(vtxprt) {
    mf::LogVerbatim("CC")<<"FindVertices plane "<<plane<<" pass "<<pass;
    PrintClusters();
  }

      float es1 = 0, eth1 = 0, ew1 = 0, et1 = 0;
      float bs1 = 0, bth1 = 0, bw1 = 0, bt1 = 0;
      float es2 = 0, eth2 = 0, ew2 = 0, et2 = 0;
      float bs2 = 0, bth2 = 0, bw2 = 0, bt2 = 0;
      float dth = 0,  dsl = 0, fvw = 0, fvt = 0;
      // Flags set true if the cluster end has low nearby charge
//      bool ec1, bc1, ec2, bc2;
      float angcut = 0;
      unsigned short vw = 0, pass1, pass2;
      bool SigOK = false;
      for(unsigned short ii1 = 0; ii1 < clulens.size() - 1; ++ii1) {
        unsigned short it1 = clulens[ii1].index;
        es1 = tcl[it1].EndSlp;
        eth1 = tcl[it1].EndAng;
        ew1 = tcl[it1].EndWir;
        et1 = tcl[it1].EndTim;
//        ec1 = (tcl[it1].EndChgNear < fChgNearCut);
        bs1 = tcl[it1].BeginSlp;
        bth1 = tcl[it1].BeginAng;
        bw1 = tcl[it1].BeginWir;
        bt1 = tcl[it1].BeginTim;
//        bc1 = (tcl[it1].BeginChgNear < fChgNearCut);
        pass1 = tcl[it1].ProcCode - 10 *(tcl[it1].ProcCode / 10);
        for(unsigned short ii2 = ii1 + 1; ii2 < clulens.size(); ++ii2) {
          unsigned short it2 = clulens[ii2].index;
          // try to attach cluster to existing vertices at either end
          ClusterVertex(it2);
          // ignore if both clusters are short
          if(tcl[it1].tclhits.size() < 5 &&
             tcl[it2].tclhits.size() < 5) continue;
          es2 = tcl[it2].EndSlp;
          eth2 = tcl[it2].EndAng;
          ew2 = tcl[it2].EndWir;
          et2 = tcl[it2].EndTim;
//          ec2 = (tcl[it2].EndChgNear < fChgNearCut);
          bs2 = tcl[it2].BeginSlp;
          bth2 = tcl[it2].BeginAng;
          bw2 = tcl[it2].BeginWir;
          bt2 = tcl[it2].BeginTim;
//          bc2 = (tcl[it2].BeginChgNear < fChgNearCut);
          pass2 = tcl[it2].ProcCode - 10 *(tcl[it2].ProcCode / 10);
          if(pass1 < pass2) {
            angcut = fKinkAngCut[pass2];
          } else {
            angcut = fKinkAngCut[pass1];
          }
  // topo 1: check for vtx US of the ends of both clusters
          dth = fabs(eth1 - eth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0 && dth > 0.1) {
//            if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0 && dth > angcut && ec1 && ec2) {
            dsl = es2 - es1;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            fvw = (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
            // vertex wire in the detector?
            if(fvw > 0. && fvw < nwires) {
              // require vtx in the range of wires with hits AND
              vw = (0.5 + fvw);
              // vtx US of both clusters AND
              // vtx not too far US of both clusters
              if(vw >= fFirstWire - 1 &&
                 fvw <= ew1 + 1    && fvw <= ew2 + 1 &&
                 fvw > (ew1 - 10)  && fvw > (ew2 - 10) ) {
                fvt = et1 + (fvw - ew1) * es1;
                if(vtxprt) mf::LogVerbatim("CC")<<"Chk clusters topo1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<"  vtx wire "<<vw<<" time "<<(int)fvt<<" dth "<<dth;
                if(fvt > 0 && fvt < maxtime) {
                  // Vertex wire US of cluster ends and time in the detector.
                  // Check for signal at the vertex position and adjust the vertex by 1 wire
                  // if necessary
                  SigOK = ChkSignal(vw, fvt, vw, fvt);
                  if(!SigOK) {
                    fvw += 1.;
                    vw = (0.5 + fvw);
                    SigOK = ChkSignal(vw, fvt, vw, fvt);
                  }
                  // Check this against existing vertices and update
                  if(SigOK) ChkVertex(fvw, fvt, it1, it2, 1);
                } // fvt in detector
              } // vw topo 1 check
            } // fvw in detector
          } // topo 1
  // topo 2: check for vtx US of it1 and DS of it2
          dth = std::abs(eth1 - bth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0 && dth > angcut) {
//            if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0 && dth > angcut && ec1 && bc2) {
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
              if(fvw <= ew1 + 2  && fvw >= bw2 - 2) {
                fvt = et1 + (vw - ew1) * es1;
                if(vtxprt) mf::LogVerbatim("CC")<<"Chk clusters topo2 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" vtx wire "<<vw<<" time "<<(int)fvt<<" dth "<<dth;
                if(fvt > 0. && fvt < maxtime) {
                  ChkVertex(fvw, fvt, it1, it2, 2);
                } // fvt in detector
              } // vw topo 2 check
            } // fvw in detector
          } // topo 2
  // topo 3: check for vtx DS of it1 and US of it2
          dth = std::abs(bth1 - eth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0 && dth > angcut) {
//            if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0 && dth > angcut && bc1 && ec2) {
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
              if(fvw <= ew2 + 2 && fvw >= bw1 - 2) {
                fvt = et2 + (fvw - ew2) * es2;
                if(vtxprt) mf::LogVerbatim("CC")<<"Chk clusters topo3 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" vtx wire "<<vw<<" time "<<(int)fvt<<" dth "<<dth;
                if(fvt > 0. && fvt < maxtime) {
                  ChkVertex(fvw, fvt, it1, it2, 3);
                } // fvt in detector
              } // vw topo 3 check
            } // fvw in detector
          } // topo 3
  // topo 4: check for vtx DS of it1 and DS of it2
          dth = std::abs(bth1 - bth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0 && dth > 0.1) {
//            if(tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0 && dth > angcut && bc1 && bc2) {
            dsl = bs2 - bs1;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            // convert to integer if within the detector (vw, vt)
            fvw = (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl;
//            if(vtxprt) mf::LogVerbatim("CC")<<"Chk clusters topo4 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" fvw "<<fvw<<" nwires "<<nwires;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              vw = (0.5 + fvw);
              // require vtx in the range of wires with hits AND vtx DS of both clusters AND
              // vtx not too far DS of both clusters
              // Cuts are looser here to try to pick up low energy decay daughter tracks that may be going backward
              if(fvw <= fLastWire + 1 &&
                 fvw >= bw1 - 4 && fvw >= bw2 - 4 &&
                 fvw <  bw2 + 10 && fvw <  bw1 + 10) {
                fvt = bt1 + (fvw - bw1) * bs1;
                if(vtxprt) mf::LogVerbatim("CC")<<"Chk clusters topo4 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" vtx wire "<<vw<<" time "<<fvt<<" dth "<<dth;
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check for signal at the vertex position and adjust the vertex by 1 wire
                  // if necessary
                  SigOK = ChkSignal(vw, fvt, vw, fvt);
                  if(!SigOK) {
                    fvw -= 1.;
                    vw = (0.5 + fvw);
                    SigOK = ChkSignal(vw, fvt, vw, fvt);
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

      if(vtx.size() > vtxSizeIn) FitAllVtx(clCTP);
      
      // "delete" any vertices that have only one cluster
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(vtx[ivx].CTP != clCTP) continue;
        if(vtx[ivx].NClusters == 1) {
          vtx[ivx].NClusters = 0;
          for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
            if(tcl[icl].CTP != clCTP) continue;
            if(tcl[icl].BeginVtx == ivx) {
              tcl[icl].BeginVtx = -99;
              break;
            }
            if(tcl[icl].EndVtx == ivx) {
              tcl[icl].EndVtx = -99;
              break;
            }
          } // icl
        } // vtx[ivx].NClusters == 1
      } // ivx
      
      
    } // FindVertices()

/////////////////////////////////////////
    void ClusterCrawlerAlg::ClusterVertex(unsigned short it)
    {
      // try to attach cluster it to an existing vertex
      
      if(vtx.size() == 0) return;
      
      unsigned short iv, jv;
      short dwib, dwjb, dwie, dwje;
      bool matchUS, matchDS;
      
      for(iv = 0; iv < vtx.size(); ++iv) {
        // ignore vertices in the wrong cryostat/TPC/Plane
        if(vtx[iv].CTP != clCTP) continue;
        // ignore deleted vertices
        if(vtx[iv].NClusters == 0) continue;
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
          if(ClusterVertexChi(it, 1, iv) < 3 && ChkSignal(vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim)) {
            // good match
            tcl[it].EndVtx = iv;
            // re-fit it
            FitVtx(iv);
            if(vtxprt) mf::LogVerbatim("CC")<<"Add Begin "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<vtx[iv].ChiDOF;
            if(vtx[iv].ChiDOF < 5 && vtx[iv].WireErr < 2) return;
            tcl[it].EndVtx = -99;
            FitVtx(iv);
          } // tChi < 3
        } // matchUS
        
        if(matchDS) {
          if(ClusterVertexChi(it, 0, iv) < 3 && ChkSignal(vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim)) {
            // good match
            tcl[it].BeginVtx = iv;
            // re-fit it
            if(vtxprt) mf::LogVerbatim("CC")<<"Add End "<<tcl[it].ID<<" to vtx "<<iv<<" chi "<<vtx[iv].ChiDOF;
            if(vtx[iv].ChiDOF < 5 && vtx[iv].WireErr < 2) return;
            tcl[it].EndVtx = -99;
            FitVtx(iv);
          } // tChi < 3
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
        
  if(vtxprt) mf::LogVerbatim("CC")
      <<" ChkVertex "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" topo "<<topo<<" fvw "<<fvw<<" fvt "<<fvt;

        unsigned short vw = (unsigned short)(0.5 + fvw);
        unsigned int iht;

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

        // check vertex and clusters for proximity to existing vertices
        unsigned short nFitOK = 0;
        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if(vtx[iv].CTP != clCTP) continue;
          // make this a really loose cut since the errors on the prospective vertex aren't known yet
          if(PointVertexChi(fvw, fvt, iv) > 300) continue;
          if(vtxprt) mf::LogVerbatim("CC")<<" vtx "<<iv<<" PointVertexChi "<<PointVertexChi(fvw, fvt, iv);
          // got a match. Check the appropriate cluster end and attach
          if( (topo == 1 || topo == 2) && tcl[it1].EndVtx < 0) {
            if(ChkSignal(vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim)) {
              // try to fit
//              std::cout<<"ChkVertex attach cluster "<<tcl[it1].ID<<" to vtx "<<iv<<"\n";
              tcl[it1].EndVtx = iv;
              FitVtx(iv);
              if(vtxprt) mf::LogVerbatim("CC")<<" FitVtx "<<iv<<" WireErr "<<vtx[iv].WireErr<<" ChiDOF "<<vtx[iv].ChiDOF;
              if(vtx[iv].WireErr < 2 && vtx[iv].ChiDOF < 5) {
                ++nFitOK;
              } else {
                // bad fit
                tcl[it1].EndVtx = -99;
                FitVtx(iv);
              } // check fit
            } // ChkSignal
//            if(vtxprt)  mf::LogVerbatim("CC")<<" 12 Attach cl "<<tcl[it1].ID<<" to vtx? "<<iv;
          } else if( (topo == 3 || topo == 4) && tcl[it1].BeginVtx < 0) {
            if(ChkSignal(vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim)) {
              tcl[it1].BeginVtx = iv;
              FitVtx(iv);
              if(vtxprt) mf::LogVerbatim("CC")<<" FitVtx "<<iv<<" WireErr "<<vtx[iv].WireErr<<" ChiDOF "<<vtx[iv].ChiDOF;
              if(vtx[iv].WireErr < 2 && vtx[iv].ChiDOF < 5) {
                ++nFitOK;
              } else {
                // bad fit
                tcl[it1].BeginVtx = -99;
                FitVtx(iv);
              } // bad fit
            } // ChkSignal
//            if(vtxprt)  mf::LogVerbatim("CC")<<" 34 Attach cl "<<tcl[it1].ID<<" to vtx? "<<iv;
          } // cluster it2
          if( (topo == 1 || topo == 3) && tcl[it2].EndVtx < 0) {
            if(ChkSignal(vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim)) {
//              std::cout<<"ChkVertex attach cluster "<<tcl[it2].ID<<" to vtx "<<iv<<"\n";
              tcl[it2].EndVtx = iv;
              FitVtx(iv);
              if(vtxprt) mf::LogVerbatim("CC")<<" FitVtx "<<iv<<" WireErr "<<vtx[iv].WireErr<<" ChiDOF "<<vtx[iv].ChiDOF;
              if(vtx[iv].WireErr < 2 && vtx[iv].ChiDOF < 5) {
                ++nFitOK;
              } else {
                // bad fit
                tcl[it2].EndVtx = -99;
                FitVtx(iv);
              } // bad fit
            } // ChkSignal
//            if(vtxprt) mf::LogVerbatim("CC")<<" 13 Attach cl "<<tcl[it2].ID<<" to vtx? "<<iv;
          } else if( (topo == 2 || topo == 4) && tcl[it2].BeginVtx < 0) {
            if(ChkSignal(vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim)) {
              tcl[it2].BeginVtx = iv;
              FitVtx(iv);
              if(vtxprt) mf::LogVerbatim("CC")<<" FitVtx "<<iv<<" WireErr "<<vtx[iv].WireErr<<" ChiDOF "<<vtx[iv].ChiDOF;
              if(vtx[iv].WireErr < 2 && vtx[iv].ChiDOF < 5) {
                ++nFitOK;
              } else {
                // bad fit
                tcl[it2].BeginVtx = -99;
                FitVtx(iv);
              } // bad fit
            } // ChkSignal
//            if(vtxprt) mf::LogVerbatim("CC")<<" 24 Attach cl "<<tcl[it2].ID<<" to vtx? "<<iv;
          } // cluster it2
          if(nFitOK > 0) {
            if(vtxprt) mf::LogVerbatim("CC")<<" Attached "<<nFitOK<<" clusters to vertex "<<iv;
            return;
          }
        } // iv

        // no match to existing vertices. Ensure that there is a wire signal between
        // the vertex and the appropriate ends of the clusters
        bool SigOK = false;
        if(topo == 1 || topo == 2) {
          SigOK = ChkSignal(vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim);
        } else {
          SigOK = ChkSignal(vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim);
        }
        if(!SigOK) return;
        
        if(topo == 1 || topo == 3) {
          SigOK = ChkSignal(vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim);
        } else {
          SigOK = ChkSignal(vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim);
        }
        if(!SigOK) return;

        VtxStore newvx;
        newvx.Wire = vw;
        newvx.Time = fvt;
        newvx.Topo = topo;
        newvx.CTP = clCTP;
        vtx.push_back(newvx);
        unsigned short iv = vtx.size() - 1;
        if(topo == 1 || topo == 2) {
          if(tcl[it1].EndVtx >= 0) {
            mf::LogError("CC")<<"ChkVertex: Coding error trying to make new vtx "<<iv<<"\n";
            vtx.pop_back();
            return;
          }
//          std::cout<<"ChkVertex attach cluster "<<tcl[it1].ID<<" to newvtx "<<iv<<"\n";
          tcl[it1].EndVtx = iv;
        } else {
          if(tcl[it1].BeginVtx >= 0) {
            mf::LogError("CC")<<"ChkVertex: Coding error trying to make new vtx "<<iv<<"\n";
            vtx.pop_back();
            return;
          }
          tcl[it1].BeginVtx = iv;
        }
        if(topo == 1 || topo == 3) {
          if(tcl[it2].EndVtx >= 0) {
            mf::LogError("CC")<<"ChkVertex: Coding error trying to make new vtx "<<iv<<"\n";
            vtx.pop_back();
            return;
          }
//          std::cout<<"ChkVertex attach cluster "<<tcl[it1].ID<<" to newvtx "<<iv<<"\n";
          tcl[it2].EndVtx = iv;
        } else {
          if(tcl[it2].BeginVtx >= 0) {
            mf::LogError("CC")<<"ChkVertex: Coding error trying to make new vtx "<<iv<<"\n";
            vtx.pop_back();
            return;
          }
          tcl[it2].BeginVtx = iv;
        }
        // fit it
        FitVtx(iv);
        // reject it if the fit is bad
        if(vtx[iv].ChiDOF < 5 && vtx[iv].WireErr < 5 && vtx[iv].TimeErr < 10) {
          if(vtxprt) mf::LogVerbatim("CC")<<" New vtx "<<iv<<" in plane "<<plane<<" topo "<<topo<<" cls "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" W:T "<<(int)vw<<":"<<(int)fvt<<" NClusters "<<vtx[iv].NClusters;
        } else {
          if(vtxprt) mf::LogVerbatim("CC")<<" Bad vtx fit "<<vtx[iv].ChiDOF<<" wire err "<<vtx[iv].WireErr<<" TimeErr "<<vtx[iv].TimeErr;
          // clobber the vertex and references to it
          vtx.pop_back();
          if(tcl[it1].BeginVtx == iv) tcl[it1].BeginVtx = -99;
          if(tcl[it1].EndVtx == iv) tcl[it1].EndVtx = -99;
          if(tcl[it2].BeginVtx == iv) tcl[it2].BeginVtx = -99;
          if(tcl[it2].EndVtx == iv) tcl[it2].EndVtx = -99;
        } // bad fit

      } // ChkVertex()

/////////////////////////////////////////
    bool ClusterCrawlerAlg::ChkSignal
      (unsigned short wire1, float time1, unsigned short wire2, float time2)
    {
      // returns  true if there is a signal on the line between
      // (wire1, time1) and (wire2, time2).
      // Be sure to set time1 < time2 if you are checking for signals on a single wire

      // Gaussian amplitude in bins of size 0.15
      const float gausAmp[20] = {1, 0.99, 0.96, 0.90, 0.84, 0.75, 0.67, 0.58, 0.49, 0.40, 0.32, 0.26, 0.20, 0.15, 0.11, 0.08, 0.06, 0.04, 0.03, 0.02};
      
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
      if(wiree < fFirstWire || wiree > fLastWire) return false;
      if(wireb < fFirstWire || wireb > fLastWire) return false;
      
      short wire0 = wiree;
      // checking a single wire?
      float slope = 0;
      bool oneWire = false;
      float prTime, prTimeLo = 0, prTimeHi = 0;
      if(wireb == wiree) {
        oneWire = true;
        if(time1 < time2) {
          prTimeLo = time1;
          prTimeHi = time2;
        } else {
          prTimeLo = time2;
          prTimeHi = time1;
        }
      } else {
        slope = (timeb - timee) / (wireb - wiree);
      }
      
      int bin;
      for(short wire = wiree; wire < wireb + 1; ++wire) {
        prTime = timee + (wire - wire0) * slope;
        unsigned short index = wire - fFirstWire;
        // skip dead wires
        if(WireHitRange[index].first == -1) continue;
        // no hits on this wire
        if(WireHitRange[index].first == -2) return false;
        unsigned int firsthit = WireHitRange[index].first;
        unsigned int lasthit = WireHitRange[index].second;
        float amp = 0;
        for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
          if(oneWire) {
            if(prTimeHi > fHits[khit].EndTick()) continue;
            if(prTimeLo < fHits[khit].StartTick()) continue;
            amp = fMinAmp;
            break;
          } else {
            // skip checking if we are far away from prTime on the positive side
            if(fHits[khit].PeakTime() - prTime > 500) continue;
            bin = fabs(fHits[khit].PeakTime() - prTime) / fHits[khit].RMS();
            if(bin > 19) continue;
	    if(bin < 0) continue;
            // add amplitude from all hits
            amp += fHits[khit].PeakAmplitude() * gausAmp[bin];
          }
        } // khit
        if(amp < fMinAmp) return false;
      } // wire
      return true;
    
    } // ChkSignal()

/////////////////////////////////////////
    bool ClusterCrawlerAlg::SplitCluster
      (unsigned short icl, unsigned short pos, unsigned short ivx)
    {
      // split cluster icl into two clusters

      if(tcl[icl].ID < 0) return false;

      if(pos < 3 || pos > tcl[icl].tclhits.size() - 3) return false;
      
      // declare icl obsolete
      MakeClusterObsolete(icl);

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
      hitNear.clear();
      chgNear.clear();
      for(unsigned short ii = 0; ii < pos; ++ii) {
        unsigned int iht = tcl[icl].tclhits[ii];
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
      if(!TmpStore()) {
        RestoreObsoleteCluster(icl);
        return false;
      }
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
      hitNear.clear();
      chgNear.clear();
      bool didFit = false;
      for(unsigned short ii = pos; ii < tcl[icl].tclhits.size(); ++ii) {
        unsigned int iht = tcl[icl].tclhits[ii];
        if(inClus[iht] != 0) {
          RestoreObsoleteCluster(icl);
          return false;
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
      if(!TmpStore()) {
        // clobber the previously store cluster
        MakeClusterObsolete(tcl.size() - 1);
        RestoreObsoleteCluster(icl);
        return false;
      }
      // associate the End with the supplied vertex
      iclnew = tcl.size() - 1;
      tcl[iclnew].BeginVtx = ivx;
      tcl[iclnew].EndVtx = tcl[icl].EndVtx;
			return true;
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
      bool bothLong, NoVtx;
      
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
    mf::LogVerbatim("CC")
      <<"Chk1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(ew1 - bw2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          if(NoVtx && bw2 < ew1 && (ew1 - bw2)  < skipcut && dth < angcut) {
            chgrat = 2 * fabs(bc2 - ec1) / (bc2 + ec1);
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.05) chgrat = 0.;
            // project bw2,bt2 to ew1
            dtim = fabs(bt2 + (ew1-bw2)*bs2 - et1);
  if(prt) {
    mf::LogVerbatim("CC")
    <<" dtim "<<dtim<<" timecut "<<(int)timecut
    <<" ec1 "<<ec1<<" bc2 "<<bc2
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut<<" es1 "<<es1<<" slpcut "<<fLAClusSlopeCut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              // ensure there is a signal between cluster ends
              if(fuBCode) {
                DoMerge(it2, it1, 10);
                tclsize = tcl.size();
                break;
              } else if(ChkSignal(ew1,et1,bw2,bt2)) {
                DoMerge(it2, it1, 10);
                tclsize = tcl.size();
                break;
              }
            } // chgrat < chgcut ...
          } // US cluster 2 merge with DS cluster 1?
          
          // look for US and DS broken clusters w similar angle
          // US cluster 1 merge with DS cluster 2?
          dth = fabs(bth1 - eth2);
  if(prt && bw1 < ew2 && (ew2 - bw1)  < skipcut) {
    mf::LogVerbatim("CC")
      <<"Chk2 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(bw1 - ew2)
      <<" skipcut "<<skipcut
      <<" dth "<<dth<<" angcut "<<angcut;
  }
          // require no vertex between
          NoVtx = (tcl[it2].EndVtx < 0) && (tcl[it1].BeginVtx < 0);
          if(NoVtx && bw1 < ew2 && (ew2 - bw1)  < skipcut && dth < angcut ) {
            chgrat = 2 * fabs((bc1 - ec2) / (bc1 + ec2));
            // ignore the charge cut for long tracks with small dth
            if(bothLong && dth < 0.05) chgrat = 0.;
            // project sw1,st1 to ew2
            dtim = std::abs(bt1 + (ew2-bw1)*bs1 - et2);
  if(prt) {
    mf::LogVerbatim("CC")
    <<" dtim "<<dtim<<" err "<<dtim<<" timecut "<<(int)timecut
    <<" chgrat "<<chgrat<<" chgcut "<<chgcut;
  }
            if(chgrat < chgcut && dtim < timecut) {
              if(fuBCode ) {
                DoMerge(it1, it2, 10);
                tclsize = tcl.size();
                break;
              } else if(ChkSignal(bw1,bt1,ew2,et2)) {
                DoMerge(it1, it2, 10);
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
    mf::LogVerbatim("CC")
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
    mf::LogVerbatim("CC")
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
    
  if(prt) mf::LogVerbatim("CC")<<"ChkMerge12 "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    
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
      unsigned int hit = cl1.tclhits[iht];
      unsigned short wire = fHits[hit].WireID().Wire;
      if(wire - ew1 < 0 || wire - ew1 > (short)cl1hits.size()) {
        mf::LogError("CC")<<"ChkMerge12 bad wire "<<(wire-ew1);
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
  if(prt) mf::LogVerbatim("CC")<<"chk next US wire "<<wiron1<<" missed "<<nmiss;
    if(wiron1 == 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    
    // compare the wires with hits on cluster 2 with the gap in cluster 1
    // the number of hit wires that fit in the gap
    unsigned short nfit = 0;
    for(unsigned short iht = 0; iht < tcl[it2].tclhits.size(); ++iht) {
      unsigned int hiton2 = tcl[it2].tclhits[iht];
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
      mf::LogError("CC")<<"ChkMerge12 bad hiton1 "<<hiton1;
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
  if(prt) mf::LogVerbatim("CC")<<"US dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut. fAveChg was calculated in FitClusterMid
    float chgrat = 2 * std::abs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
  if(prt) mf::LogVerbatim("CC")<<"US chgrat "<<chgrat<<" cut "<<fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    bool SigOK;
    if(fuBCode) {
      SigOK = true;
    } else {
      SigOK = ChkSignal(wiron1, timon1, ew2, et2);
    }
  if(prt) mf::LogVerbatim("CC")<<"US SigOK? "<<SigOK;
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
      mf::LogError("CC")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    timon1 = fHits[hiton1].PeakTime();
    dtim = std::abs(bt2 - (wiron1 - bw2) *tcl[it2].BeginSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    FitClusterMid(it1, hiton1, -3);
    if(clChisq > 20.) return;
    // check for angle consistency
    dth = std::abs(std::atan(fScaleF * clpar[1]) - std::atan(fScaleF * tcl[it2].BeginSlp));
  if(prt) mf::LogVerbatim("CC")
      <<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass];
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * std::abs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(prt) mf::LogVerbatim("CC")
      <<"DS chgrat "<<chgrat<<" cut "<<fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    if(fuBCode) {
      SigOK = true;
    } else {
      SigOK = ChkSignal(wiron1, timon1, bw2, bt2);
    }
  if(prt) mf::LogVerbatim("CC")<<"DS SigOK? "<<SigOK;
    if( !SigOK ) return;

  if(prt) mf::LogVerbatim("CC")<<"Merge em";
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
    MakeClusterObsolete(it1);
    MakeClusterObsolete(it2);
    
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
    std::vector<int> wirehit;
    for(unsigned short wire = lowire; wire < hiwire + 2; ++wire) {
      wirehit.push_back(-1);
    }
    // put in the hit IDs for cluster 2
    for(unsigned short iht = 0; iht < cl2.tclhits.size(); ++iht) {
      unsigned int hit = cl2.tclhits[iht];
      // un-assign the hit
      inClus[hit] = 0;
      unsigned short wire = fHits[hit].WireID().Wire;
      unsigned short index = wire - lowire;
      wirehit[index] = hit;
/*
  if(myprt) mf::LogVerbatim("CC")
    <<"Cl2 hit "<<wire<<":"<<(int)fHits[hit].PeakTime()
    <<" wire index "<<index<<" hit index "<<hit;
*/
     } // iht
    // now cluster 1
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned int hit = cl1.tclhits[iht];
      unsigned short wire = fHits[hit].WireID().Wire;
      unsigned short index = wire - lowire;
      inClus[hit] = 0;
      // TODO: Should merge hits instead of clobbering one of them?
      wirehit[index] = hit;
/*
  if(myprt) mf::LogVerbatim("CC")
    <<"Cl1 hit "<<wire<<":"<<(int)fHits[hit].PeakTime()
    <<" wire index "<<index<<" hit index "<<hit;
*/
    } // iht
    // make the new cluster
    fcl2hits.clear();
    for(std::vector<int>::reverse_iterator rit = wirehit.rbegin();
        rit != wirehit.rend(); ++rit) {
      if(*rit < 0) continue;
      unsigned int hit = *rit;
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
      clBeginChgNear = cl1.BeginChgNear;
      begVtx = cl1.BeginVtx;
      del1Vtx = cl1.EndVtx;
      // and cluster 2 End info
      clEndSlp = cl2.EndSlp;
      clEndSlpErr = cl2.EndSlpErr;
      clEndAng = cl2.EndAng;
      clEndWir = cl2.EndWir;
      clEndTim = cl2.EndTim;
      clEndChg = cl2.EndChg;
      clEndChgNear = cl2.EndChgNear;
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
      clBeginChgNear = cl2.BeginChgNear;
      begVtx = cl2.BeginVtx;
      del2Vtx = cl2.EndVtx;
      // and cluster 1 End info
      clEndSlp = cl1.EndSlp;
      clEndSlpErr = cl1.EndSlpErr;
      clEndWir = cl1.EndWir;
      clEndTim = cl1.EndTim;
      clEndChg = cl1.EndChg;
      clEndChgNear = cl1.EndChgNear;
      endVtx = cl1.EndVtx;
      del1Vtx = cl1.BeginVtx;
      clStopCode = cl1.StopCode;
    }

    // append it to the tcl vector
    clCTP = cl1.CTP;
    if(!TmpStore()) return;
    unsigned short itnew = tcl.size()-1;
  if(prt) mf::LogVerbatim("CC")<<"DoMerge "<<cl1.ID<<" "<<cl2.ID<<" -> "<<tcl[itnew].ID;
    // stuff the processor code with the current pass
    tcl[itnew].ProcCode = inProcCode + pass;
    // transfer the vertex info
    // delete a vertex between these two?
    if(del1Vtx >= 0 && del1Vtx == del2Vtx) vtx[del1Vtx].NClusters = 0;
    // preserve the vertex assignments
    tcl[itnew].BeginVtx = begVtx;
//    std::cout<<"DoMerge attach cluster "<<tcl[itnew].ID<<" to vtx "<<endVtx<<"\n";
    tcl[itnew].EndVtx = endVtx;
  } // DoMerge

/////////////////////////////////////////
  void ClusterCrawlerAlg::PrintClusters()
  {

    // prints clusters to the screen for code development
    mf::LogVerbatim myprt("CC");

    if(vtx3.size() > 0) {
      // print out 3D vertices
      myprt<<"****** 3D vertices ******************************************__2DVtx_Indx__*******\n";
      myprt<<"Vtx  Cstat  TPC Proc     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for(unsigned short iv = 0; iv < vtx3.size(); ++iv) {
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(7)<<vtx3[iv].CStat;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].TPC;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].ProcCode;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].X;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].Y;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].Z;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].XErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].YErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].ZErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[0];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[1];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[2];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Wire;
        if(vtx3[iv].Wire < 0) {
          myprt<<"    Matched in 3 planes";
        } else {
          myprt<<"    Incomplete";
        }
        myprt<<"\n";
      }
    } // vtx3.size
    
    if(vtx.size() > 0) {
      // print out 2D vertices
      myprt<<"************ 2D vertices ************\n";
      myprt<<"Vtx   CTP    wire     error   tick     error  ChiDOF  NCl  topo  cluster IDs\n";
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(fDebugPlane < 3 && fDebugPlane != (int)vtx[iv].CTP) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<vtx[iv].CTP;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Wire<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].WireErr;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Time<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].TimeErr;
        myprt<<std::right<<std::setw(8)<<vtx[iv].ChiDOF;
        myprt<<std::right<<std::setw(5)<<vtx[iv].NClusters;
        myprt<<std::right<<std::setw(6)<<vtx[iv].Topo;
        myprt<<"    ";
        // display the cluster IDs
        for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
          if(fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
          if(tcl[ii].ID < 0) continue;
          if(tcl[ii].BeginVtx == (short)iv) myprt<<std::right<<std::setw(4)<<tcl[ii].ID;
          if(tcl[ii].EndVtx == (short)iv) myprt<<std::right<<std::setw(4)<<tcl[ii].ID;
        }
        myprt<<"\n";
      } // iv
    } // vtx.size
    
    float aveRMS = 0;
    myprt<<"*************************************** Clusters *****************************************************\n";
    myprt<<"  ID CTP nht Stop  Proc  beg_W:T    bAng   bSlp bChg  bCN   end_W:T    eAng   eSlp eChg  eCN bVx  eVx aveRMS\n";
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      // print clusters in all planes (fDebugPlane = 3) or in a selected plane
      if(fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
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
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].BeginSlp;
      myprt<<std::right<<std::setw(5)<<(int)tcl[ii].BeginChg;
      myprt<<std::right<<std::setw(5)<<std::fixed<<std::setprecision(1)<<tcl[ii].BeginChgNear;
      iTime = tcl[ii].EndTim;
      myprt<<std::right<<std::setw(6)<<tcl[ii].EndWir<<":"<<iTime;
      if(iTime < 10) {
        myprt<<"   ";
      } else if(iTime < 100) {
        myprt<<"  ";
      } else if(iTime < 1000) myprt<<" ";
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].EndAng;
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].EndSlp;
      myprt<<std::right<<std::setw(5)<<(int)tcl[ii].EndChg;
      myprt<<std::right<<std::setw(5)<<std::fixed<<std::setprecision(1)<<tcl[ii].EndChgNear;
      myprt<<std::right<<std::setw(5)<<tcl[ii].BeginVtx;
      myprt<<std::right<<std::setw(5)<<tcl[ii].EndVtx;
      aveRMS = 0;
      unsigned int iht = 0;
      for(unsigned short jj = 0; jj < tcl[ii].tclhits.size(); ++jj) {
        iht = tcl[ii].tclhits[jj];
        aveRMS += fHits[iht].RMS();
      }
      aveRMS /= (float)tcl[ii].tclhits.size();
      myprt<<std::right<<std::setw(5)<<std::fixed<<std::setprecision(1)<<aveRMS;
      myprt<<"\n";
      // TEMP
//      if(tcl[ii].BeginVtx >= 0) myprt<<"Begin vtx Chi "<<ClusterVertexChi(ii, 0, tcl[ii].BeginVtx)<<"\n";
//      if(tcl[ii].EndVtx >= 0) myprt<<"End vtx Chi "<<ClusterVertexChi(ii, 1, tcl[ii].EndVtx)<<"\n";
    } // ii
    
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
      clBeginChgNear = tcl[it1].BeginChgNear;
      clEndSlp = tcl[it1].EndSlp;
      clEndSlpErr = tcl[it1].EndSlpErr;
      clEndAng = tcl[it1].EndAng;
      clEndWir = tcl[it1].EndWir;
      clEndTim = tcl[it1].EndTim;
      clEndChg = tcl[it1].EndChg;
      clEndChgNear = tcl[it1].EndChgNear;
      clStopCode = tcl[it1].StopCode;
      clProcCode = tcl[it1].ProcCode;
      clCTP = tcl[it1].CTP;
      fcl2hits = tcl[it1].tclhits;
    }


/////////////////////////////////////////
  bool ClusterCrawlerAlg::TmpStore()
  {

    if(fcl2hits.size() < 3) return false;
    
    if(fcl2hits.size() == UINT_MAX) return false;
    
    if(NClusters == USHRT_MAX) return false;
    
    ++NClusters;
    
    unsigned int hit0 = fcl2hits[0];
    unsigned int tCST = fHits[hit0].WireID().Cryostat;
    unsigned int tTPC = fHits[hit0].WireID().TPC;
    unsigned int tPlane = fHits[hit0].WireID().Plane;
    unsigned short lastWire = 0;
    
    for(unsigned short it = 0; it < fcl2hits.size(); ++it) {
      unsigned int hit = fcl2hits[it];
      if(inClus[hit] != 0) {
        mf::LogError("CC")<<"TmpStore: Trying to use obsolete/used hit "<<hit<<" inClus "<<inClus[hit]
          <<" on wire "<<fHits[hit].WireID().Wire<<" on cluster "<<NClusters
          <<" in plane "<<plane<<" ProcCode "<<clProcCode;
        --NClusters;
        return false;
      }
      // check for WireID() consistency
      if(fHits[hit].WireID().Cryostat != tCST || fHits[hit].WireID().TPC != tTPC || fHits[hit].WireID().Plane != tPlane) {
        mf::LogError("CC")<<"TmpStore: CTP mis-match "<<hit<<" WireID().TPC "<<fHits[hit].WireID().TPC
        <<" WireID().Plane "<<fHits[hit].WireID().Plane<<" tCST "<<tCST<<" tTPC "<<tTPC<<" tPlane "<<tPlane
        <<" on cluster "<<NClusters<<" ProcCode "<<clProcCode;
        --NClusters;
        return false;
      }
      // check for decreasing wire number
      if(it > 0 && fHits[hit].WireID().Wire >= lastWire) {
        mf::LogError("CC")<<"TmpStore: Hits not in correct wire order "<<NClusters;
        --NClusters;
        return false;
      }
      lastWire = fHits[hit].WireID().Wire;
      inClus[hit] = NClusters;
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
    
    std::vector<unsigned int>::const_iterator ibg = fcl2hits.begin();
    unsigned short hitb = *ibg;
    std::vector<unsigned int>::const_iterator iend = fcl2hits.end() - 1;
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
    clstr.BeginChgNear = clBeginChgNear;
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndAng      = std::atan(fScaleF * clEndSlp);
    clstr.EndWir      = fHits[hite].WireID().Wire;
    clstr.EndTim      = fHits[hite].PeakTime();
    clstr.EndChg      = clEndChg;
    clstr.EndChgNear  = clEndChgNear;
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.BeginVtx    = -99;
    clstr.EndVtx      = -99;
    clstr.CTP         = EncodeCTP(tCST, tTPC, tPlane);
    clstr.tclhits     = fcl2hits;
    tcl.push_back(clstr);
    return true;
  } // TmpStore()

/////////////////////////////////////////
  void ClusterCrawlerAlg::LACrawlUS() {
    // Crawl a large angle cluster upstream. Similar to CrawlUS but require
    // that a hit be added on each wire


    unsigned int dhit = fcl2hits[0];
    short dwir = fHits[dhit].WireID().Wire;
    prt = false;
  if(fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
    prt = std::abs(fHits[dhit].PeakTime() - fDebugHit) < 20;

  if(prt) {
    mf::LogVerbatim myprt("CC");
    myprt<<"******************* LACrawlUS PASS "<<pass<<" Hits ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned int iht = fcl2hits[fcl2hits.size() - 1 - ii];
      myprt<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" ";
    }
  }

    // Merge all the hits if this is a very large angle cluster
    if(std::abs(clpar[1]) > 5) {
      unsigned short ii, nmult = 0;
      for(ii = 0; ii < fcl2hits.size(); ++ii) 
        if(fHits[fcl2hits[ii]].Multiplicity() > 1) ++nmult;
      if(nmult > 0.5 * fcl2hits.size()) {
        bool didMerge;
        for(ii = 0; ii < fcl2hits.size(); ++ii) {
          MergeHits(fcl2hits[ii], didMerge);
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
    unsigned int it = fcl2hits.size() - 1;
    unsigned int lasthit = fcl2hits[it];
    unsigned short lastwire = fHits[lasthit].WireID().Wire;
    bool ChkCharge = false;
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) mf::LogVerbatim("CC")<<"LACrawlUS: next wire "<<nextwire;
      // stop crawling if there is a nearby vertex
      if(CrawlVtxChk(nextwire)) {
  if(prt) mf::LogVerbatim("CC")<<"LACrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // AddLAHit will merge the hit on nextwire if necessary
      AddLAHit(nextwire, ChkCharge, HitOK, SigOK);
  if(prt) mf::LogVerbatim("CC")<<"LACrawlUS: HitOK "<<HitOK<<" SigOK "<<SigOK<<" nmissed "<<nmissed;
      if(fuBCode) {
        SigOK = true;
      } else {
        if(!SigOK) break;
      }
      if(HitOK) nmissed = 0;
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
        bool didMerge;
        for(unsigned short kk = 0; kk< fcl2hits.size()-1; ++kk) {
          unsigned int hit = fcl2hits[kk];
          MergeHits(hit, didMerge);
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
        mf::LogError("CC")
          <<"LACrawlUS: chifits size error "<<chifits.size()<<" "<<fcl2hits.size();
        return;
      }
  if(prt) {
    mf::LogVerbatim("CC")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
      if( chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
          chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
        // find the kink angle (crudely) from the 0th and 2nd hit
        unsigned int ih0 = fcl2hits.size() - 1;
        unsigned int hit0 = fcl2hits[ih0];
        unsigned int ih2 = ih0 - 2;
        unsigned int hit2 = fcl2hits[ih2];
        float dt02 = fHits[hit2].PeakTime() - fHits[hit0].PeakTime();
        float dw02 = fHits[hit2].WireID().Wire - fHits[hit0].WireID().Wire;
        float th02 = std::atan( fScaleF * dt02 / dw02);
        // and the 3rd and 5th hit
        unsigned int ih3 = ih2 - 1;
        unsigned int hit3 = fcl2hits[ih3];
        unsigned int ih5 = ih3 - 2;
        unsigned int hit5 = fcl2hits[ih5];
        float dt35 = fHits[hit5].PeakTime() - fHits[hit3].PeakTime();
        float dw35 = fHits[hit5].WireID().Wire - fHits[hit3].WireID().Wire;
        float th35 = std::atan(fScaleF * dt35 / dw35);
        float dth = std::abs(th02 - th35);
  if(prt) mf::LogVerbatim("CC")<<" Kink angle "<<std::setprecision(3)<<dth<<" cut "<<fKinkAngCut[pass];
        if(dth > fKinkAngCut[pass]) {
          // hit a kink. Lop of the first 3 hits, refit and keep crawling?
          for(short jj = 0; jj < 3; ++jj) {
            fcl2hits.pop_back();
            chifits.pop_back();
            hitNear.pop_back();
            chgNear.pop_back();
          }
          FitCluster();
          // See if this is a second kink and it is close to the first
          // kink (which had hits removed).
          if(kinkOnWire > 0) {
            if(kinkOnWire - nextwire < 4) {
  if(prt) mf::LogVerbatim("CC")
    <<"Hit a second kink. kinkOnWire = "<<kinkOnWire<<" Stopping";
              // set the kink stop code
              clStopCode = 3;
              break;
            }
          }
          kinkOnWire = nextwire;
  if(prt) mf::LogVerbatim("CC")<<"Removed kink hits";
        } // kinkang check
      } // chifits test
    } // nextwire

    CheckClusterHitFrac(prt);

    clProcCode += 300;
  if(prt) mf::LogVerbatim("CC")
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
    mf::LogVerbatim myprt("CC");
    myprt<<"******************* CrawlUS PASS "<<pass<<" Hits: ";
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      unsigned int iht = fcl2hits[fcl2hits.size() - 1 - ii];
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
    unsigned int it = fcl2hits.size() - 1;
    unsigned int lasthit = fcl2hits[it];
    if(lasthit > fHits.size() - 1) {
      mf::LogError("CC")<<"CrawlUS bad lasthit "<<lasthit;
    }
    unsigned short lastwire = fHits[lasthit].WireID().Wire;
    if(prt) mf::LogVerbatim("CC")<<"CrawlUS: last wire "<<lastwire<<" hit "<<lasthit;
    
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
      if(prt) mf::LogVerbatim("CC")<<"CrawlUS: next wire "<<nextwire;
      // stop crawling if there is a nearby vertex
      if(CrawlVtxChk(nextwire)) {
        if(prt) mf::LogVerbatim("CC")<<"CrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // add hits and check for PH and width consistency
      AddHit(nextwire, HitOK, SigOK);
      if(prt) mf::LogVerbatim("CC")<<"CrawlUS: HitOK "<<HitOK<<" SigOK "<<SigOK<<" nmissed "<<nmissed;
      if(!HitOK) {
        // no hit on this wire. Was there a signal or dead wire?
        if(SigOK) {
        // no hit on the wire but there is a signal
          ++nmissed;
          // stop if too many missed wires
          if(prt) mf::LogVerbatim("CC")<<" nmissed "<<nmissed<<" cut "<<fMaxWirSkip[pass];
          if(nmissed > fMaxWirSkip[pass]) {
            clStopCode = 1;
            break;
          }
/*
          // see if we are in the PostSkip phase and missed more than 1 wire
// is this an error?
//          if(PostSkip && nmissed > 1) {
          if(PostSkip && nmissed > fMinWirAfterSkip[pass]) {
            // cluster is really short
            if((int)(fcl2hits.size() - nHitAfterSkip) < 4) {
              fcl2hits.clear();
              return;
            }
  if(prt) mf::LogVerbatim("CC")<<" PostSkip && nmissed = "<<nmissed;
            clStopCode = 2;
            for(short jj = 0; jj < nHitAfterSkip; ++jj) {
              fcl2hits.pop_back();
              chifits.pop_back();
              hitNear.pop_back();
              chgNear.pop_back();
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
*/
        } // SigOK
        else {
          // SigOK is false
          clStopCode = 0;
          if(prt) mf::LogVerbatim("CC")<<"No hit or signal on wire "<<nextwire;
          break;
        } // else SigOK false
      } // !HitOK
      else {
        if(clChisq > 99.) {
          if(fcl2hits.size() < 3) return;
          // a fit error occurred. Lop off the leading hit and refit
  if(prt) mf::LogVerbatim("CC")<<"Fit failed ";
          fcl2hits.pop_back();
          chifits.pop_back();
          hitNear.pop_back();
          chgNear.pop_back();
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
            mf::LogError("CC")
              <<"CrawlUS: chifits size error "<<chifits.size()<<" "<<fcl2hits.size();
            return;
          }
          unsigned short chsiz = chifits.size()-1;
  if(prt) {
    mf::LogVerbatim("CC")
    <<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" "
    <<chifits[chsiz-2]<<" "<<chifits[chsiz-3];
  }
          if( chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
              chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
            if(fcl2hits.size() != chifits.size()) {
              mf::LogError("CC")
              <<"bad kink check size "<<chifits.size()<<" "<<fcl2hits.size()
              <<" plane "<<plane<<" cluster "<<dwir<<":"<<dhit;
              continue;
            }
            // find the kink angle (crudely) from the 0th and 2nd hit
            unsigned short ih0 = fcl2hits.size() - 1;
            unsigned int hit0 = fcl2hits[ih0];
            unsigned short ih2 = ih0 - 2;
            unsigned int hit2 = fcl2hits[ih2];
            float dt02 = fHits[hit2].PeakTime() - fHits[hit0].PeakTime();
            float dw02 = fHits[hit2].WireID().Wire - fHits[hit0].WireID().Wire;
            float th02 = std::atan( fScaleF * dt02 / dw02);
            // and the 3rd and 5th hit
            unsigned short ih3 = ih2 - 1;
            unsigned int hit3 = fcl2hits[ih3];
            unsigned short ih5 = ih3 - 2;
            unsigned int hit5 = fcl2hits[ih5];
            float dt35 = fHits[hit5].PeakTime() - fHits[hit3].PeakTime();
            float dw35 = fHits[hit5].WireID().Wire - fHits[hit3].WireID().Wire;
            float th35 = std::atan(fScaleF * dt35 / dw35);
            float dth = std::abs(th02 - th35);
  if(prt) mf::LogVerbatim("CC")<<" Kink angle "<<std::setprecision(3)<<dth<<" cut "<<fKinkAngCut[pass];
            // cut on the allowed kink angle
            if(dth > fKinkAngCut[pass]) {
  if(prt) mf::LogVerbatim("CC")<<" Stopped tracking ";
              // kill the last 3 hits, refit and return
              for(short jj = 0; jj < 3; ++jj) {
                fcl2hits.pop_back();
                chifits.pop_back();
                hitNear.pop_back();
                chgNear.pop_back();
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
  if(prt) mf::LogVerbatim("CC")<<" Set clBeginChg "<<clBeginChg;
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
  if(prt) mf::LogVerbatim("CC")<<"chirat "<<chirat
    <<" last hit "<<fcl2hits[fcl2hits.size()-1];
            if(chirat < 1.2) break;
            fcl2hits.pop_back();
            chifits.pop_back();
            hitNear.pop_back();
            chgNear.pop_back();
            if(fcl2hits.size() < 6) break;
            if(chifits.size() < 6) break;
          } // nlop
          if(fcl2hits.size() < 6) {
            clStopCode = 4;
  if(prt) mf::LogVerbatim("CC")<<"Bad fit chisq - short cluster. Break";
            break;
          }
          FitCluster();
          FitClusterChg();
  if(prt) mf::LogVerbatim("CC")<<"Bad fit chisq - removed hits. Continue...";
        } // clChisq > fChiCut[pass]
      } // !HitOK check
    } // nextwire

    // count the number of hits on adjacent wires at the leading edge and
    // ensure that the count is consistent with fMinWirAfterSkip
    bool reFit = false;
    if((unsigned short)fcl2hits.size() > fMinWirAfterSkip[pass] + 3) {
      unsigned short ih0 = fcl2hits.size() - 1;
      unsigned int hit0 = fcl2hits[ih0];
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
  if(prt) mf::LogVerbatim("CC")<<"Lop off "<<nAdjHit<<" hits ";
        for(unsigned short ii = 0; ii < nAdjHit + 1; ++ii) {
          fcl2hits.pop_back();
          chifits.pop_back();
          hitNear.pop_back();
          chgNear.pop_back();
          reFit = true;
        }
      }
    } // fcl2hits.size() > fMinWirAfterSkip[pass] + 3
    if(reFit) {
      FitCluster();
      FitClusterChg();
    }
    CheckClusterHitFrac(prt);
    if(prt) mf::LogVerbatim("CC")<<"CheckClusterHitFrac "<<fcl2hits.size()<<" min length for this pass "<<fMinHits[pass];

    prt = false;
  } // CrawlUS()

/////////////////////////////////////////
  void ClusterCrawlerAlg::CheckClusterHitFrac(bool prt)
  {


    // Find the fraction of the wires on the cluster that have
    // hits.
    unsigned int iht = fcl2hits[fcl2hits.size() - 1];
    clEndWir = fHits[iht].WireID().Wire;
    float hitFrac = (float)fcl2hits.size() / (float)(clBeginWir - clEndWir + 1);

    if(hitFrac < fMinHitFrac) {
      fcl2hits.clear();
      if(prt) mf::LogVerbatim("CC")
        <<"CheckClusterHitFrac: Poor hit fraction "<<hitFrac;
      return;
    } // hitFrac
    
    // lop off the last hit if it is part of a hit multiplet
    if(fHits[iht].Multiplicity() > 1) {
      fcl2hits.resize(fcl2hits.size() - 1);
    }

    // check for short track ghosts
    if(fcl2hits.size() < 5) {
      unsigned short nsing = 0;
      for(unsigned short iht = 0; iht < fcl2hits.size(); ++iht) if(fHits[fcl2hits[iht]].Multiplicity() == 1) ++nsing;
      hitFrac = (float)nsing / (float)fcl2hits.size();
      if(hitFrac < fMinHitFrac) {
        fcl2hits.clear();
        if(prt) mf::LogVerbatim("CC")<<"CheckClusterHitFrac: Poor short track hit fraction "<<hitFrac;
        return;
      } // hitFrac
    } // short ghost track check
    
    // analyze the pattern of nearby charge
    // will need the cluster charge so calculate it here if it isn't defined yet
    if(clBeginChg <= 0) {
      unsigned int iht, nht = 0;
      for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
        iht = fcl2hits[ii];
        clBeginChg += fHits[iht].Integral();
        ++nht;
        if(nht == 8) break;
      }
      clBeginChg /= (float)nht;
    } // clBeginChg == 0
    // handle short vs long clusters
    unsigned short cnt = chgNear.size()/2;
    // get the average charge from <= 30 hits at each end
    if(chgNear.size() > 60) cnt = 30;
    clBeginChgNear = 0;
    clEndChgNear = 0;
    for(unsigned short ids = 0; ids < cnt; ++ids) {
      clBeginChgNear += chgNear[ids];
      clEndChgNear += chgNear[chgNear.size() - 1 - ids];
    }
    clBeginChgNear /= (float)cnt;
    clEndChgNear /= (float)cnt;
    
  } // CheckClusterHitFrac()

/////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterMid(unsigned short it1, unsigned int ihtin, short nhit)
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
        unsigned int ihit = cls.tclhits[it];
        if(ihit > fHits.size()-1) {
          mf::LogError("CC")<<"FitClusterMid bad ihit "<<ihit;
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
          mf::LogVerbatim("CC")<<"FitClusterMid bad ihit "<<ihit;
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
    clpar[2] = wire0;
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
    
    if(pass > fNumPass - 1) {
      mf::LogError("CC")<<"FitCluster called on invalid pass "<<pass;
      return;
    }
    
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
    for(std::vector<unsigned int>::reverse_iterator it = fcl2hits.rbegin();
       it != fcl2hits.rend(); ++it) {
      unsigned int ihit = *it;
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
    mf::LogVerbatim myprt("CC");
    myprt<<"FitCluster W:T ";
    unsigned short cnt = 0;
    for(std::vector<unsigned int>::reverse_iterator it = fcl2hits.rbegin();
       it != fcl2hits.rend(); ++it) {
      unsigned int ihit = *it;
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
    clpar[2] = wire0;
    clparerr[0] = intcpterr;
    clparerr[1] = slopeerr;

  if(prt) mf::LogVerbatim("CC")
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
      mf::LogError("CC")<<"FitClusterChg bad pass "<<pass;
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
      for(unsigned short ii = fcl2hits.size() - 1; ii > imlast; --ii) {
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
  if(prt) mf::LogVerbatim("CC")<<"FitClusterChg wire "<<wire0
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
    // check for signal on a known good wire
    if(fuBCode) {
      // List of uB bad wires not defined yet so assume all
      // wires have a signal
      SigOK = true;
    } else {
      // skip bad wire, but assume the track was there
      if(WireHitRange[index].first == -1) {
        SigOK = true;
        return;
      }
      // return SigOK false if no hit on a good wire
      if(WireHitRange[index].first == -2) return;
    }

    unsigned int firsthit = WireHitRange[index].first;
    unsigned int lasthit = WireHitRange[index].second;
    if(fuBCode) {
      if(prt && firsthit > fHits.size()) mf::LogError("CC")<<"AddLAHit: Bad firsthit "<<firsthit<<" size "<<fHits.size();
      if(firsthit > fHits.size()) return;
    } // fuBCode
    
    // max allowable time difference between projected cluster and a hit
    float timeDiff = 40 * AngleFactor(clpar[1]);
    float dtime;
    
    // the last hit added to the cluster
    unsigned int lastClHit = fcl2hits[fcl2hits.size()-1];
    unsigned short wire0 = fHits[lastClHit].WireID().Wire;
    
    if(fuBCode && (wire0 - kwire) > fMaxWirSkip[pass]) return;

    // the projected time of the cluster on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    float prtimeErr = (kwire - wire0) * clparerr[1];
    float chgWinLo = prtime - fChgNearWindow;
    float chgWinHi = prtime + fChgNearWindow;
    float chgrat;
    float cnear = 0;
    if(prt) mf::LogVerbatim("CC")<<"AddLAHit: wire "<<kwire<<" prtime "<<(int)prtime<<" prtimeErr "<<prtimeErr;
    unsigned int imbest = 0;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(inClus[khit] < 0) continue;
      chgrat = fHits[khit].Integral() / fHits[lastClHit].Integral();
  if(prt) mf::LogVerbatim("CC")
    <<" Chk W:T "<<kwire<<":"<<(short)fHits[khit].PeakTime()
    <<" Charge "<<(short)fHits[khit].Integral()
    <<" chgrat "<<std::setw(8)<<std::setprecision(2)<<chgrat
    <<" inClus "<<inClus[khit]
    <<" mult "<<fHits[khit].Multiplicity()
    <<" RMS "<<std::setprecision(2)<<fHits[khit].RMS()
    <<" Chi2 "<<std::fixed<<std::setw(10)<<std::setprecision(2)<<fHits[khit].GoodnessOfFit()
    <<" LoT "<<(int)fHits[khit].StartTick()
    <<" HiT "<<(int)fHits[khit].EndTick()
    <<" index " << khit;
      // count charge in the window
      if(fHits[khit].PeakTime() > chgWinLo && fHits[khit].PeakTime() < chgWinHi) cnear += fHits[khit].Integral();
      // projected time outside the Signal time window?
      if(prtime < fHits[khit].StartTick() - 20) continue;
      if(prtime > fHits[khit].EndTick() + 20) continue;
      SigOK = true;
      // hit used?
      if(inClus[khit] > 0) continue;
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
    
    if(prt && !HitOK) mf::LogVerbatim("CC")<<" no hit found ";

    if(!HitOK) return;

    if(prt) mf::LogVerbatim("CC")<<" Pick hit time "<<(int)fHits[imbest].PeakTime()<<" hit index "<<imbest;
    
    // merge hits in a multiplet?
    short hnear = 0;
    if(fHits[imbest].Multiplicity() > 1) {
      bool doMerge = true;
      if(!fuBCode) {
        // Standard code
        // don't merge if we are close to a vertex
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].CTP != clCTP) continue;
          if(prt) mf::LogVerbatim("CC")
            <<" close vtx chk W:T "<<vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time;
          if(std::abs(kwire - vtx[ivx].Wire) < 5 &&
             std::abs(int(fHits[imbest].PeakTime() - vtx[ivx].Time)) < 20 ) {
            if(prt) mf::LogVerbatim("CC")
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
            if(inClus[jht] < 0) continue;
            //          if(!isHitPresent(jht)) continue;
            // count used hits
            //          if(isHitFree(jht)) multipletChg += fHits[jht].Integral();
            if(inClus[jht] == 0) multipletChg += fHits[jht].Integral();
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
            if(!doMerge) mf::LogVerbatim("CC")
              <<" Hits are well separated. Don't merge them ";
          }
          if(doMerge && nused == 0) {
            // compare the charge with the last hit added?
            if(ChkCharge) {
              // there is a nearby hit
              hnear = 1;
              float chgrat = multipletChg / fHits[lastClHit].Integral();
              if(prt) mf::LogVerbatim("CC")<<" merge hits charge check "
                <<(int)multipletChg<<" Previous hit charge "<<(int)fHits[lastClHit].Integral();
              if(chgrat > 2) doMerge = false;
            }
          } // doMerge && nused == 0
        } // doMerge true
      } // !fuBCode
      if(doMerge) {
        // there is a nearby hit and it will be merged
        hnear = -1;
        bool didMerge;
        MergeHits(imbest, didMerge);
      } // doMerge
    } // Hits[imbest].Multiplicity() > 1
    
    // attach to the cluster and fit
    fcl2hits.push_back(imbest);
    FitCluster();
    chifits.push_back(clChisq);
    hitNear.push_back(hnear);
    // remove the charge of the just added hit
    cnear -= fHits[imbest].Integral();
    if(cnear < 0) cnear = 0;
    // divide by the just added hit charge
    cnear /= fHits[imbest].Integral();
    chgNear.push_back(cnear);
    if(prt) mf::LogVerbatim("CC")<<" >>LADD"<<pass<<" W:T "<<kwire<<":"<<(int)fHits[imbest].PeakTime()<<std::setprecision(3)<<" clChisq "<<clChisq<<" charge "<<(int)fHits[imbest].Integral()<<" hitNear "<<hitNear[hitNear.size()-1]<<" chgNear "<<std::setprecision(2)<<cnear<<" fcl2hits size "<<fcl2hits.size();
    // decide what to do with a bad fit
    if(clChisq > fChiCut[pass]) {
      fcl2hits.pop_back();
      chifits.pop_back();
      hitNear.pop_back();
      chgNear.pop_back();
      FitCluster();
      HitOK = false;
  if(prt) mf::LogVerbatim("CC")<<" Bad fit. Removed hit. New clChisq "
    <<std::setprecision(3)<<clChisq<<" nhits "<<fcl2hits.size();
    }
      // stop tracking if previous chisq values were close to the cut.
      // this is an indicator that the track is wandering too much for this pass
      if(chifits.size() > 2 && 
        chifits[chifits.size()-2] > 0.8 * fChiCut[pass]) SigOK = false;
  if(prt) mf::LogVerbatim("CC")<<"  Set SigOK = "<<SigOK;

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
    unsigned int lastClHit = fcl2hits[fcl2hits.size()-1];
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

    unsigned int firsthit = WireHitRange[index].first;
    unsigned int lasthit = WireHitRange[index].second;

    
    // the projected time of the cluster on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    // Find the projected time error including the projection error and the
    // error from the last hit added
    float prtimerr2, nerr;
    if(fuBCode) {
      prtimerr2 = 10 * std::abs(kwire-wire0)*clparerr[1];
      nerr = 20;
    } else {
      prtimerr2 = std::abs(kwire-wire0)*clparerr[1]*clparerr[1];
      nerr = 6;
    }
    // apply an angle dependent scale factor to the hit error
    float hiterr = AngleFactor(clpar[1]) * fHitErrFac * fHits[lastClHit].RMS();
    float err = std::sqrt(prtimerr2 + hiterr * hiterr);
    // Time window for accepting a hit.
    float prtimeLo = prtime - nerr * err;
    float prtimeHi = prtime + nerr * err;
    float chgWinLo = prtime - fChgNearWindow;
    float chgWinHi = prtime + fChgNearWindow;
    if(prt) mf::LogVerbatim("CC")<<"AddHit: wire "<<kwire<<" prtime Lo "<<(int)prtimeLo<<" Hi "<<(int)prtimeHi<<" prtimerr "<<sqrt(prtimerr2)
      <<" fAveChg "<<(int)fAveChg;

    // loop through the hits
    unsigned int imbest = INT_MAX;
    float best = 9999.;
    float cnear = 0;
    for(unsigned int khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if(inClus[khit] < 0) continue;
  if(prt) mf::LogVerbatim("CC")
    <<" Chk W:T "<<kwire<<":"<<(short)fHits[khit].PeakTime()
    <<" InClus "<<inClus[khit]
    <<" mult "<<fHits[khit].Multiplicity()
    <<" RMS "<<std::setprecision(2)<<fHits[khit].RMS()
    <<" Chi2 "<<std::fixed<<std::setprecision(2)<<fHits[khit].GoodnessOfFit()
    <<" Charge "<<(int)fHits[khit].Integral()
    <<" LoT "<<(int)fHits[khit].StartTick()
    <<" HiT "<<(int)fHits[khit].EndTick()
    <<" index " << khit;
      // count charge in the window
      if(fHits[khit].StartTick() > chgWinLo && fHits[khit].EndTick() < chgWinHi) cnear += fHits[khit].Integral();
      // check for signal
      if(prtimeHi < fHits[khit].StartTick()) continue;
      if(prtimeLo > fHits[khit].EndTick()) continue;
      SigOK = true;
      // check for good hit
      if(fHits[khit].PeakTime() < prtimeLo) continue;
      if(fHits[khit].PeakTime() > prtimeHi) continue;
      // hit used?
      if(inClus[khit] > 0) continue;
      float dtime = std::abs(fHits[khit].PeakTime() - prtime);
      if(prt) mf::LogVerbatim("CC")<<" dtime "<<dtime<<" best "<<best<<" imbest "<<imbest;
      if(dtime < best) {
        best = dtime;
        imbest = khit;
      }
    } // khit
    
    if(!SigOK) {
      if(fAllowNoHitWire == 0) return;
      if(prt) mf::LogVerbatim("CC")<<" wire0 "<<wire0<<" kwire "<<kwire<<" max "<<fAllowNoHitWire<<" imbest "<<imbest;
      if((wire0 - kwire) > fAllowNoHitWire) return;
      SigOK = true;
    }

    if(imbest == INT_MAX) return;

    recob::Hit const& hit = fHits[imbest];
    
    if(prt) mf::LogVerbatim("CC")<<" Best hit time "<<(int)hit.PeakTime();
    
    short hnear = 0;
    // merge hits in a doublet?
    bool didMerge = false;
    if(fHitMergeChiCut > 0 && hit.Multiplicity() == 2) {
      bool doMerge = true;
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(std::abs(kwire - vtx[ivx].Wire) < 10 &&
           std::abs(int(hit.PeakTime() - vtx[ivx].Time)) < 20 )
        {
          doMerge = false;
          break;
        }
      } // ivx
      if (doMerge) {
        // find the neighbor hit
        unsigned int oht;
        if(hit.LocalIndex() == 0) {
          oht = imbest + 1;
        } else {
          oht = imbest - 1;
        } // hit.LocalIndex() == 0
        // check the hit time separation
        recob::Hit const& other_hit = fHits[oht];
        float hitSep = std::abs(hit.PeakTime() - other_hit.PeakTime());
        hitSep /= hit.RMS();
        // check the the charge similarity
        float totChg = hit.Integral() + other_hit.Integral();
        float lastHitChg = fAveChg;
        if(lastHitChg < 0) lastHitChg = fHits[lastClHit].Integral();
        hnear = 1;
        if(inClus[oht] == 0 && hitSep < fHitMergeChiCut && std::abs(totChg - lastHitChg) < std::abs(hit.Integral() - lastHitChg)) {
          if(prt) mf::LogVerbatim("CC")<<"Merging hit doublet "<<imbest;
          MergeHits(imbest, didMerge);
        } // not in a cluster, hitSep OK, total charge OK
      } // doMerge
    } // fHitMergeChiCut > 0 && hit.Multiplicity() == 2
    
    // Make a charge similarity cut if the average charge is defined
    bool fitChg = true;
    if(fAveChg > 0.) {

      float chgrat = (hit.Integral() - fAveChg) / fAveChg;
      if(prt) mf::LogVerbatim("CC")<<" Chgrat "<<std::setprecision(2)<<chgrat;

      // charge is way too high?
      if(chgrat > 2 * fChgCut[pass]) {
        if(prt) mf::LogVerbatim("CC")<<" fails high charge cut";
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
        if(prt) mf::LogVerbatim("CC")<<" fails low charge cut. Stop crawling.";
        SigOK = false;
        return;
      } // lasthitlow
    
      // the last hit was high charge and this one is also
      if(lasthitbig && chgrat > fChgCut[pass]) {
        if(prt) mf::LogVerbatim("CC")<<" fails 2nd high charge cut";
        return;
      } // lasthitbig

    
      // require that large charge hits have a very good projection error
      if(chgrat > fChgCut[pass]) {
        if(best > 2 * err) {
          if(prt) mf::LogVerbatim("CC")<<" high charge && bad dT= "
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
    hitNear.push_back(hnear);
    // remove the charge of the just added hit
    cnear -= fHits[imbest].Integral();
    if(cnear < 0) cnear = 0;
    // divide by the just added hit charge
    cnear /= fHits[imbest].Integral();
    chgNear.push_back(cnear);
    // nearby hit check
    ChkClusterNearbyHits(prt);
    HitOK = true;
    
    if(chgNear.size() != fcl2hits.size()) {
      mf::LogError("CC")<<"AddHit: Bad length";
      return;
    }

    if(prt) mf::LogVerbatim("CC")<<" >>ADD"<<pass<<" W:T "<<kwire<<":"<<(short)fHits[imbest].PeakTime()<<" dT "<<best<<std::setprecision(2)<<" Chisq "<<clChisq<<" Chg "<<(int)fHits[imbest].Integral()<<" AveChg "<<(int)fAveChg<<" hitNear "<<hitNear[hitNear.size()-1]<<" chgNear "<<cnear<<" fcl2hits size "<<fcl2hits.size();

    if(!fitChg) return;
  if(prt) mf::LogVerbatim("CC")<<" Fit charge ";
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
      
      if(hitNear.size() != fcl2hits.size()) {
        mf::LogWarning("CC")<<"Coding error: hitNear size != fcl2hits";
        return;
      }
      
      // Analyze the last 6 hits added but don't consider the first few hits
      if(hitNear.size() < 12) return;
      
      // TODO move into loops
      unsigned short ii, indx;
      unsigned short merged = 0;
      unsigned short notmerged = 0;
      for(ii = 0; ii < 6; ++ii) {
        indx = hitNear.size() - 1 - ii;
        if(hitNear[indx] > 0) ++notmerged;
        if(hitNear[indx] < 0) ++merged;
      }
      
//  if(prt) mf::LogVerbatim("CC")
//    <<"ChkClusterNearbyHits: nearby hits merged "<<merged
//    <<" not merged "<<notmerged;

      if(notmerged < 2) return;
      
      // a number of nearby hits were not merged while crawling, so the 
      // average charge is probably wrong. Look at the last 6 hits added
      // and merge them if they are close
      bool didMerge;
      for(ii = 0; ii < 6; ++ii) {
        indx = fcl2hits.size() - 1 - ii;
        const unsigned int iht = fcl2hits[indx];
        recob::Hit const& hit = fHits[iht];
        if(hit.Multiplicity() == 2) {
          // hit doublet. Get the index of the other hit
          unsigned int oht;
          if(hit.LocalIndex() == 0) {
//            oht = NextHitPresent(iht);
            oht = iht + 1;
          } else {
//            oht = PrevHitPresent(iht);
            oht = iht - 1;
          } // hit.LocalIndex() == 0
          recob::Hit const& other_hit = fHits[oht];
          // TODO use Hit::TimeDistanceAsRMS()
          float hitSep = std::abs(hit.PeakTime() - other_hit.PeakTime());
          hitSep /= hit.RMS();
//          if(hitSep < fHitMergeChiCut && isHitFree(oht)) {
          if(hitSep < fHitMergeChiCut && inClus[oht] == 0) {
            // check charge assuming the hits are merged
            float chgRat = hit.Integral() + other_hit.Integral();
            chgRat /= hit.Integral();
            // merge em
            if(chgRat < 4) {
      if(prt) mf::LogVerbatim("CC")<<"Merging hit doublet "<<iht;
              MergeHits(iht, didMerge);
              if(didMerge) hitNear[indx] = -1;
            }
          } // hitSep OK and not in a cluster
        } // hit doublet
      } // ii
      
      // now re-fit
      FitCluster();
      FitClusterChg();

      if(prt) mf::LogVerbatim("CC")<<"ChkClusterNearbyHits refit cluster. fAveChg= "<<fAveChg;
      
    } // ChkClusterHitNear()

//////////////////////////////////////
    void ClusterCrawlerAlg::FitVtx(unsigned short iv)
    {
      std::vector<float> x;
      std::vector<float> y;
      std::vector<float> ey2;
      float arg;

      // Set this large in case something bad happens
      vtx[iv].ChiDOF = 99;

      // make a list of clusters
      unsigned short icl;
      double darg, slp = 0.1, sumerr = 0;
      std::vector<unsigned short> vcl;
      for(icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != vtx[iv].CTP) continue;
        if(tcl[icl].EndVtx == iv) {
          vcl.push_back(icl);
          slp = fabs(tcl[icl].EndSlp);
          if(slp < 0.01) slp = 0.01;
          darg = tcl[icl].EndSlpErr / slp;
          sumerr += darg * darg;
        } // EndVtx
        if(tcl[icl].BeginVtx == iv) {
          vcl.push_back(icl);
          slp = fabs(tcl[icl].BeginSlp);
          if(slp < 0.01) slp = 0.01;
          darg = tcl[icl].BeginSlpErr / slp;
          sumerr += darg * darg;
        } // BeginVtx
      }
      
      vtx[iv].NClusters = vcl.size();
      
      if(vcl.size() == 0) return;
      
      if(vcl.size() == 1) {
        unsigned int hit;
        unsigned short indx;
        icl = vcl[0];
        // Put the vertex at the appropriate end of the cluster
        if(tcl[icl].EndVtx == iv) {
          vtx[iv].Wire = tcl[icl].EndWir;
          vtx[iv].WireErr = 1;
          vtx[iv].Time = tcl[icl].EndTim;
          // set the vertex time error to the hit error used for fitting
          indx = tcl[icl].tclhits.size() - 1;
          hit = tcl[icl].tclhits[indx];
          vtx[iv].TimeErr = fHitErrFac * fHits[hit].RMS() * fHits[hit].Multiplicity();
          vtx[iv].ChiDOF = 0;
        }
        if(tcl[icl].BeginVtx == iv) {
          vtx[iv].Wire = tcl[icl].BeginWir;
          vtx[iv].WireErr = 1;
          vtx[iv].Time = tcl[icl].BeginTim;
          // set the vertex time error to the hit error used for fitting
          hit = tcl[icl].tclhits[0];
          vtx[iv].TimeErr = fHitErrFac * fHits[hit].RMS() * fHits[hit].Multiplicity();
          vtx[iv].ChiDOF = 0;
        }
        return;
      } // size 1

      for(unsigned short ii = 0; ii < vcl.size(); ++ii) {
        icl = vcl[ii];
        if(tcl[icl].EndVtx == iv) {
          x.push_back(tcl[icl].EndSlp);
          arg = tcl[icl].EndSlp * tcl[icl].EndWir - tcl[icl].EndTim;
          y.push_back(arg);
          if(tcl[icl].EndSlpErr > 0.) {
            arg = tcl[icl].EndSlpErr * tcl[icl].EndWir;
          } else {
            arg = .1 * tcl[icl].EndWir;
          }
          ey2.push_back(arg * arg);
        } else if(tcl[icl].BeginVtx == iv) {
          x.push_back(tcl[icl].BeginSlp);
          float arg = tcl[icl].BeginSlp * tcl[icl].BeginWir - tcl[icl].BeginTim;
          y.push_back(arg);
          if(tcl[icl].BeginSlpErr > 0.) {
            arg = tcl[icl].BeginSlpErr * tcl[icl].BeginWir;
          } else {
            arg = .1 * tcl[icl].BeginWir;
          }
          ey2.push_back(arg * arg);
        }
      } // ii
      if(x.size() < 2) return;
      
      float vTime = 0.;
      float vTimeErr = 0.;
      float vWire = 0.;
      float vWireErr = 0.;
      float chiDOF;
      fLinFitAlg.LinFit(x, y, ey2, vTime, vWire, vTimeErr, vWireErr, chiDOF);
      if(chiDOF > 900) return;
      vTime = -vTime;
      // a crazy time from the fit?
      if(vTime < 0 || vTime > fMaxTime) return;
      // a crazy wire from the fit?
      geo::PlaneID iplID = DecodeCTP(vtx[iv].CTP);
      if(vWire < 0 || vWire > geom->Nwires(iplID.Plane, iplID.TPC, iplID.Cryostat)) return;
      vtx[iv].ChiDOF = chiDOF;
      vtx[iv].Wire = vWire;
      vtx[iv].Time = vTime;
      if(x.size() > 2) {
        // Use errors from the fit
        // set minimum wire error to 1/2 wire
        vtx[iv].WireErr = std::max(vWireErr, (float)1);
        // and time error to one tick
        vtx[iv].TimeErr = std::max(vTimeErr, (float)1);
      } else {
        // Use approximate errors
        sumerr = sqrt(sumerr);
        vtx[iv].WireErr = sumerr;
        vtx[iv].TimeErr = sumerr * fabs(slp);
      }
      // set minimum wire error to 1 wire
      if(vtx[iv].WireErr < 1) vtx[iv].WireErr = 1;
      // set minimum time error to one tick
      if(vtx[iv].TimeErr < 1) vtx[iv].TimeErr = 1;
      
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
        float theTime;
        short dwb, dwe;

	const dataprov::DetectorProperties* detprop = art::ServiceHandle<util::DetectorPropertiesService>()->getDetectorProperties();

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
          // Create a new 2D vertex and see how many clusters we can attach to it
          VtxStore vnew;
          vnew.Wire = theWire;
          vnew.Time = theTime;
          vnew.CTP = clCTP;
          vnew.Topo = 7;
          vtx.push_back(vnew);
          unsigned short ivnew = vtx.size() -1;
          std::vector<short> vclIndex;
          for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
            if(tcl[icl].ID < 0) continue;
            if(tcl[icl].CTP != clCTP) continue;
            dwb = std::abs(theWire - tcl[icl].BeginWir);
            dwe = std::abs(theWire - tcl[icl].EndWir);
            // rough cut to start
            if(dwb > 10 && dwe > 10) continue;
            if(dwb < dwe && dwb < 10 && tcl[icl].BeginVtx < 0) {
              // cluster begin is closer
              if(theWire < tcl[icl].BeginWir + 5) continue;
              if(ClusterVertexChi(icl, 0, ivnew) > fVertex3DCut) continue;
              tcl[icl].BeginVtx = ivnew;
              vclIndex.push_back(icl);
            } else if(dwe < 10 && tcl[icl].EndVtx < 0) {
              // cluster end is closer
              if(theWire > tcl[icl].EndWir - 5) continue;
              if(ClusterVertexChi(icl, 1, ivnew) > fVertex3DCut) continue;
              tcl[icl].EndVtx = ivnew;
              vclIndex.push_back(icl);
            } // dwb/dwe check
          } // icl
          bool goodVtx = false;
          if(vclIndex.size() > 0) {
            FitVtx(ivnew);
            goodVtx = (vtx[ivnew].ChiDOF < fVertex3DCut);
            vtx3[ivx].Ptr2D[thePlane] = ivnew;
          }
          if(goodVtx) {
            vtx3[ivx].Ptr2D[thePlane] = ivnew;
            vtx3[ivx].Wire = -1;
          } else {
            // clobber the vertex
            vtx.pop_back();
            for(unsigned short ii = 0; ii < vclIndex.size(); ++ii) {
              unsigned short icl = vclIndex[ii];
              if(tcl[icl].BeginVtx == ivnew) tcl[icl].BeginVtx = -99;
              if(tcl[icl].EndVtx == ivnew) tcl[icl].EndVtx = -99;
            } // ii
          }
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

	const dataprov::DetectorProperties* detprop = art::ServiceHandle<util::DetectorPropertiesService>()->getDetectorProperties();

        for(unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
          if(vtx3[ivx].CStat != cstat || vtx3[ivx].TPC != tpc) continue;
          // Complete 3D vertex with matching 2D vertices in all planes?
  if(vtxprt) mf::LogVerbatim("CC")<<"Vtx3ClusterSplit ivx "<<ivx
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
          // make a list of clusters that have hits near this point on nearby wires
          std::vector<short> clIDs;
          if(theWire > fFirstWire + 5) { loWire = theWire - 5; } else { loWire = fFirstWire; }
          if(theWire < fLastWire  - 5) { hiWire = theWire + 5; } else { hiWire = fLastWire; }
    if(vtxprt) mf::LogVerbatim("CC")<<"3DVtx "<<ivx
      <<" look for cluster hits near P:W:T "<<thePlane<<":"<<theWire<<":"<<(int)theTime
      <<" Wire range "<<loWire<<" to "<<hiWire;
          for(unsigned short wire = loWire; wire < hiWire; ++wire) {
            unsigned short index = wire - fFirstWire;
            // ignore dead wires or wires with no hits
            if(WireHitRange[index].first < 0) continue;
            unsigned int firsthit = WireHitRange[index].first;
            unsigned int lasthit = WireHitRange[index].second;
            for(unsigned int khit = firsthit; khit < lasthit; ++khit) {
              // ignore obsolete and un-assigned hits
              if(inClus[khit] <= 0) continue;
              if((unsigned short)inClus[khit] > tcl.size() + 1) {
                mf::LogError("CC")<<"Invalid hit InClus. "<<khit<<" "<<inClus[khit];
                continue;
              }
              // check an expanded time range
              if(theTime < fHits[khit].StartTick() - 10) continue;
              if(theTime > fHits[khit].EndTick() + 10) continue;
              kclID = inClus[khit];
              kcl = kclID - 1;
              // ignore obsolete clusters
              if(tcl[kcl].ID < 0) continue;
              // ignore short clusters
              if(tcl[kcl].tclhits.size() < 6) continue;
              
              // put the cluster in the list if it's not there already
    if(vtxprt) mf::LogVerbatim("CC")<<"Bingo "<<ivx<<" plane "<<thePlane
      <<" wire "<<wire<<" hit "<<fHits[khit].WireID().Wire<<":"<<(int)fHits[khit].PeakTime()
        <<" inClus "<<inClus[khit]
        <<" P:W:T "<<thePlane<<":"<<tcl[kcl].BeginWir<<":"<<(int)tcl[kcl].BeginTim;
              if(std::find(clIDs.begin(), clIDs.end(), kclID) == clIDs.end()) {
                // ignore long straight clusters
                if(tcl[kcl].tclhits.size() > 100 ) {
                  dth = tcl[kcl].BeginAng - tcl[kcl].EndAng;
    if(vtxprt) mf::LogVerbatim("CC")<<"Long straight check: nhits "
      <<tcl[kcl].tclhits.size()<<" dth "<<dth;
                  if(std::abs(dth) < 0.05) continue;
                } // tcl[kcl].tclhits.size() > 100
                clIDs.push_back(kclID);
              } // std::find
            } // khit
          } // wire
          if(clIDs.size() == 0) continue;
    if(vtxprt) {
      for(unsigned int ii = 0; ii < clIDs.size(); ++ii) mf::LogVerbatim("CC")<<" clIDs "<<clIDs[ii];
    }

          unsigned short ii, icl, jj;
          unsigned int iht;
          short nhitfit;
          bool didit;
          // find a reasonable time error using the 2D vertices that comprise this
          // incomplete 3D vertex
          float tErr = 1;
          unsigned short i2Dvx = 0;
          for(ii = 0; ii < 3; ++ii) {
            if(ii == thePlane) continue;
            i2Dvx = vtx3[ivx].Ptr2D[ii];
            if(i2Dvx > vtx.size() - 1) {
              mf::LogError("CC")<<"Vtx3ClusterSplit: Coding error";
              return;
            }
            if(vtx[i2Dvx].TimeErr > tErr) tErr = vtx[i2Dvx].TimeErr;
          } // ii -> plane

          // do a local fit near the crossing point and make a tighter cut
          for(ii = 0; ii < clIDs.size(); ++ii) {
            icl = clIDs[ii] - 1;
            didit = false;
            for(jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
              iht = tcl[icl].tclhits[jj];
              if(fHits[iht].WireID().Wire < theWire) {
                nhitfit = 3;
                if(jj > 3) nhitfit = -3;
                FitClusterMid(icl, iht, nhitfit);
                float doca = DoCA(-1, 1, theWire, theTime);
    if(vtxprt) mf::LogVerbatim("CC")<<" cls "<<icl<<" dt "<<dth<<" DoCA "<<doca<<" tErr "<<tErr;
                if((doca / tErr) > 2) clIDs[ii] = -1;
                didit = true;
                break;
              } // fHits[iht].WireID().Wire < theWire
              if(didit) break;
            } // jj
            if(didit) break;
          } // ii
    if(vtxprt) {
      mf::LogVerbatim("CC")<<"clIDs after fit "<<clIDs.size();
      for(ii = 0; ii < clIDs.size(); ++ii) mf::LogVerbatim("CC")<<" clIDs "<<clIDs[ii];
    }

          // see if any candidates remain
          unsigned short nok = 0;
          for(ii = 0; ii < clIDs.size(); ++ii) if(clIDs[ii] >= 0) ++nok;
          if(nok == 0) continue;

          // make a new 2D vertex
          VtxStore vnew;
          vnew.Wire = theWire;
          vnew.WireErr = 1;
          vnew.Time = theTime;
          vnew.TimeErr = 1;
          vnew.Topo = 8;
          vnew.CTP = clCTP;
          vtx.push_back(vnew);
          // update the 2D -> 3D vertex pointer
          unsigned short ivnew = vtx.size() - 1;
          if(vtxprt) mf::LogVerbatim("CC")<<"Make new 2D vtx "<<ivnew<<" in plane "<<thePlane<<" from 3D vtx "<<ivx;
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
  if(vtxprt) mf::LogVerbatim("CC")<<"Attach to Begin "<<icl;
            } else if(pos > tcl[icl].tclhits.size()) {
              // vertex is US of the cluster Eend
              tcl[icl].EndVtx = ivnew;
  if(vtxprt) mf::LogVerbatim("CC")<<"Attach to End "<<icl;
            } else {
              // vertex is in the middle of the cluster
              if(vtxprt) mf::LogVerbatim("CC")<<"Split cluster "<<clIDs[ii]<<" at pos "<<pos;
              if(!SplitCluster(icl, pos, ivnew)) {
                if(vtxprt) mf::LogVerbatim("CC")<<"SplitCluster failed";
                continue;
              }
              tcl[icl].ProcCode += 10000;
              tcl[tcl.size()-1].ProcCode += 10000;
            } // pos check
          } // ii
          // Fit the vertex position
          FitVtx(ivnew);
          if(vtx[ivnew].ChiDOF < 5 && vtx[ivnew].WireErr < 2) {
            // mark the 3D vertex as complete
            vtx3[ivx].Wire = -1;
          } else {
  if(vtxprt) mf::LogVerbatim("CC")<<"Bad vtx fit "<<ivnew<<" ChiDOF "<<vtx[ivnew].ChiDOF<<" WireErr "<<vtx[ivnew].WireErr;
            // Recover (partially) from a bad fit. Leave the ProcCode as-is to trace this problem
            vtx.pop_back();
            vtx3[ivx].Ptr2D[thePlane] = -1;
            // find the cluster - vertex associations
            for(jj = 0; jj < tcl.size(); ++jj) {
              if(tcl[jj].BeginVtx == ivnew) tcl[jj].BeginVtx = -99;
              if(tcl[jj].EndVtx == ivnew) tcl[jj].EndVtx = -99;
            } // jj
          }
        } // ivx
        
      } // Vtx3ClusterSplit()

  
  //////////////////////////////////////
		void ClusterCrawlerAlg::FindHammerClusters()
  {
    // look for a long cluster that stops at a short cluster in two views. This can occur in a CCmu
    // interaction where two protons are emitted back-to-back and are therefore reconstructed as one cluster
    // This routine looks for this signature, and if found, splits the short clusters and creates a new 3D vertex.
    // This routine only considers the case where the long cluster intersects the short cluster at the US (End) end.
    
    unsigned short nPln = geom->Cryostat(cstat).TPC(tpc).Nplanes();
    if(nPln != 3) return;
    
    float ew1, ew2, bw2, fvw;
    
    struct Hammer {
      bool Used;
      unsigned int Wire;  // intersection point of the long cluster and the short cluster
      float Tick;  // intersection point of the long cluster and the short cluster
      float X;
      unsigned short longClIndex;
      unsigned short shortClIndex;
      unsigned short splitPos;
    };
    std::array< std::vector<Hammer>, 3> hamrVec;
    
    const dataprov::DetectorProperties* detprop = art::ServiceHandle<util::DetectorPropertiesService>()->getDetectorProperties();

    unsigned short ipl;
    bool useit = false;
    for(ipl = 0; ipl < 3; ++ipl) {
      clCTP = EncodeCTP(cstat, tpc, ipl);
      for(unsigned short ic1 = 0; ic1 < tcl.size(); ++ic1) {
        if(tcl[ic1].ID < 0) continue;
        // require a long cluster
        if(tcl[ic1].tclhits.size() < 20) continue;
        if(tcl[ic1].CTP != clCTP) continue;
        // ignore long clusters with an End vertex assignment
        if(tcl[ic1].EndVtx >= 0) continue;
        ew1 = tcl[ic1].EndWir;
        for(unsigned short ic2 = 0; ic2 < tcl.size(); ++ic2) {
          if(tcl[ic2].ID < 0) continue;
          // require a short cluster
          if(tcl[ic2].tclhits.size() > 20) continue;
          // but not too short cluster
          if(tcl[ic2].tclhits.size() < 6) continue;
          if(tcl[ic2].CTP != clCTP) continue;
          ew2 = tcl[ic2].EndWir;
          bw2 = tcl[ic2].BeginWir;
          // check the US end. The End of the long cluster must lie between the Begin and End wire of the
          // short cluster
          if(ew1 < ew2 || ew1 > bw2 ) continue;
          // look for intersection of the two clusters
          float best = 10;
          short ibst = -1;
          unsigned short spos = 0;
          for(unsigned short ii = 0; ii < tcl[ic2].tclhits.size(); ++ii) {
            unsigned int iht = tcl[ic2].tclhits[ii];
            float dw = fHits[iht].WireID().Wire - tcl[ic1].EndWir;
            float dt = fabs(fHits[iht].PeakTime() - tcl[ic1].EndTim - tcl[ic1].EndSlp * dw);
            if(dt < best) {
              best = dt;
              ibst = iht;
              spos = ii;
            }
          } // ii
          if(ibst < 0) continue;
          fvw = 0.5 + fHits[ibst].WireID().Wire;
          Hammer aHam;
          aHam.Used = false;
          aHam.Wire = (0.5 + fvw);
          aHam.Tick = fHits[ibst].PeakTime();
          aHam.X = detprop->ConvertTicksToX((double)aHam.Tick, (int)ipl, (int)tpc, (int)cstat);
//          std::cout<<"aHam "<<" fvw "<<fvw<<" X "<<aHam.X<<" ic1 ID "<<tcl[ic1].ID<<" ic2 ID "<<tcl[ic2].ID<<"\n";
          aHam.longClIndex = ic1;
          aHam.shortClIndex = ic2;
          aHam.splitPos = spos;
          unsigned short indx = hamrVec[ipl].size();
          hamrVec[ipl].resize(indx + 1);
          hamrVec[ipl][indx] = aHam;
//          hamrVec[ipl].push_back(aHam);
          useit = true;
        } // ic2
        if(useit) break;
      } // ic1
    } // ipl
    
    unsigned short noham = 0;
    for(ipl = 0; ipl < 3; ++ipl) if(hamrVec[ipl].size() == 0) ++noham;
    if(noham > 1) return;
    
    // Y,Z limits of the detector
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    
    const geo::TPCGeo &thetpc = geom->TPC(tpc, cstat);
    thetpc.LocalToWorld(local,world);
    float YLo = world[1]-geom->DetHalfHeight(tpc,cstat) + 1;
    float YHi = world[1]+geom->DetHalfHeight(tpc,cstat) - 1;
    float ZLo = world[2]-geom->DetLength(tpc,cstat)/2 + 1;
    float ZHi = world[2]+geom->DetLength(tpc,cstat)/2 - 1;

    // Match in 3D
    float dX;
    double y, z;
    unsigned short icl, jpl, jcl, kpl, splitPos;
    for(ipl = 0; ipl < 3; ++ipl) {
      if(hamrVec[ipl].size() == 0) continue;
      jpl = (ipl + 1) % nPln;
      kpl = (jpl + 1) % nPln;
      for(unsigned short ii = 0; ii < hamrVec[ipl].size(); ++ii) {
        if(hamrVec[ipl][ii].Used) continue;
        for(unsigned short jj = 0; jj < hamrVec[jpl].size(); ++jj) {
          if(hamrVec[jpl][jj].Used) continue;
          dX = hamrVec[ipl][ii].X - hamrVec[jpl][jj].X;
          if(fabs(dX) > fVertex3DCut) continue;
          geom->IntersectionPoint(hamrVec[ipl][ii].Wire, hamrVec[jpl][jj].Wire, ipl, jpl, cstat, tpc, y, z);
          if(y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
          // mark them used
          hamrVec[ipl][ii].Used = true;
          hamrVec[jpl][jj].Used = true;
          // make a new 3D vertex
          Vtx3Store newVtx3;
          newVtx3.ProcCode = 7;
          newVtx3.X = 0.5 * (hamrVec[ipl][ii].X + hamrVec[jpl][jj].X);
          // TODO: do this correctly;
          newVtx3.XErr = fabs(hamrVec[ipl][ii].X - hamrVec[jpl][jj].X);
          newVtx3.Y = y;
          newVtx3.YErr = 1; // TODO
          newVtx3.Z = z;
          newVtx3.ZErr = 1;  // TODO
          newVtx3.CStat = cstat;
          newVtx3.TPC = tpc;
          
          // make 2D vertex in ipl
          VtxStore newVtx2;
          newVtx2.Wire = hamrVec[ipl][ii].Wire;
          newVtx2.WireErr = 2;
          newVtx2.Time = hamrVec[ipl][ii].Tick;
          newVtx2.TimeErr = 5;
          newVtx2.Topo = 6;
          icl = hamrVec[ipl][ii].longClIndex;
          newVtx2.CTP = tcl[icl].CTP;
          vtx.push_back(newVtx2);
          unsigned short ivnew = vtx.size() - 1;
          // associate the new vertex with the long cluster
          tcl[icl].EndVtx = ivnew;
          FitVtx(ivnew);
          // stash the index in newVtx3
          newVtx3.Ptr2D[ipl] = (short)ivnew;
          // split the short cluster and associate the new clusters with the new vtx
          icl = hamrVec[ipl][ii].shortClIndex;
          splitPos = hamrVec[ipl][ii].splitPos;
          if(!SplitCluster(icl, splitPos, ivnew)) return;
          
          // make 2D vertex in jpl
          newVtx2.Wire = hamrVec[jpl][jj].Wire;
          newVtx2.Time = hamrVec[jpl][jj].Tick;
          newVtx2.Topo = 6;
          jcl = hamrVec[jpl][jj].longClIndex;
          newVtx2.CTP = tcl[jcl].CTP;
          vtx.push_back(newVtx2);
          ivnew = vtx.size() - 1;
//          std::cout<<"newVtx2 jpl index "<<vtx.size()-1<<" CTP "<<newVtx2.CTP<<"\n";
          // associate the new vertex with the long cluster
          tcl[jcl].EndVtx = ivnew;
          // stash the index in newVtx3
          newVtx3.Ptr2D[jpl] = (short)(vtx.size() - 1);
          // split the short cluster and associate the new clusters with the new vtx
          jcl = hamrVec[jpl][jj].shortClIndex;
          splitPos = hamrVec[jpl][jj].splitPos;
          if(!SplitCluster(jcl, splitPos, vtx.size() - 1)) return;
          FitVtx(ivnew);
//          std::cout<<"Split jpl "<<jpl<<" cl ID "<<tcl[jcl].ID<<" pos "<<splitPos<<" vtx "<<vtx.size()-1<<"\n";
          
          // set the kpl 2D vertex index < 0. Let follow-on code find the 3rd plane vertex
          newVtx3.Ptr2D[kpl] = -1;
          double WPos[3] = {0, y, z};
          newVtx3.Wire = geom->NearestWire(WPos, kpl, tpc, cstat);
//          std::cout<<"3D Vtx "<<ipl<<":"<<ii<<" "<<jpl<<":"<<jj<<" kpl "<<kpl<<" wire "<<newVtx3.Wire<<"\n";
          vtx3.push_back(newVtx3);
        } // jj
      } // ii
    }
    
  } // FindHammerClusters


//////////////////////////////////////
    void ClusterCrawlerAlg::VtxMatch(geo::TPCID const& tpcid)
    {
      // Create 3D vertices from 2D vertices. 3D vertices that are matched
      // in all three planes have Ptr2D >= 0 for all planes
      
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      
      const dataprov::DetectorProperties* detprop = art::ServiceHandle<util::DetectorPropertiesService>()->getDetectorProperties();
      
      const unsigned int cstat = tpcid.Cryostat;
      const unsigned int tpc = tpcid.TPC;
      
      // Y,Z limits of the detector
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      
      const geo::TPCGeo &thetpc = geom->TPC(tpc, cstat);
      thetpc.LocalToWorld(local,world);
      // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
      float YLo = world[1]-geom->DetHalfHeight(tpc,cstat) + 1;
      float YHi = world[1]+geom->DetHalfHeight(tpc,cstat) - 1;
      float ZLo = world[2]-geom->DetLength(tpc,cstat)/2 + 1;
      float ZHi = world[2]+geom->DetLength(tpc,cstat)/2 - 1;
      
      vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);
      
      // wire spacing in cm
      float wirePitch = geom->WirePitch(0, 1, 0, tpcid.TPC, tpcid.Cryostat);
      
      // fill temp vectors of 2D vertex X and X errors
      std::vector<float> vX(vtx.size());
      std::vector<float> vXsigma(vtx.size());
      float vXp;
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(vtx[ivx].NClusters == 0) continue;
        geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
        if(iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
        // Convert 2D vertex time error to X error
        vX[ivx]  = detprop->ConvertTicksToX((double)vtx[ivx].Time, (int)iplID.Plane, (int)tpc, (int)cstat);
        vXp = detprop->ConvertTicksToX((double)(vtx[ivx].Time + vtx[ivx].TimeErr), (int)iplID.Plane, (int)tpc, (int)cstat);
        vXsigma[ivx] = fabs(vXp - vX[ivx]);
      } // ivx
      
      // create a array/vector of 2D vertex indices in each plane
      std::array<std::vector<unsigned short>, 3> vIndex;
      unsigned short indx, ipl;
      for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if(vtx[ivx].NClusters == 0) continue;
        geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
        if(iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
        ipl = iplID.Plane;
        if(ipl > 2) continue;
        indx = vIndex[ipl].size();
        vIndex[ipl].resize(indx + 1);
        vIndex[ipl][indx] = ivx;
      }
      
      // vector of 2D vertices -> 3D vertices.
      std::vector<short> vPtr;
      for(unsigned short ii = 0; ii < vtx.size(); ++ii) vPtr.push_back(-1);
      
      // temp vector of all 2D vertex matches
      std::vector<Vtx3Store> v3temp;
      
      double y = 0, z = 0;
      TVector3 WPos = {0, 0, 0};
      // i, j, k indicates 3 different wire planes
      unsigned short ii, jpl, jj, kpl, kk, ivx, jvx, kvx;
      unsigned int iWire, jWire;
      unsigned short v3dBest = 0;
      float xbest = 0, ybest = 0, zbest = 0;
      float kX, kWire;
      // compare vertices in each view
      bool gotit = false;
      for(ipl = 0; ipl < 2; ++ipl) {
        for(ii = 0; ii < vIndex[ipl].size(); ++ii) {
          ivx = vIndex[ipl][ii];
          if(ivx > vtx.size() - 1) {
            mf::LogError("CC")<<"VtxMatch: bad ivx "<<ivx;
            return;
          }
          // vertex has been matched already
          if(vPtr[ivx] >= 0) continue;
          iWire = vtx[ivx].Wire;
          float best = fVertex3DCut;
          // temp array of 2D vertex indices in each plane
          std::array<short, 3> t2dIndex = {-1, -1, -1};
          std::array<short, 3> tmpIndex = {-1, -1, -1};
          for(jpl = ipl + 1; jpl < 3; ++jpl) {
            for(jj = 0; jj < vIndex[jpl].size(); ++jj) {
              jvx = vIndex[jpl][jj];
              if(jvx > vtx.size() - 1) {
                mf::LogError("CC")<<"VtxMatch: bad jvx "<<jvx;
                return;
              }
              // vertex has been matched already
              if(vPtr[jvx] >= 0) continue;
//              jWire = (0.5 + vtx[jvx].Wire);
              jWire = vtx[jvx].Wire;
              // new stuff
              float dX = fabs(vX[ivx] - vX[jvx]);
              float dXSigma = sqrt(vXsigma[ivx] * vXsigma[ivx] + vXsigma[jvx] * vXsigma[jvx]);
              float dXChi = dX / dXSigma;
              
              if(vtxprt) mf::LogVerbatim("CC")<<"VtxMatch: ipl "<<ipl<<" ivx "<<ivx
                <<" jpl "<<jpl<<" jvx "<<jvx<<" W:T "<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time<<" dXChi "<<dXChi;
              
              if(dXChi > fVertex3DCut) continue;
              geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
              if(y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
              WPos[1] = y;
              WPos[2] = z;
              kpl = 3 - ipl - jpl;
              kX = 0.5 * (vX[ivx] + vX[jvx]);
              kWire = -1;
              if(TPC.Nplanes() > 2) kWire = geom->NearestWire(WPos, kpl, tpc, cstat);
              kpl = 3 - ipl - jpl;
              // save this incomplete 3D vertex
              Vtx3Store v3d;
              v3d.ProcCode = 1;
              tmpIndex[ipl] = ivx;
              tmpIndex[jpl] = jvx;
              tmpIndex[kpl] = -1;
              v3d.Ptr2D = tmpIndex;
              v3d.X = kX;
              v3d.XErr = dXSigma;
              v3d.Y = y;
              float yzSigma = wirePitch * sqrt(vtx[ivx].WireErr * vtx[ivx].WireErr + vtx[jvx].WireErr * vtx[jvx].WireErr);
              v3d.YErr = yzSigma;
              v3d.Z = z;
              v3d.ZErr = yzSigma;
              v3d.Wire = kWire;
              v3d.CStat = cstat;
              v3d.TPC = tpc;
              v3temp.push_back(v3d);

  if(vtxprt) mf::LogVerbatim("CC")
    <<"VtxMatch: 2 Plane match ivx "<<ivx<<" P:W:T "<<ipl<<":"<<(int)vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time
    <<" jvx "<<jvx<<" P:W:T "<<jpl<<":"<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time<<" dXChi "<<dXChi<<" yzSigma "<<yzSigma;

              if(TPC.Nplanes() == 2) continue;
              // look for a 3 plane match
              best = fVertex3DCut;
              for(kk = 0; kk < vIndex[kpl].size(); ++kk) {
                kvx = vIndex[kpl][kk];
                if(vPtr[kvx] >= 0) continue;
//                float kvxX = detprop->ConvertTicksToX((double)vtx[kvx].Time, (int)kpl, (int)tpc, (int)cstat);
                // Wire difference error
                float dW = wirePitch * (vtx[kvx].Wire - kWire) / yzSigma;
                // X difference error
                float dX = (vX[kvx] - kX) / dXSigma;
                float kChi = 0.5 * sqrt(dW * dW + dX * dX);
                if(kChi < best) {
                  best = kChi;
                  xbest = (vX[kvx] + 2 * kX) / 3;
                  ybest = y;
                  zbest = z;
                  t2dIndex[ipl] = ivx;
                  t2dIndex[jpl] = jvx;
                  t2dIndex[kpl] = kvx;
                  v3dBest = v3temp.size() - 1;
                }

  if(vtxprt) mf::LogVerbatim("CC")<<" kvx "<<kvx<<" kpl "<<kpl
    <<" wire "<<(int)vtx[kvx].Wire<<" kTime "<<(int)vtx[kvx].Time<<" kChi "<<kChi<<" best "<<best<<" dW "<<vtx[kvx].Wire - kWire;

              } // kk
              if(vtxprt) mf::LogVerbatim("CC")<<" done best = "<<best<<" fVertex3DCut "<<fVertex3DCut;
              if(TPC.Nplanes() > 2 && best < fVertex3DCut) {
                // create a real 3D vertex using the previously entered incomplete 3D vertex as a template
                if(v3dBest > v3temp.size() - 1) {
                  mf::LogError("CC")<<"VtxMatch: bad v3dBest "<<v3dBest;
                  return;
                }
                Vtx3Store v3d = v3temp[v3dBest];
                v3d.Ptr2D = t2dIndex;
                v3d.Wire = -1;
                // TODO need to average ybest and zbest here with error weighting
                v3d.X = xbest;
                v3d.Y = ybest;
                v3d.Z = zbest;
                vtx3.push_back(v3d);
                gotit = true;
                // mark the 2D vertices as used
                for(unsigned short jj = 0; jj < 3; ++jj) if(t2dIndex[jj] >= 0) vPtr[t2dIndex[jj]] = vtx3.size() - 1;
                
                if(vtxprt) mf::LogVerbatim("CC")<<"New 3D vtx "<<vtx3.size()<<" X "<<v3d.X<<" Y "<<v3d.Y<<" Z "<<v3d.Z
                  <<" t2dIndex "<<t2dIndex[0]<<" "<<t2dIndex[1]<<" "<<t2dIndex[2]<<" best Chi "<<best;
                
              } // best < dRCut
              if(gotit) break;
            } // jj
            if(gotit) break;
          } // jpl
          if(gotit) break;
        } // ii
      } // ipl
      
      // Store incomplete 3D vertices but ignore those that are part of a complete 3D vertex
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
      
      // Modify Ptr2D for 2-plane detector
      if(TPC.Nplanes() == 2) {
        for(unsigned short iv3 = 0; iv3 < vtx3.size(); ++iv3) {
          vtx3[iv3].Ptr2D[2] = 666;
        } //iv3
      } // 2 planes
      
  if(vtxprt) {
    for(unsigned short it = 0; it < vtx3.size(); ++it) {
      mf::LogVerbatim("CC")<<"vtx3 "<<it<<" Ptr2D "<<vtx3[it].Ptr2D[0]<<" "<<vtx3[it].Ptr2D[1]<<" "<<vtx3[it].Ptr2D[2]
        <<" wire "<<vtx3[it].Wire;
    }
  }

    } // VtxMatch
  
//////////////////////////////////
    void ClusterCrawlerAlg::GetHitRange(CTP_t CTP, std::vector< std::pair<int, int> >& WireHitRange,
      unsigned short& firstwire, unsigned short& lastwire)
    {
      // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
      fFirstWire = fLastWire = 0;
			fFirstHit = 0;
			unsigned int lastHit = 0;
      WireHitRange.clear();
      bool first = true;
      geo::PlaneID planeID = DecodeCTP(clCTP);
      // find the first and last wire with a hit
      for(unsigned int hit = 0; hit < fHits.size(); ++hit) {
        // TODO change to planeID()
        if(fHits[hit].WireID().Plane != planeID.Plane) continue;
        if(fHits[hit].WireID().TPC != planeID.TPC) continue;
        if(fHits[hit].WireID().Cryostat != planeID.Cryostat) continue;
        unsigned short theWireNum = fHits[hit].WireID().Wire;
        if(first) {
          fFirstHit = hit;
          fFirstWire = theWireNum;
          first = false;
        }
        fLastWire = theWireNum;
        lastHit = hit;
      } //hit
      
      if (first) return; // we collected nothing

      // now we can define the WireHitRange vector.
      // start by defining the "no hits on wire" condition
      short sflag = -2;
      for(unsigned short wire = fFirstWire; wire <= fLastWire; ++wire) {
        WireHitRange.push_back(std::make_pair(sflag, sflag));
      }
      // overwrite with the "dead wires" condition
      filter::ChannelFilter cf;
      sflag = -1;
      for(unsigned short wire = fFirstWire+1; wire < fLastWire; ++wire) {
        raw::ChannelID_t chan = geom->PlaneWireToChannel
          ((int)planeID.Plane,(int)wire,(int)planeID.TPC,(int)planeID.Cryostat);
        // remember to offset references to WireHitRange by the FirstWire
        unsigned short indx = wire - fFirstWire;
        if(cf.BadChannel(chan)) WireHitRange[indx] = std::make_pair(sflag, sflag);
      }
      fLastWire = fFirstWire;
      unsigned int thishit = fFirstHit;
      unsigned int lastfirsthit = fFirstHit;
      // next overwrite with the index of the first/last hit on each wire
      for(unsigned int hit = fFirstHit; hit <= lastHit; ++hit) {
        recob::Hit const& theHit = fHits[hit];
        if(theHit.WireID().Plane != planeID.Plane) continue;
        if(theHit.WireID().TPC != planeID.TPC) continue;
        if(theHit.WireID().Cryostat != planeID.Cryostat) continue;
	
        unsigned short thiswire = theHit.WireID().Wire;
        if(thiswire > fLastWire) {
          unsigned short indx = fLastWire - fFirstWire;
          int itmp1 = lastfirsthit;
          int itmp2 = thishit;
          if (itmp1<0||itmp2<0) {
            throw art::Exception(art::errors::LogicError)<<"ClusterCrawlerAlg::GetHitRange() indx= "<<indx<<" itmp1= "<<itmp1<<" itmp2= "<<itmp2<<" lastfirsthit= "<<lastfirsthit<<" thishit= "<<thishit;
          }
          WireHitRange[indx] = std::make_pair(itmp1,itmp2);
          fLastWire = thiswire;
          lastfirsthit = thishit;
        } else if(thiswire < fLastWire) {
          mf::LogError("CC")<<"ERROR: Hits not sorted!!";
          return;
        }
        ++thishit;
      } //hit
      // define for the last wire
      unsigned short indx = fLastWire - fFirstWire;
      int itmp1 = lastfirsthit;
      int itmp2 = thishit;
      WireHitRange[indx] = std::make_pair(itmp1,itmp2);
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
    
    
    bool ClusterCrawlerAlg::CheckHitDuplicates
      (std::string location, std::string marker /* = "" */) const
    {
      // currently unused, only for debug
      unsigned int nDuplicates = 0;
      std::set<unsigned int> hits;
      for (unsigned int hit: fcl2hits) {
        if (hits.count(hit)) {
          ++nDuplicates;
          mf::LogProblem log("CC");
          log << "Hit #" << hit
            << " being included twice in the future cluster (ID="
            << (tcl.size() + 1) << "?) at location: " << location;
         if (!marker.empty()) log << " (marker: '" << marker << "')";
        }
        hits.insert(hit);
      } // for
      return nDuplicates > 0;
    } // ClusterCrawlerAlg::CheckHitDuplicates()

  /////////////////////////////////////////
  float ClusterCrawlerAlg::DoCA(short icl, unsigned short end, float vwire, float vtick)
  {
    // Find the Distance of Closest Approach betweeen a cluster and a point (vwire, vtick). The
    // DoCA is returned in Tick units.
    
    if(icl > (short)tcl.size()) return 9999;
    
    float cwire, cslp, ctick;
    // figure out which cluster to use
    if(icl < 0) {
      if(fcl2hits.size() == 0) return 9999;
      // cluster under construction
      if(end == 0) {
        cwire = clBeginWir;
        cslp = clBeginSlp;
        ctick = clBeginTim;
      } else {
        cwire = clpar[2];
        cslp = clpar[1];
        ctick = clpar[0];
      } // end
    } else {
      // tcl cluster
      if(end == 0) {
        cwire = tcl[icl].BeginWir;
        cslp = tcl[icl].BeginSlp;
        ctick = tcl[icl].BeginTim;
      } else {
        cwire = tcl[icl].EndWir;
        cslp = tcl[icl].EndSlp;
        ctick = tcl[icl].EndTim;
      } // end
    }
    
    // Closest approach wire
    float docaW = (vwire + cslp * (vtick - ctick) + cwire * cslp * cslp) / (1 + cslp * cslp);
    float dW = docaW - vwire;
    float dT = ctick + (vwire - cwire) * cslp - vtick;
    return sqrt(dW * dW + dT * dT);
    
  } // DoCA
  
  /////////////////////////////////////////
  float ClusterCrawlerAlg::ClusterVertexChi(short icl, unsigned short end, unsigned short ivx)
  {
    // Returns the chisq/DOF between a cluster and a vertex
    
    if(icl > (short)tcl.size()) return 9999;
    if(ivx > vtx.size()) return 9999;
    
    float cwire, cslp, cslpErr, ctick;
    // figure out which cluster to use
    if(icl < 0) {
      if(fcl2hits.size() == 0) return 9999;
      // cluster under construction
      if(end == 0) {
        cwire = clBeginWir;
        cslp = clBeginSlp;
        cslpErr = clBeginSlpErr;
        ctick = clBeginTim;
      } else {
        cwire = clpar[2];
        cslp = clpar[1];
        cslpErr = clparerr[1];
        ctick = clpar[0];
      } // end
    } else {
      // tcl cluster
      if(end == 0) {
        cwire = tcl[icl].BeginWir;
        cslp = tcl[icl].BeginSlp;
        cslpErr = tcl[icl].BeginSlpErr;
        ctick = tcl[icl].BeginTim;
      } else {
        cwire = tcl[icl].EndWir;
        cslp = tcl[icl].EndSlp;
        cslpErr = tcl[icl].EndSlpErr;
        ctick = tcl[icl].EndTim;
      } // end
    }
    
    // Closest approach wire
    float docaW = (vtx[ivx].Wire + cslp * (vtx[ivx].Time - ctick) + cwire * cslp * cslp) / (1 + cslp * cslp);
    float dW = docaW - vtx[ivx].Wire;
    float chi = dW / vtx[ivx].WireErr;
    float totChi = chi * chi;
    dW = vtx[ivx].Wire - cwire;
    float dT = ctick + dW * cslp - vtx[ivx].Time;
    if(cslpErr < 1E-3) cslpErr = 1E-3;
    // increase slope error for large angle clusters
    cslpErr *= AngleFactor(cslp);
    // cluster slope projection error
    float dTErr = dW * cslpErr;
    // squared
    dTErr *= dTErr;
    // add the vertex time error^2 to the cluster projection error^2
    dTErr += vtx[ivx].TimeErr * vtx[ivx].TimeErr;
    if(dTErr < 1E-3) dTErr = 1E-3;
    totChi += dT * dT / dTErr;
    totChi /= 2;
    
    return totChi;
    
  } // ClusterVertexChi
  
  /////////////////////////////////////////
  float ClusterCrawlerAlg::PointVertexChi(float wire, float tick, unsigned short ivx)
  {
    // Returns the Chisq/DOF between a (Wire, Tick) point and a vertex
    
    if(ivx > vtx.size()) return 9999;
    
    float dW = wire - vtx[ivx].Wire;
    float chi = dW / vtx[ivx].WireErr;
    float totChi = chi * chi;
    float dT = tick - vtx[ivx].Time;
    chi = dT / vtx[ivx].TimeErr;
    totChi += chi * chi;
    
    return totChi;
    
  } // PointVertexChi
  
    
} // namespace cluster
