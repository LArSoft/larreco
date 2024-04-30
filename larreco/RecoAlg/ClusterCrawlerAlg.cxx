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
#include <algorithm> // std::fill(), std::find(), std::sort()...
#include <cmath>
#include <iomanip>
#include <iostream>

// framework libraries
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larreco/RecoAlg/ClusterCrawlerAlg.h"

namespace {
  struct CluLen {
    int index;
    int length;
  };

  bool greaterThan(CluLen c1, CluLen c2)
  {
    return c1.length > c2.length;
  }
}

namespace cluster {

  //------------------------------------------------------------------------------
  ClusterCrawlerAlg::ClusterCrawlerAlg(fhicl::ParameterSet const& pset)
  {
    fNumPass = pset.get<unsigned short>("NumPass", 0);
    fMaxHitsFit = pset.get<std::vector<unsigned short>>("MaxHitsFit");
    fMinHits = pset.get<std::vector<unsigned short>>("MinHits");
    fNHitsAve = pset.get<std::vector<unsigned short>>("NHitsAve");
    fChgCut = pset.get<std::vector<float>>("ChgCut");
    fChiCut = pset.get<std::vector<float>>("ChiCut");
    fMaxWirSkip = pset.get<std::vector<unsigned short>>("MaxWirSkip");
    fMinWirAfterSkip = pset.get<std::vector<unsigned short>>("MinWirAfterSkip");
    fKinkChiRat = pset.get<std::vector<float>>("KinkChiRat");
    fKinkAngCut = pset.get<std::vector<float>>("KinkAngCut");
    fDoMerge = pset.get<std::vector<bool>>("DoMerge");
    fTimeDelta = pset.get<std::vector<float>>("TimeDelta");
    fMergeChgCut = pset.get<std::vector<float>>("MergeChgCut");
    fFindVertices = pset.get<std::vector<bool>>("FindVertices");
    fLACrawl = pset.get<std::vector<bool>>("LACrawl");
    fMinAmp = pset.get<std::vector<float>>("MinAmp", {5, 5, 5});
    fChgNearWindow = pset.get<float>("ChgNearWindow");
    fChgNearCut = pset.get<float>("ChgNearCut");

    fChkClusterDS = pset.get<bool>("ChkClusterDS", false);
    fVtxClusterSplit = pset.get<bool>("VtxClusterSplit", false);
    fFindStarVertices = pset.get<bool>("FindStarVertices", false);
    if (pset.has_key("HammerCluster")) {
      mf::LogWarning("CC")
        << "fcl setting HammerCluster is replaced by FindHammerClusters. Ignoring...";
    }
    fFindHammerClusters = pset.get<bool>("FindHammerClusters", false);
    fKillGarbageClusters = pset.get<float>("KillGarbageClusters", 0);
    fRefineVertexClusters = pset.get<bool>("RefineVertexClusters", false);
    fHitErrFac = pset.get<float>("HitErrFac", 0.2);
    fHitMinAmp = pset.get<float>("HitMinAmp", 0.2);
    fClProjErrFac = pset.get<float>("ClProjErrFac", 4);
    fMinHitFrac = pset.get<float>("MinHitFrac", 0.6);

    fLAClusAngleCut = pset.get<float>("LAClusAngleCut", 45);
    fLAClusMaxHitsFit = pset.get<unsigned short>("LAClusMaxHitsFit");
    fMergeAllHits = pset.get<bool>("MergeAllHits", false);
    fHitMergeChiCut = pset.get<float>("HitMergeChiCut", 2.5);
    fMergeOverlapAngCut = pset.get<float>("MergeOverlapAngCut");
    fAllowNoHitWire = pset.get<unsigned short>("AllowNoHitWire", 0);
    fVertex2DCut = pset.get<float>("Vertex2DCut", 5);
    fVertex2DWireErrCut = pset.get<float>("Vertex2DWireErrCut", 5);
    fVertex3DCut = pset.get<float>("Vertex3DCut", 5);

    fDebugPlane = pset.get<int>("DebugPlane", -1);
    fDebugWire = pset.get<int>("DebugWire", -1);
    fDebugHit = pset.get<int>("DebugHit", -1);

    // some error checking
    bool badinput = false;
    if (fNumPass > fMaxHitsFit.size()) badinput = true;
    if (fNumPass > fMinHits.size()) badinput = true;
    if (fNumPass > fNHitsAve.size()) badinput = true;
    if (fNumPass > fChgCut.size()) badinput = true;
    if (fNumPass > fChiCut.size()) badinput = true;
    if (fNumPass > fMaxWirSkip.size()) badinput = true;
    if (fNumPass > fMinWirAfterSkip.size()) badinput = true;
    if (fNumPass > fKinkChiRat.size()) badinput = true;
    if (fNumPass > fKinkAngCut.size()) badinput = true;
    if (fNumPass > fDoMerge.size()) badinput = true;
    if (fNumPass > fTimeDelta.size()) badinput = true;
    if (fNumPass > fMergeChgCut.size()) badinput = true;
    if (fNumPass > fFindVertices.size()) badinput = true;
    if (fNumPass > fLACrawl.size()) badinput = true;

    if (badinput)
      throw art::Exception(art::errors::Configuration)
        << "ClusterCrawlerAlg: Bad input from fcl file";

  } // reconfigure

  // used for sorting hits on wires
  bool SortByLowHit(unsigned int i, unsigned int j)
  {
    return i > j;
  }

  bool ClusterCrawlerAlg::SortByMultiplet(recob::Hit const& a, recob::Hit const& b)
  {
    // compare the wire IDs first:
    int cmp_res = a.WireID().cmp(b.WireID());
    if (cmp_res != 0) return cmp_res < 0; // order is decided, unless equal
    // decide by start time
    if (a.StartTick() != b.StartTick()) return a.StartTick() < b.StartTick();
    // if still undecided, resolve by local index
    return a.LocalIndex() < b.LocalIndex(); // if still unresolved, it's a bug!
  }                                         // ClusterCrawlerAlg::SortByMultiplet()

  //------------------------------------------------------------------------------
  void ClusterCrawlerAlg::ClearResults()
  {
    fHits.clear();
    tcl.clear();
    vtx.clear();
    vtx3.clear();
    inClus.clear();
  } // ClusterCrawlerAlg::ClearResults()

  //------------------------------------------------------------------------------
  void ClusterCrawlerAlg::CrawlInit()
  {
    prt = false;
    vtxprt = false;
    NClusters = 0;
    clBeginSlp = 0;
    clBeginSlpErr = 0;
    clBeginTim = 0;
    clBeginWir = 0;
    clBeginChg = 0;
    clBeginChgNear = 0;
    clEndSlp = 0;
    clEndSlpErr = 0;
    clEndTim = 0;
    clEndWir = 0;
    clEndChg = 0;
    clEndChgNear = 0;
    clChisq = 0;
    clStopCode = 0;
    clProcCode = 0;
    fFirstWire = 0;
    fLastWire = 0;
    fAveChg = 0.;
    fChgSlp = 0.;
    pass = 0;
    fScaleF = 0;
    WireHitRange.clear();

    ClearResults();
  }

  //------------------------------------------------------------------------------
  void ClusterCrawlerAlg::ClusterInit()
  {
    fcl2hits.clear();
    chifits.clear();
    hitNear.clear();
    chgNear.clear();
    fAveChg = -1.;
    fAveHitWidth = -1;
    clEndChg = -1.;
    clStopCode = 0;
    clProcCode = pass;
  }

  //------------------------------------------------------------------------------
  void ClusterCrawlerAlg::RunCrawler(detinfo::DetectorClocksData const& clock_data,
                                     detinfo::DetectorPropertiesData const& det_prop,
                                     std::vector<recob::Hit> const& srchits)
  {
    // Run the ClusterCrawler algorithm - creating seed clusters and crawling upstream.

    CrawlInit();

    fHits = srchits; // plain copy of the sources; it's the base of our hit result

    if (fHits.size() < 3) return;
    if (fHits.size() > UINT_MAX) {
      mf::LogWarning("CC") << "Too many hits for ClusterCrawler " << fHits.size();
      return;
    }

    // don't do anything...
    if (fNumPass == 0) return;

    // sort it as needed;
    // that is, sorted by wire ID number,
    // then by start of the region of interest in time, then by the multiplet
    std::sort(fHits.begin(), fHits.end(), &SortByMultiplet);

    inClus.resize(fHits.size());
    mergeAvailable.resize(fHits.size());
    for (unsigned int iht = 0; iht < inClus.size(); ++iht) {
      inClus[iht] = 0;
      mergeAvailable[iht] = false;
    }

    // FIXME (KJK): The 'cstat', 'tpc', and 'plane' class variables should be removed.
    for (geo::TPCID const& tpcid : geom->Iterate<geo::TPCID>()) {
      cstat = tpcid.Cryostat;
      tpc = tpcid.TPC;
      for (geo::PlaneID const& planeid : geom->Iterate<geo::PlaneID>(tpcid)) {
        plane = planeid.Plane;
        WireHitRange.clear();
        // define a code to ensure clusters are compared within the same plane
        clCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, planeid.Plane);
        cstat = tpcid.Cryostat;
        tpc = tpcid.TPC;
        // fill the WireHitRange vector with first/last hit on each wire
        // dead wires and wires with no hits are flagged < 0
        GetHitRange(clCTP);

        // sanity check
        if (WireHitRange.empty() || (fFirstWire == fLastWire)) continue;
        raw::ChannelID_t channel = fHits[fFirstHit].Channel();
        // get the scale factor to convert dTick/dWire to dX/dU. This is used
        // to make the kink and merging cuts
        float wirePitch = geom->WirePitch(geom->View(channel));
        float tickToDist = det_prop.DriftVelocity(det_prop.Efield(), det_prop.Temperature());
        tickToDist *= 1.e-3 * sampling_rate(clock_data); // 1e-3 is conversion of 1/us to 1/ns
        fScaleF = tickToDist / wirePitch;
        // convert Large Angle Cluster crawling cut to a slope cut
        if (fLAClusAngleCut > 0)
          fLAClusSlopeCut = std::tan(3.142 * fLAClusAngleCut / 180.) / fScaleF;
        fMaxTime = det_prop.NumberTimeSamples();
        fNumWires = geom->Nwires(planeid);
        // look for clusters
        if (fNumPass > 0) ClusterLoop();
      } // plane
      if (fVertex3DCut > 0) {
        // Match vertices in 3 planes
        VtxMatch(det_prop, tpcid);
        Vtx3ClusterMatch(det_prop, tpcid);
        if (fFindHammerClusters) FindHammerClusters(det_prop);
        // split clusters using 3D vertices
        Vtx3ClusterSplit(det_prop, tpcid);
      }
      if (fDebugPlane >= 0) {
        mf::LogVerbatim("CC") << "Clustering done in TPC ";
        PrintClusters();
      }
    } // for all tpcs

    // clean up
    WireHitRange.clear();
    fcl2hits.clear();
    chifits.clear();
    hitNear.clear();
    chgNear.clear();

    // remove the hits that have become obsolete
    RemoveObsoleteHits();

  } // RunCrawler

  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::ClusterLoop()
  {
    // looks for seed clusters in a plane and crawls along a trail of hits

    unsigned int nHitsUsed = 0, iwire, jwire, kwire;
    bool AllDone = false, SigOK = false, HitOK = false;
    unsigned int ihit, jhit;
    for (unsigned short thispass = 0; thispass < fNumPass; ++thispass) {
      pass = thispass;
      // look for a starting cluster that spans a block of wires
      unsigned int span = 3;
      if (fMinHits[pass] < span) span = fMinHits[pass];
      for (iwire = fLastWire; iwire > fFirstWire + span; --iwire) {
        // skip bad wires or no hits on the wire
        if (WireHitRange[iwire].first < 0) continue;
        auto ifirsthit = (unsigned int)WireHitRange[iwire].first;
        auto ilasthit = (unsigned int)WireHitRange[iwire].second;
        for (ihit = ifirsthit; ihit < ilasthit; ++ihit) {
          bool ClusterAdded = false;
          recob::Hit const& hit = fHits[ihit];
          // skip used hits
          if (ihit > fHits.size() - 1) {
            mf::LogError("CC") << "ClusterLoop bad ihit " << ihit << " fHits size " << fHits.size();
            return;
          }
          // skip used and obsolete hits
          if (inClus[ihit] != 0) continue;
          // skip narrow hits
          if (fHits[ihit].PeakAmplitude() < fHitMinAmp) continue;
          if ((iwire + 1) < span) continue;
          jwire = iwire - span + 1;
          if (prt)
            mf::LogVerbatim("CC") << "Found debug hit " << PrintHit(ihit) << " on pass" << pass;
          // skip if good wire and no hit
          if (WireHitRange[jwire].first == -2) continue;
          if (WireHitRange[jwire].first == -1) {
            // Found a dead jwire. Keep looking upstream until we find a good wire
            unsigned int nmissed = 0;
            while (WireHitRange[jwire].first == -1 && jwire > 1 && nmissed < fMaxWirSkip[pass]) {
              --jwire;
              ++nmissed;
            }
            if (prt)
              mf::LogVerbatim("CC")
                << " new jwire " << jwire << " dead? " << WireHitRange[jwire].first;
            if (WireHitRange[jwire].first < 0) continue;
          } // dead jwire
          // Find the hit on wire jwire that best matches a line between
          // a nearby vertex and hit ihit. No constraint if useHit < 0
          unsigned int useHit = 0;
          bool doConstrain = false;
          VtxConstraint(iwire, ihit, jwire, useHit, doConstrain);
          unsigned int jfirsthit = (unsigned int)WireHitRange[jwire].first;
          unsigned int jlasthit = (unsigned int)WireHitRange[jwire].second;
          if (jfirsthit > fHits.size() - 1 || jfirsthit > fHits.size() - 1)
            throw art::Exception(art::errors::LogicError)
              << "ClusterLoop jwire " << jwire << " bad firsthit " << jfirsthit << " lasthit "
              << jlasthit << " fhits size " << fHits.size();
          for (jhit = jfirsthit; jhit < jlasthit; ++jhit) {
            if (jhit > fHits.size() - 1)
              throw art::Exception(art::errors::LogicError)
                << "ClusterLoop bad jhit " << jhit << " firsthit " << jfirsthit << " lasthit "
                << jlasthit << " fhits size" << fHits.size();
            // Constraint?
            if (doConstrain && jhit != useHit) continue;
            recob::Hit const& other_hit = fHits[jhit];
            // skip used and obsolete hits
            if (inClus[jhit] != 0) continue;
            // skip narrow hits
            if (fHits[jhit].PeakAmplitude() < fHitMinAmp) continue;
            // start a cluster with these two hits
            ClusterInit();
            fcl2hits.push_back(ihit);
            chifits.push_back(0.);
            hitNear.push_back(0);
            chgNear.push_back(0); // These will be defined if the cluster survives the cuts
            // enter the jhit
            fcl2hits.push_back(jhit);
            chifits.push_back(0.);
            hitNear.push_back(0);
            chgNear.push_back(0);
            clLA = false;
            clpar[0] = other_hit.PeakTime();
            clpar[1] = (hit.PeakTime() - other_hit.PeakTime()) / (iwire - jwire);
            // increase slope errors for large angle clusters
            clparerr[1] = 0.2 * std::abs(clpar[1]);
            clpar[2] = fHits[jhit].WireID().Wire;
            clChisq = 0;
            // now look for hits to add on the intervening wires
            bool clok = false;
            for (kwire = jwire + 1; kwire < iwire; ++kwire) {
              // ensure this cluster doesn't cross a vertex
              if (CrawlVtxChk(kwire)) {
                clok = false;
                break;
              }
              AddHit(kwire, HitOK, SigOK);
              if (prt)
                mf::LogVerbatim("CC") << " HitOK " << HitOK << " clChisq " << clChisq << " cut "
                                      << fChiCut[pass] << " ClusterHitsOK " << ClusterHitsOK(-1);
              // No hit found
              if (!HitOK) break;
              // This should be a really good chisq
              if (clChisq > 2) break;
              // hit widths & overlap not consistent
              if (!ClusterHitsOK(-1)) continue;
              clok = true;
            }
            // drop it?
            if (!clok) continue;
            prt = (fDebugPlane == (int)plane && (int)iwire == fDebugWire &&
                   std::abs((int)hit.PeakTime() - fDebugHit) < 20);
            if (prt)
              mf::LogVerbatim("CC")
                << "ADD >>>>> Starting cluster with hits " << PrintHit(fcl2hits[0]) << " "
                << PrintHit(fcl2hits[1]) << " " << PrintHit(fcl2hits[2]) << " on pass " << pass;
            // save the cluster begin info
            clBeginWir = iwire;
            clBeginTim = hit.PeakTime();
            clBeginSlp = clpar[1];
            // don't do a small angle crawl if the cluster slope is too large
            // and Large Angle crawling is NOT requested on this pass
            if (!fLACrawl[pass] && std::abs(clBeginSlp) > fLAClusSlopeCut) continue;
            // See if we are trying to start a cluster between a vertex
            // and a cluster that is associated to that vertex. If so, skip it
            if (CrawlVtxChk2()) continue;
            clBeginSlpErr = clparerr[1];
            clBeginChg = 0;
            // Calculate the average width
            fAveHitWidth = 0;
            for (unsigned short kk = 0; kk < fcl2hits.size(); ++kk) {
              fAveHitWidth += fHits[fcl2hits[kk]].EndTick() - fHits[fcl2hits[kk]].StartTick();
            }
            fAveHitWidth /= (float)fcl2hits.size();
            // decide whether to crawl a large angle cluster. Requirements are:
            // 1) the user has set the LACluster angle cut > 0, AND
            // 2) the cluster slope exceeds the cut
            // Note that if condition 1 is met, normal cluster crawling is done
            // only if the slope is less than the cut
            if (fLACrawl[pass] && fLAClusSlopeCut > 0) {
              // LA cluster crawling requested
              if (std::abs(clBeginSlp) > fLAClusSlopeCut) { LACrawlUS(); }
              else {
                CrawlUS();
              } // std::abs(clBeginSlp) > fLAClusAngleCut
            }
            else {
              // allow clusters of any angle
              CrawlUS();
            } // fLAClusSlopeCut > 0
            if (fcl2hits.size() >= fMinHits[pass]) {
              // it's long enough so save it
              clEndSlp = clpar[1]; // save the slope at the end
              clEndSlpErr = clparerr[1];
              // store the cluster
              if (!TmpStore()) {
                mf::LogError("CC") << "Failed to store cluster in plane " << plane;
                continue;
              }
              ClusterAdded = true;
              nHitsUsed += fcl2hits.size();
              AllDone = (nHitsUsed == fHits.size());
              break;
            }
            else {
              // abandon it
              if (prt) mf::LogVerbatim("CC") << "ClusterLoop: dropped the cluster";
            }
            if (ClusterAdded || AllDone) break;
          } // jhit
          if (AllDone) break;
        } // ihit
        if (AllDone) break;
      } // iwire

      // try to merge clusters
      if (fDoMerge[pass]) ChkMerge();
      // form 2D vertices
      if (fFindVertices[pass]) FindVertices();

      if (AllDone) break;

    } // pass

    // Kill Garbage clusters
    if (fKillGarbageClusters > 0 && !tcl.empty()) KillGarbageClusters();
    // Merge overlapping clusters
    if (fMergeOverlapAngCut > 0) MergeOverlap();
    // Check the DS end of clusters
    if (fChkClusterDS) ChkClusterDS();
    // split clusters using vertices
    if (fVtxClusterSplit) {
      bool didSomething = VtxClusterSplit();
      // iterate once to handle the case where a cluster crosses two vertices
      if (didSomething) VtxClusterSplit();
    }
    // Look for 2D vertices with star topology - short, back-to-back clusters
    if (fFindStarVertices) FindStarVertices();

    if (fDebugPlane == (int)plane) {
      mf::LogVerbatim("CC") << "Clustering done in plane " << plane;
      PrintClusters();
    }

    CheckHitClusterAssociations();

  } // ClusterLoop()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::KillGarbageClusters()
  {
    // Ghost Clusters:

    if (tcl.size() < 2) return;

    unsigned short icl, jcl;
    // This code preferentially selects icl clusters that were
    // reconstructed on an early pass (unless they were split)
    std::vector<float> iHits, jHits;
    unsigned int indx;
    // find the average hit width on the first pass and construct
    // a hit separation cut
    float sepcut = 0, iFrac, jFrac;
    bool first = true;
    unsigned short iLoIndx, jLoIndx, olapSize, iop, ii, jj;
    unsigned short nclose;
    float iChg, jChg;
    // vecrtor of clusters that will be killed after all is done
    std::vector<unsigned short> killMe;
    bool doKill;
    for (icl = 0; icl < tcl.size() - 1; ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != clCTP) continue;
      // put the hits into a wire ordered vector
      iHits.clear();
      // initialize to a large positive value
      iHits.resize(tcl[icl].BeginWir - tcl[icl].EndWir + 1, 1000);
      iChg = 0;
      for (auto iht : tcl[icl].tclhits) {
        indx = fHits[iht].WireID().Wire - tcl[icl].EndWir;
        if (indx > iHits.size() - 1) {
          mf::LogWarning("CC") << "KillGarbageClusters: icl ID " << tcl[icl].ID << " Bad indx "
                               << indx << " " << iHits.size() << "\n";
          continue;
        }
        iHits[indx] = fHits[iht].PeakTime();
        iChg += fHits[iht].Integral();
        if (first) sepcut += fHits[iht].RMS();
      } // iht
      if (first) {
        sepcut /= (float)tcl[icl].tclhits.size();
        // clusters are consider ghost candidates if many hits
        // are within sepcut of each other on the same wire
        sepcut *= 10;
        first = false;
      } // first
      for (jcl = icl + 1; jcl < tcl.size(); ++jcl) {
        if (tcl[jcl].ID < 0) continue;
        if (tcl[jcl].CTP != clCTP) continue;
        // ignore if there is no overlap
        if (tcl[icl].BeginWir < tcl[jcl].EndWir) continue;
        if (tcl[icl].EndWir > tcl[jcl].BeginWir) continue;
        // require similar angle
        if (std::abs(tcl[icl].BeginAng - tcl[jcl].BeginAng) > fKillGarbageClusters) continue;
        // find the overlap region
        if (tcl[icl].EndWir < tcl[jcl].EndWir) {
          //  icl E-----------....
          //  jcl    E----------....
          // olap    xxxxxxxxxx...
          iLoIndx = tcl[jcl].EndWir - tcl[icl].EndWir;
          jLoIndx = 0;
          if (tcl[icl].BeginWir < tcl[jcl].BeginWir) {
            //  icl E-----------B
            //  jcl    E------------B
            // olap    xxxxxxxxxx
            olapSize = tcl[icl].BeginWir - tcl[jcl].EndWir + 1;
          }
          else {
            //  icl E-----------B
            //  jcl    E-----B
            // olap    xxxxxxx
            olapSize = tcl[jcl].BeginWir - tcl[jcl].EndWir + 1;
          } // iBegin
        }   // iEnd < jEnd
        else {
          //  icl    E-----------....
          //  jcl E----------....
          // olap    xxxxxxxxxx...
          iLoIndx = 0;
          jLoIndx = tcl[icl].EndWir - tcl[icl].EndWir;
          if (tcl[icl].BeginWir < tcl[jcl].BeginWir) {
            //  icl    E-----B
            //  jcl E-----------B
            // olap    xxxxxxx
            olapSize = tcl[icl].BeginWir - tcl[icl].EndWir + 1;
          }
          else {
            //  icl    E-----------B
            //  jcl E----------B
            // olap    xxxxxxxxx
            olapSize = tcl[jcl].BeginWir - tcl[icl].EndWir + 1;
          }
        } // iEnd > jEnd
        jHits.clear();
        // initialize to a large negative value
        jHits.resize(tcl[jcl].BeginWir - tcl[jcl].EndWir + 1, -1000);
        jChg = 0;
        for (auto jht : tcl[jcl].tclhits) {
          indx = fHits[jht].WireID().Wire - tcl[jcl].EndWir;
          if (indx > jHits.size() - 1) {
            mf::LogWarning("CC") << "KillGarbageClusters: jcl ID " << tcl[jcl].ID << " Bad indx "
                                 << indx << " " << jHits.size() << "\n";
            continue;
          }
          jHits[indx] = fHits[jht].PeakTime();
          jChg += fHits[jht].Integral();
        } // jht
        // count the number of close hits
        nclose = 0;
        for (iop = 0; iop < olapSize; ++iop) {
          ii = iLoIndx + iop;
          if (ii > iHits.size() - 1) continue;
          jj = jLoIndx + iop;
          if (jj > jHits.size() - 1) continue;
          if (std::abs(iHits[ii] - jHits[jj]) < sepcut) ++nclose;
        } // iop
        iFrac = (float)nclose / (float)iHits.size();
        jFrac = (float)nclose / (float)jHits.size();
        if (iFrac < 0.5 && jFrac < 0.5) continue;
        doKill = (iFrac < jFrac && iChg > jChg);
        if (doKill) killMe.push_back(jcl);
      } // jcl
    }   // icl

    if (killMe.size() == 0) return;
    for (auto icl : killMe) {
      // killing time
      if (tcl[icl].ID < 0) continue;
      tcl[icl].ProcCode = 666;
      MakeClusterObsolete(icl);
    } // icl

  } // KillGarbageClusters

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

    bool chkprt = (fDebugWire == 666);
    if (chkprt) mf::LogVerbatim("CC") << "Inside MergeOverlap using clCTP " << clCTP;

    unsigned short minLen = 6;
    unsigned short minOvrLap = 2;

    // use the loosest dTick cut which is probably the last one
    float maxDTick = fTimeDelta[fTimeDelta.size() - 1];

    unsigned int overlapSize, ii, indx, bWire, eWire;
    unsigned int iht, jht;
    float dang, prtime, dTick;
    for (icl = 0; icl < tcl.size(); ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != clCTP) continue;
      prt = chkprt && fDebugPlane == (int)clCTP;
      if (tcl[icl].BeginVtx >= 0) continue;
      if (tcl[icl].tclhits.size() < minLen) continue;
      for (jcl = 0; jcl < tcl.size(); ++jcl) {
        if (icl == jcl) continue;
        if (tcl[jcl].ID < 0) continue;
        if (tcl[jcl].CTP != clCTP) continue;
        if (tcl[jcl].EndVtx >= 0) continue;
        if (tcl[jcl].tclhits.size() < minLen) continue;
        // icl Begin is not far enough DS from the end of jcl
        if (tcl[icl].BeginWir < tcl[jcl].EndWir + minOvrLap) continue;
        // and it doesn't end within the wire boundaries of jcl
        if (tcl[icl].BeginWir > tcl[jcl].BeginWir - minOvrLap) continue;
        // jcl End isn't far enough US from the end of icl
        if (tcl[jcl].EndWir < tcl[icl].EndWir + minOvrLap) continue;
        dang = std::abs(tcl[icl].BeginAng - tcl[jcl].EndAng);
        if (prt)
          mf::LogVerbatim("CC") << "MergeOverlap icl ID " << tcl[icl].ID << " jcl ID "
                                << tcl[jcl].ID << " dang " << dang;
        if (dang > 0.5) continue;
        overlapSize = tcl[icl].BeginWir - tcl[jcl].EndWir + 1;
        eWire = tcl[jcl].EndWir;
        bWire = tcl[icl].BeginWir;
        if (prt)
          mf::LogVerbatim("CC") << " Candidate icl ID " << tcl[icl].ID << " " << tcl[icl].EndWir
                                << "-" << tcl[icl].BeginWir << " jcl ID " << tcl[jcl].ID << " "
                                << tcl[jcl].EndWir << "-" << tcl[jcl].BeginWir << " overlapSize "
                                << overlapSize << " bWire " << bWire << " eWire " << eWire;
        iht = 0;
        jht = 0;
        for (ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
          iht = tcl[icl].tclhits[ii];
          if (fHits[iht].WireID().Wire < eWire) break;
        } // ii
        // require that the ends be similar in time
        dTick = std::abs(fHits[iht].PeakTime() - tcl[jcl].EndTim);
        if (dTick > maxDTick) continue;
        if (prt)
          mf::LogVerbatim("CC") << " dTick icl iht time " << PrintHit(iht) << " jcl EndTim "
                                << tcl[jcl].EndTim << " dTick " << dTick;
        for (ii = 0; ii < tcl[jcl].tclhits.size(); ++ii) {
          jht = tcl[jcl].tclhits[tcl[jcl].tclhits.size() - ii - 1];
          if (fHits[jht].WireID().Wire > bWire) break;
        } // ii
        dTick = std::abs(fHits[jht].PeakTime() - tcl[icl].BeginTim);
        if (dTick > maxDTick) continue;
        if (prt)
          mf::LogVerbatim("CC") << " dTick jcl jht time " << PrintHit(jht) << " icl BeginTim "
                                << tcl[icl].BeginTim << " dTick " << dTick;
        // Calculate the line between iht and jht
        clpar[0] = fHits[iht].PeakTime();
        clpar[2] = fHits[iht].WireID().Wire;
        clpar[1] = (fHits[jht].PeakTime() - fHits[iht].PeakTime()) /
                   ((float)fHits[jht].WireID().Wire - clpar[2]);
        // put the hits in the overlap region into a vector if they are close to the line
        std::vector<unsigned int> oWireHits(overlapSize, INT_MAX);
        std::vector<float> delta(overlapSize, maxDTick);
        for (ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
          iht = tcl[icl].tclhits[ii];
          if (fHits[iht].WireID().Wire < eWire) break;
          prtime = clpar[0] + clpar[1] * ((float)fHits[iht].WireID().Wire - clpar[2]);
          dTick = std::abs(fHits[iht].PeakTime() - prtime);
          indx = fHits[iht].WireID().Wire - eWire;
          if (dTick > delta[indx]) continue;
          delta[indx] = dTick;
          oWireHits[indx] = iht;
        } // ii
        // enter the second set of hits
        for (ii = 0; ii < tcl[jcl].tclhits.size(); ++ii) {
          jht = tcl[jcl].tclhits[tcl[jcl].tclhits.size() - ii - 1];
          if (fHits[jht].WireID().Wire > bWire) break;
          prtime = clpar[0] + clpar[1] * ((float)fHits[jht].WireID().Wire - clpar[2]);
          dTick = std::abs(fHits[jht].PeakTime() - prtime);
          indx = fHits[jht].WireID().Wire - eWire;
          if (dTick > delta[indx]) continue;
          delta[indx] = dTick;
          oWireHits[indx] = jht;
        } // ii
        // stuff them into fcl2hits
        fcl2hits.clear();
        for (ii = 0; ii < oWireHits.size(); ++ii) {
          if (oWireHits[ii] == INT_MAX) continue;
          iht = oWireHits[ii];
          fcl2hits.push_back(iht);
          if (prt) mf::LogVerbatim("CC") << "hit " << PrintHit(iht);
        } // ii
        if (fcl2hits.size() < 0.5 * overlapSize) continue;
        if (fcl2hits.size() < 3) continue;
        std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
        FitCluster();
        if (prt)
          mf::LogVerbatim("CC") << " Overlap size " << overlapSize << " fit chisq " << clChisq
                                << " nhits " << fcl2hits.size();
        if (clChisq > 20) continue;
        // save these hits so we can paste them back on fcl2hits when merging
        std::vector<unsigned int> oHits = fcl2hits;
        // prepare to make a new cluster
        TmpGet(jcl);
        // resize it
        unsigned short jclNewSize;
        for (jclNewSize = 0; jclNewSize < fcl2hits.size(); ++jclNewSize) {
          iht = fcl2hits[jclNewSize];
          if (fHits[iht].WireID().Wire <= bWire) break;
        } // jclNewSize
        if (prt) {
          mf::LogVerbatim("CC") << "jcl old size " << fcl2hits.size() << " newSize " << jclNewSize;
          iht = fcl2hits[fcl2hits.size() - 1];
          unsigned int iiht = fcl2hits[jclNewSize - 1];
          mf::LogVerbatim("CC") << "jcl old last wire " << fHits[iht].WireID().Wire
                                << " After resize last wire " << fHits[iiht].WireID().Wire;
        }
        fcl2hits.resize(jclNewSize);
        // append the hits in the overlap region
        fcl2hits.insert(fcl2hits.end(), oHits.begin(), oHits.end());
        // now paste in the icl hits that are US of the overlap region
        for (ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
          iht = tcl[icl].tclhits[ii];
          if ((unsigned int)fHits[iht].WireID().Wire >= eWire) continue;
          fcl2hits.insert(fcl2hits.end(), tcl[icl].tclhits.begin() + ii, tcl[icl].tclhits.end());
          break;
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
        MakeClusterObsolete(icl);
        MakeClusterObsolete(jcl);
        if (!TmpStore()) {
          // Merged cluster is fubar. Try to recover
          RestoreObsoleteCluster(icl);
          RestoreObsoleteCluster(jcl);
          CheckHitClusterAssociations();
          continue;
        }
        CheckHitClusterAssociations();
        tcl[tcl.size() - 1].BeginVtx = tcl[jcl].BeginVtx;
        tcl[tcl.size() - 1].EndVtx = tcl[icl].BeginVtx;
        // after this point any failure should result in a jcl loop break
        if (prt)
          mf::LogVerbatim("CC") << "MergeOverlap new long cluster ID " << tcl[tcl.size() - 1].ID
                                << " in clCTP " << clCTP;
        // icl cluster made obsolete so we must have done something
        if (tcl[icl].ID < 0) break;
      } // jcl
    }   // icl

  } // MergeOverlap()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::MakeClusterObsolete(unsigned short icl)
  {
    short& ID = tcl[icl].ID;
    if (ID <= 0) {
      mf::LogError("CC") << "Trying to make already-obsolete cluster obsolete ID = " << ID;
      return; // already obsolete
    }
    ID = -ID; // mark the cluster as obsolete

    // release the hits
    for (unsigned int iht = 0; iht < tcl[icl].tclhits.size(); ++iht)
      inClus[tcl[icl].tclhits[iht]] = 0;

  } // ClusterCrawlerAlg::MakeClusterObsolete()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::RestoreObsoleteCluster(unsigned short icl)
  {
    short& ID = tcl[icl].ID;
    if (ID > 0) {
      mf::LogError("CC") << "Trying to restore non-obsolete cluster ID = " << ID;
      return;
    }
    ID = -ID;

    for (unsigned short iht = 0; iht < tcl[icl].tclhits.size(); ++iht)
      inClus[tcl[icl].tclhits[iht]] = ID;

  } // RestoreObsoleteCluster()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::FclTrimUS(unsigned short nTrim)
  {

    // Trims nTrim hits off the UpStream end of the fcl2hits, etc vectors.
    if (nTrim == 0) return;

    if (nTrim >= fcl2hits.size()) nTrim = fcl2hits.size();

    //    RestoreUnMergedClusterHits((short)nTrim);
    for (unsigned short ii = 0; ii < nTrim; ++ii) {
      fcl2hits.pop_back();
      chifits.pop_back();
      hitNear.pop_back();
      chgNear.pop_back();
    } // ii

  } // FclTrim

  //////////////////////////////////////////
  void ClusterCrawlerAlg::CheckHitClusterAssociations()
  {
    // check hit - cluster associations

    if (fHits.size() != inClus.size()) {
      mf::LogError("CC") << "CHCA: Sizes wrong " << fHits.size() << " " << inClus.size();
      return;
    }

    unsigned int iht, nErr = 0;
    short clID;

    // check cluster -> hit association
    for (unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if (tcl[icl].ID < 0) continue;
      clID = tcl[icl].ID;
      for (unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
        iht = tcl[icl].tclhits[ii];
        if (iht > fHits.size() - 1) {
          mf::LogError("CC") << "CHCA: Bad tclhits index " << iht << " fHits size " << fHits.size();
          return;
        } // iht > fHits.size() - 1
        if (inClus[iht] != clID) {
          mf::LogError("CC") << "CHCA: Bad cluster -> hit association. clID " << clID
                             << " hit inClus " << inClus[iht] << " ProcCode " << tcl[icl].ProcCode
                             << " CTP " << tcl[icl].CTP;
          ++nErr;
          if (nErr > 10) return;
        }
      } // ii
    }   // icl

    // check hit -> cluster association
    unsigned short icl;
    for (iht = 0; iht < fHits.size(); ++iht) {
      if (inClus[iht] <= 0) continue;
      icl = inClus[iht] - 1;
      // see if the cluster is obsolete
      if (tcl[icl].ID < 0) {
        mf::LogError("CC") << "CHCA: Hit associated with an obsolete cluster. hit W:T "
                           << fHits[iht].WireID().Wire << ":" << (int)fHits[iht].PeakTime()
                           << " tcl[icl].ID " << tcl[icl].ID;
        ++nErr;
        if (nErr > 10) return;
      }
      if (std::find(tcl[icl].tclhits.begin(), tcl[icl].tclhits.end(), iht) ==
          tcl[icl].tclhits.end()) {
        mf::LogError("CC") << "CHCA: Bad hit -> cluster association. hit index " << iht << " W:T "
                           << fHits[iht].WireID().Wire << ":" << (int)fHits[iht].PeakTime()
                           << " inClus " << inClus[iht];
        ++nErr;
        if (nErr > 10) return;
      }
    } // iht

  } // CheckHitClusterAssociations()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::RemoveObsoleteHits()
  {

    unsigned int destHit = 0;

    if (fHits.size() != inClus.size()) {
      mf::LogError("CC") << "RemoveObsoleteHits size mis-match " << fHits.size() << " "
                         << inClus.size();
      return;
    }

    unsigned short icl;
    for (unsigned int srcHit = 0; srcHit < fHits.size(); ++srcHit) {
      if (inClus[srcHit] < 0) continue;
      if (srcHit != destHit) {
        fHits[destHit] = std::move(fHits[srcHit]);
        inClus[destHit] = inClus[srcHit];
        if (inClus[destHit] > 0) {
          // hit is in a cluster. Find it and change the index
          icl = inClus[destHit] - 1;
          auto& hits = tcl[icl].tclhits;
          auto iHitIndex = std::find(hits.begin(), hits.end(), srcHit);
          if (iHitIndex == hits.end()) {
            mf::LogError("CC") << "RemoveObsoleteHits: Hit #" << srcHit
                               << " not found in cluster ID " << inClus[destHit];
          }
          else {
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
  void ClusterCrawlerAlg::ChkClusterDS()
  {
    // Try to extend clusters DS by a few wires.
    // DS hits may not have been  included in a cluster if they have high
    // multiplicity or high charge.
    // Ref ClusterLoop cuts for starting a seed cluster.

    prt = (fDebugPlane == 3);

    // use the most generous kink angle cut
    float dThCut = fKinkAngCut[fNumPass - 1];

    if (prt)
      mf::LogVerbatim("CC") << "ChkClusterDS clCTP " << clCTP << " kink angle cut " << dThCut;

    const unsigned short tclsize = tcl.size();
    bool didMerge, skipit;
    unsigned short icl, ii, nhm, iv;
    unsigned int iht;

    // first merge any hits on the DS end of clusters
    for (icl = 0; icl < tclsize; ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != clCTP) continue;
      // ignore clusters that have a Begin vertex
      if (tcl[icl].BeginVtx >= 0) continue;
      // and clusters near a vertex
      skipit = false;
      for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if (vtx[iv].CTP != clCTP) continue;
        if (std::abs(vtx[iv].Wire - tcl[icl].BeginWir) > 4) continue;
        if (std::abs(vtx[iv].Time - tcl[icl].BeginTim) > 20) continue;
        skipit = true;
        break;
      }
      if (skipit) continue;
      // check the first few hits
      nhm = 0;
      for (ii = 0; ii < 3; ++ii) {
        iht = fcl2hits[ii];
        if (fHits[iht].Multiplicity() > 1) {
          MergeHits(iht, didMerge);
          if (didMerge) ++nhm;
        }
      } // ii
      if (nhm > 0) {
        // update the Begin parameters in-place
        FitClusterMid(icl, 0, 3);
        tcl[icl].BeginTim = clpar[0];
        tcl[icl].BeginSlp = clpar[1];
        tcl[icl].BeginAng = atan(fScaleF * clpar[1]);
        tcl[icl].BeginSlpErr = clparerr[1];
        tcl[icl].BeginChg = fAveChg;
        tcl[icl].ProcCode += 5000;
        if (prt) mf::LogVerbatim("CC") << "ChkClusterDS: Merge hits on cluster " << tcl[icl].ID;
      } // nhm > 0
    }   // icl

    float thhits, prevth, hitrms, rmsrat;
    bool ratOK;
    std::vector<unsigned int> dshits;
    for (icl = 0; icl < tclsize; ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != clCTP) continue;
      // ignore clusters that have a Begin vertex
      if (tcl[icl].BeginVtx >= 0) continue;
      // and clusters near a vertex
      skipit = false;
      for (iv = 0; iv < vtx.size(); ++iv) {
        if (vtx[iv].CTP != clCTP) continue;
        if (std::abs(vtx[iv].Wire - tcl[icl].BeginWir) > 4) continue;
        if (std::abs(vtx[iv].Time - tcl[icl].BeginTim) > 20) continue;
        skipit = true;
        break;
      }
      if (skipit) continue;
      // find the angle using the first 2 hits
      unsigned int ih0 = tcl[icl].tclhits[1];
      unsigned int ih1 = tcl[icl].tclhits[0];
      const float slp = (fHits[ih1].PeakTime() - fHits[ih0].PeakTime()) /
                        (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
      prevth = std::atan(fScaleF * slp);
      // move the "origin" to the first hit
      ih0 = ih1;
      unsigned int wire = fHits[ih0].WireID().Wire;
      hitrms = fHits[ih0].RMS();
      float time0 = fHits[ih0].PeakTime();
      float prtime;
      dshits.clear();
      // follow DS a few wires. Stop if any encountered
      // hit is associated with a cluster
      for (ii = 0; ii < 4; ++ii) {
        ++wire;
        if (wire > fLastWire) break;
        prtime = time0 + slp;
        if (prt)
          mf::LogVerbatim("CC") << "ChkClusterDS: Try to extend " << tcl[icl].ID << " to W:T "
                                << wire << " hitrms " << hitrms << " prevth " << prevth
                                << " prtime " << (int)prtime;
        // stop if no hits on this wire
        if (WireHitRange[wire].first == -2) break;
        unsigned int firsthit = WireHitRange[wire].first;
        unsigned int lasthit = WireHitRange[wire].second;
        bool hitAdded = false;
        for (ih1 = firsthit; ih1 < lasthit; ++ih1) {
          if (inClus[ih1] != 0) continue;
          if (prtime < fHits[ih1].PeakTimeMinusRMS(5)) continue;
          if (prtime > fHits[ih1].PeakTimePlusRMS(5)) continue;
          const float slp = (fHits[ih1].PeakTime() - fHits[ih0].PeakTime()) /
                            (fHits[ih1].WireID().Wire - fHits[ih0].WireID().Wire);
          thhits = std::atan(fScaleF * slp);
          if (prt)
            mf::LogVerbatim("CC") << " theta " << thhits << " prevth " << prevth << " cut "
                                  << dThCut;
          if (std::abs(thhits - prevth) > dThCut) continue;
          // make a hit rms cut for small angle clusters
          ratOK = true;
          if (std::abs(slp) < fLAClusSlopeCut) {
            rmsrat = fHits[ih1].RMS() / hitrms;
            ratOK = rmsrat > 0.3 && rmsrat < 3;
          }
          else {
            // merge the hits
            bool didMerge;
            MergeHits(ih1, didMerge);
          }
          if (prt) mf::LogVerbatim("CC") << " rmsrat " << rmsrat << " OK? " << ratOK;
          // require small angle and not wildly different width compared
          // to the first hit in the cluster
          // TODO handle hit multiplets here
          if (ratOK) {
            dshits.push_back(ih1);
            hitAdded = true;
            prevth = thhits;
            ih0 = ih1;
            if (prt)
              mf::LogVerbatim("CC") << " Add hit " << fHits[ih1].WireID().Wire << ":"
                                    << (int)fHits[ih1].PeakTime() << " rmsrat " << rmsrat;
            break;
          }
        } // ih1
        // stop looking if no hit was added on this wire
        if (!hitAdded) break;
      } // ii
      // Found hits not associated with a different cluster
      if (dshits.size() > 0) {
        // put the tcl cluster into the working vectors
        TmpGet(icl);
        // clobber the hits
        fcl2hits.clear();
        // sort the DS hits
        if (dshits.size() > 1) std::sort(dshits.begin(), dshits.end(), SortByLowHit);
        // stuff them into fcl2hits
        fcl2hits = dshits;
        // Append the existing hits
        for (ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
          // un-assign the hits so that TmpStore will re-assign them
          iht = tcl[icl].tclhits[ii];
          inClus[iht] = 0;
          fcl2hits.push_back(iht);
        }
        clProcCode += 5000;
        pass = fNumPass - 1;
        FitClusterChg();
        clBeginChg = fAveChg;
        // declare the old one obsolete
        MakeClusterObsolete(icl);
        // add the new one
        if (!TmpStore()) {
          mf::LogError("CC") << "ChkClusterDS TmpStore failed while extending cluster ID "
                             << tcl[icl].ID;
          continue;
        }
        const size_t newcl = tcl.size() - 1;
        if (prt) { mf::LogVerbatim("CC") << " Store " << newcl; }
        tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
        tcl[newcl].EndVtx = tcl[icl].EndVtx;
      } // dshits.size() > 0
    }   // icl
  }     // ChkClusterDS

  //////////////////////////////////////////
  bool ClusterCrawlerAlg::CrawlVtxChk2()
  {
    // returns true if the (presumably short) cluster under construction
    // resides between a vertex and another cluster that is associated with
    // that vertex

    if (vtx.size() == 0) return false;
    if (fcl2hits.size() == 0) return false;

    unsigned int iht = fcl2hits.size() - 1;
    unsigned short icl;
    float wire0 = (0.5 + fHits[fcl2hits[iht]].WireID().Wire);
    float dt;
    float thclus = std::atan(fScaleF * clpar[1]);

    for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if (vtx[iv].CTP != clCTP) continue;
      if (wire0 < vtx[iv].Wire) continue;
      if (wire0 > vtx[iv].Wire + 10) continue;
      dt = clpar[0] + (vtx[iv].Wire - wire0) * clpar[1] - vtx[iv].Time;
      if (std::abs(dt) > 10) continue;
      // cluster points to an US vertex. See if the angle is similar to
      // cluster associated with this vertex
      for (icl = 0; icl < tcl.size(); ++icl) {
        if (tcl[icl].CTP != clCTP) continue;
        if (tcl[icl].EndVtx != iv) continue;
        if (std::abs(tcl[icl].EndAng - thclus) < fKinkAngCut[pass]) return true;
      }
    }

    return false;

  } // CrawlVtxChk2()

  //////////////////////////////////////////
  bool ClusterCrawlerAlg::CrawlVtxChk(unsigned int kwire)
  {

    // returns true if the cluster is near a vertex on wire kwire
    if (vtx.size() == 0) return false;
    unsigned int iht = fcl2hits.size() - 1;
    float wire0 = (0.5 + fHits[fcl2hits[iht]].WireID().Wire);
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if (vtx[iv].CTP != clCTP) continue;
      if ((unsigned int)(0.5 + vtx[iv].Wire) != kwire) continue;
      if (std::abs(prtime - vtx[iv].Time) < 10) return true;
    }
    return false;
  } // CrawlVtxChk()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::VtxConstraint(unsigned int iwire,
                                        unsigned int ihit,
                                        unsigned int jwire,
                                        unsigned int& useHit,
                                        bool& doConstrain)
  {
    // checks hits on wire jwire to see if one is on a line between a US vertex
    // and the hit ihit on wire iwire. If one is found, doConstrain is set true
    // and the hit index is returned.
    doConstrain = false;
    if (vtx.size() == 0) return;
    // no vertices made yet on the first pass
    if (pass == 0) return;
    // skip if vertices were not requested to be made on the previous pass
    if (!fFindVertices[pass - 1]) return;

    if (jwire > WireHitRange.size() - 1) {
      mf::LogError("CC") << "VtxConstraint fed bad jwire " << jwire << " WireHitRange size "
                         << WireHitRange.size();
      return;
    }

    unsigned int jfirsthit = WireHitRange[jwire].first;
    unsigned int jlasthit = WireHitRange[jwire].second;
    for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if (vtx[iv].CTP != clCTP) continue;
      // vertex must be US of the cluster
      if (vtx[iv].Wire > jwire) continue;
      // but not too far US
      if (vtx[iv].Wire < jwire - 10) continue;
      recob::Hit const& hit = fHits[ihit];
      clpar[0] = hit.PeakTime();
      clpar[1] = (vtx[iv].Time - hit.PeakTime()) / (vtx[iv].Wire - iwire);
      clpar[2] = hit.WireID().Wire;
      float prtime = clpar[0] + clpar[1] * (jwire - iwire);
      for (unsigned int jhit = jfirsthit; jhit < jlasthit; ++jhit) {
        if (inClus[jhit] != 0) continue;
        const float tdiff = std::abs(fHits[jhit].TimeDistanceAsRMS(prtime));
        if (tdiff < 2.5) {
          useHit = jhit;
          doConstrain = true;
          return;
        }
      } // jhit
    }   // iv
  }     // VtxConstraint()

  /////////////////////////////////////////

  void ClusterCrawlerAlg::RefineVertexClusters(unsigned short iv)
  {

    // Try to attach or remove hits on the ends of vertex clusters

    std::vector<unsigned short> begClusters;
    std::vector<short> begdW;
    std::vector<unsigned short> endClusters;
    std::vector<short> enddW;

    unsigned int vWire = (unsigned int)(vtx[iv].Wire + 0.5);
    unsigned int vWireErr = (unsigned int)(2 * vtx[iv].WireErr);
    unsigned int vWireLo = vWire - vWireErr;
    unsigned int vWireHi = vWire + vWireErr;

    unsigned short icl, ii;
    short dW;
    bool needsWork = false;
    short maxdW = -100;
    short mindW = 100;
    for (icl = 0; icl < tcl.size(); ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != vtx[iv].CTP) continue;
      if (tcl[icl].BeginVtx == iv) {
        dW = vWire - tcl[icl].BeginWir;
        if (dW > maxdW) maxdW = dW;
        if (dW < mindW) mindW = dW;
        if (std::abs(dW) > 1) needsWork = true;
        // TODO: Check dTime also?
        begClusters.push_back(icl);
        begdW.push_back(dW);
      }
      if (tcl[icl].EndVtx == iv) {
        dW = vWire - tcl[icl].EndWir;
        if (dW > maxdW) maxdW = dW;
        if (dW < mindW) mindW = dW;
        if (std::abs(dW) > 1) needsWork = true;
        endClusters.push_back(icl);
        enddW.push_back(dW);
      }
    } // icl

    if (vtxprt)
      mf::LogVerbatim("CC") << "RefineVertexClusters: vertex " << iv << " needsWork " << needsWork
                            << " mindW " << mindW << " maxdW " << maxdW << " vWireErr " << vWireErr;

    if (!needsWork) return;

    // See if we can move the vertex within errors to reconcile the differences
    // without altering the clusters
    if (((unsigned int)std::abs(mindW) < vWireErr) && ((unsigned int)std::abs(maxdW) < vWireErr) &&
        std::abs(maxdW - mindW) < 2) {
      if (vtxprt) mf::LogVerbatim("CC") << " Move vtx wire " << vtx[iv].Wire;
      vtx[iv].Wire -= (float)(maxdW + mindW) / 2;
      if (vtxprt) mf::LogVerbatim("CC") << " to " << vtx[iv].Wire;
      // TODO: Fix the vertex time here if necessary
      vtx[iv].Fixed = true;
      // try to attach other clusters
      VertexCluster(iv);
      return;
    }

    // Check the vertex End clusters
    unsigned short newSize;
    for (ii = 0; ii < endClusters.size(); ++ii) {
      icl = endClusters[ii];
      if (vtxprt)
        mf::LogVerbatim("CC") << " endCluster " << tcl[icl].ID << " dW " << enddW[ii] << " vWire "
                              << vWire;
      if (tcl[icl].EndWir < vWire) {
        // vertex is DS of the cluster end -> remove hits
        TmpGet(icl);
        newSize = fcl2hits.size();
        for (auto hiter = fcl2hits.rbegin(); hiter < fcl2hits.rend(); ++hiter) {
          if (fHits[*hiter].WireID().Wire > vWire) break;
          --newSize;
        }
        // release the hits
        for (auto hiter = fcl2hits.begin(); hiter < fcl2hits.end(); ++hiter)
          inClus[*hiter] = 0;
        // shorten the cluster
        fcl2hits.resize(newSize);
        MakeClusterObsolete(icl);
        // fit
        FitCluster();
        clProcCode += 700;
        // store it
        TmpStore();
        tcl[tcl.size() - 1].EndVtx = iv;
        // update the vertex association
        if (vtxprt)
          mf::LogVerbatim("CC") << " modified cluster " << tcl[icl].ID << " -> "
                                << tcl[tcl.size() - 1].ID;
      } // tcl[icl].EndWir < vWire
      else if (tcl[icl].EndWir > vWire) {
        mf::LogVerbatim("CC") << "RefineVertexClusters: Write some EndVtx code";
      } //
    }   // ii endClusters

    if (begClusters.size() > 0)
      mf::LogVerbatim("CC") << "RefineVertexClusters: Write some BeginVtx code";

    if (mindW < 0 && maxdW > 0) {
      // vertex wire is in between the ends of the clusters
      // inspect the hits on both clusters near the vertex. The vertex should probably be on the hit
      // with the highest charge
      int vtxHit = -1;
      unsigned short clsBigChg = 0;
      float bigChg = 0;
      unsigned int iht;
      unsigned int ihit;
      // check the begClusters
      for (ii = 0; ii < begClusters.size(); ++ii) {
        icl = begClusters[ii];
        for (iht = 0; iht < tcl[icl].tclhits.size(); ++iht) {
          ihit = tcl[icl].tclhits[iht];
          if (fHits[ihit].Integral() > bigChg) {
            bigChg = fHits[ihit].Integral();
            vtxHit = ihit;
            clsBigChg = icl;
          }
          if (fHits[ihit].WireID().Wire < vWireLo) break;
        } // iht
      }   // ii
      // now check the endClusters
      for (ii = 0; ii < endClusters.size(); ++ii) {
        icl = endClusters[ii];
        for (iht = 0; iht < tcl[icl].tclhits.size(); ++iht) {
          ihit = tcl[icl].tclhits[tcl[icl].tclhits.size() - 1 - iht];
          if (fHits[ihit].Integral() > bigChg) {
            bigChg = fHits[ihit].Integral();
            vtxHit = ihit;
            clsBigChg = icl;
          }
          if (fHits[ihit].WireID().Wire > vWireHi) break;
        } // iht
      }   // ii
      if (vtxHit > 0) {
        if (vtxprt)
          mf::LogVerbatim("CC") << " moving vertex location to hit " << fHits[vtxHit].WireID().Wire
                                << ":" << (int)fHits[vtxHit].PeakTime() << " on cluster "
                                << tcl[clsBigChg].ID;
        vtx[iv].Wire = fHits[vtxHit].WireID().Wire;
        vtx[iv].Time = fHits[vtxHit].PeakTime();
        vtx[iv].Fixed = true;
      } // vtxHit > 0
    }   // mindW < 0 && maxdW > 0

    FitVtx(iv);

  } // RefineVertexClusters

  /////////////////////////////////////////
  bool ClusterCrawlerAlg::VtxClusterSplit()
  {

    // split clusters that cross vertices

    vtxprt = (fDebugPlane >= 0 && fDebugWire == 5555);

    if (vtxprt) mf::LogVerbatim("CC") << "VtxClusterSplit ";

    if (vtx.size() == 0) return false;
    unsigned short tclsize = tcl.size();
    if (tclsize < 2) return false;

    bool didit;
    bool didSomething = false;

    for (unsigned short icl = 0; icl < tclsize; ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != clCTP) continue;
      // ignore short clusters
      if (tcl[icl].tclhits.size() < 5) continue;
      // check vertices
      didit = false;
      for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if (vtx[ivx].CTP != clCTP) continue;
        // abandoned vertex?
        if (vtx[ivx].NClusters == 0) continue;
        // already assigned to this vertex?
        if (tcl[icl].BeginVtx == ivx) continue;
        if (tcl[icl].EndVtx == ivx) continue;
        // vertex wire in-between the cluster ends?
        if (vtx[ivx].Wire < tcl[icl].EndWir) continue;
        if (vtx[ivx].Wire > tcl[icl].BeginWir) continue;
        // vertex time in-between the cluster ends?
        float hiTime = tcl[icl].BeginTim;
        if (tcl[icl].EndTim > hiTime) hiTime = tcl[icl].EndTim;
        if (vtx[ivx].Time > hiTime + 5) continue;
        float loTime = tcl[icl].BeginTim;
        if (tcl[icl].EndTim < loTime) loTime = tcl[icl].EndTim;
        if (vtx[ivx].Time < loTime - 5) continue;
        // find the hit on the cluster that is closest to the vertex on the
        // DS side
        if (vtxprt)
          mf::LogVerbatim("CC") << " Chk cluster ID " << tcl[icl].ID << " with vertex " << ivx;
        short ihvx = -99;
        // nSplit is the index of the hit in the cluster where we will
        // split it if all requirements are met
        unsigned short nSplit = 0;
        unsigned short nLop = 0;
        unsigned int iht;
        for (unsigned short ii = tcl[icl].tclhits.size() - 1; ii > 0; --ii) {
          iht = tcl[icl].tclhits[ii];
          ++nLop;
          if (fHits[iht].WireID().Wire >= vtx[ivx].Wire) {
            nSplit = ii;
            if (vtxprt)
              mf::LogVerbatim("CC") << "Split cluster " << tcl[icl].ID << " at wire "
                                    << fHits[iht].WireID().Wire << " nSplit " << nSplit;
            ihvx = iht;
            break;
          }
        } // ii
        // found the wire. Now make a rough time cut
        if (ihvx < 0) continue;
        if (fabs(fHits[ihvx].PeakTime() - vtx[ivx].Time) > 10) continue;
        // check the angle between the crossing cluster icl and the
        // clusters that comprise the vertex.
        // First decide which end of cluster icl to use to define the angle
        float iclAng = 0.;
        if (nSplit > tcl[icl].tclhits.size() / 2) { iclAng = tcl[icl].EndAng; }
        else {
          iclAng = tcl[icl].BeginAng;
        }
        if (vtxprt) mf::LogVerbatim("CC") << " iclAng " << iclAng;
        // check angle wrt the the vertex clusters
        bool angOK = false;
        for (unsigned short jcl = 0; jcl < tclsize; ++jcl) {
          if (tcl[jcl].ID < 0) continue;
          if (tcl[jcl].CTP != clCTP) continue;
          if (tcl[jcl].BeginVtx == ivx) {
            if (fabs(tcl[jcl].BeginAng - iclAng) > 0.4) {
              // large angle difference. Set the flag
              angOK = true;
              break;
            }
          } // tcl[jcl].BeginVtx == ivx
          if (tcl[jcl].EndVtx == ivx) {
            if (fabs(tcl[jcl].EndAng - iclAng) > 0.4) {
              // large angle difference. Set the flag
              angOK = true;
              break;
            }
          } // tcl[jcl].EndVtx == ivx
        }   // jcl
        // time to split or chop
        if (angOK) {
          if (vtxprt) mf::LogVerbatim("CC") << "Split/Chop at pos " << nLop;
          if (nLop < 3) {
            // lop off hits at the US end
            // Put the cluster in the local arrays
            TmpGet(icl);
            for (unsigned short ii = 0; ii < nLop; ++ii) {
              unsigned int iht = fcl2hits[fcl2hits.size() - 1];
              inClus[iht] = 0;
              fcl2hits.pop_back();
            }
            // store it
            clProcCode += 1000;
            // declare this cluster obsolete
            MakeClusterObsolete(icl);
            // store the new one
            if (!TmpStore()) continue;
            unsigned short newcl = tcl.size() - 1;
            tcl[newcl].BeginVtx = tcl[icl].BeginVtx;
            tcl[newcl].EndVtx = ivx;
          }
          else {
            // split the cluster into two
            // correct the split position
            ++nSplit;
            if (SplitCluster(icl, nSplit, ivx)) {
              tcl[tcl.size() - 1].ProcCode += 1000;
              tcl[tcl.size() - 2].ProcCode += 1000;
            }
          }
          didit = true;
          didSomething = true;
        } // angOK
        if (didit) break;
      } // ivx
    }   // icl

    return didSomething;

  } // VtxClusterSplit()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::MergeHits(const unsigned int theHit, bool& didMerge)
  {
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

    if (theHit > fHits.size() - 1) { return; }

    recob::Hit const& hit = fHits[theHit];

    // don't bother trying to merge an already merged hit
    if (fHits[theHit].GoodnessOfFit() == 6666) {
      if (prt)
        mf::LogVerbatim("CC") << "MergeHits Trying to merge already merged hit "
                              << hit.WireID().Plane << ":" << hit.WireID().Wire << ":"
                              << (int)hit.PeakTime() << " Multiplicity " << hit.Multiplicity()
                              << " theHit " << theHit;
      return;
    }

    // don't merge if it isn't available
    if (!mergeAvailable[theHit]) { return; }

    if (hit.Multiplicity() < 2) return;

    // number of hits in this hit multiplet
    std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(theHit);

    // ensure that this is a high multiplicity hit:
    if (MultipletRange.second <= MultipletRange.first) return;

    // do a quick check to see how many hits are available to be merged
    unsigned short nAvailable = 0;
    unsigned short nInClus = 0;
    for (size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
      if (jht == theHit) continue;
      if (fHits[jht].GoodnessOfFit() == 6666) continue;
      if (inClus[jht] != 0) {
        ++nInClus;
        continue;
      }
      ++nAvailable;
    } // jht
    if (nAvailable == 0) return;
    // don't merge if any hit is used
    if (nInClus > 0) return;

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
    for (size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
      if (inClus[jht] < 0) continue;
      recob::Hit const& other_hit = fHits[jht];
      // error checking
      if ((other_hit.StartTick() != hit.StartTick()) || (other_hit.WireID() != hit.WireID())) {
        return;
      }
      if (other_hit.Multiplicity() != hit.Multiplicity()) { return; }
      // hit is not used by another cluster
      if (inClus[jht] != 0) continue;
      short arg = (short)(other_hit.PeakTimeMinusRMS(3));
      if (arg < loTime) loTime = arg;
      arg = (short)(other_hit.PeakTimePlusRMS(3));
      if (arg > hiTime) hiTime = arg;
      if (jht != theHit) ++nGaus;
      hitSep = std::abs(other_hit.PeakTime() - hit.PeakTime()) / other_hit.RMS();
      if (jht != theHit && hitSep < 3) ++nclose;
    } // jht
    // all hits in the multiplet other than this one used?
    if (nGaus < 2) return;

    // the hits in this multiplet will have this multiplicity from now on
    const short int NewMultiplicity = hit.Multiplicity() + 1 - nGaus;

    if (loTime < 0) loTime = 0;
    ++hiTime;
    // define a signal shape, fill it with zeros
    std::vector<double> signal(hiTime - loTime, 0.);
    // now add the Gaussians for each hit
    double chgsum = 0.;
    for (size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
      recob::Hit const& other_hit = fHits[jht];
      if (jht != theHit) {
        // hit used in another cluster
        if (inClus[jht] != 0) continue;
        // declare this hit obsolete
        inClus[jht] = -1;
      } // jht != theHit
      // add up the charge
      chgsum += other_hit.Integral();
      for (unsigned short time = loTime; time < hiTime; ++time) {
        unsigned short indx = time - loTime;
        double arg = (other_hit.PeakTime() - (double)time) / other_hit.RMS();
        signal[indx] += other_hit.PeakAmplitude() * exp(-0.5 * arg * arg);
      } // time
    }   // jj
    // find the average weighted time
    double sigsum = 0.;
    double sigsumt = 0.;
    for (unsigned short time = loTime; time < hiTime; ++time) {
      sigsum += signal[time - loTime];
      sigsumt += signal[time - loTime] * time;
    }
    if (sigsum == 0.) {
      //      mf::LogError("CC")<<"MergeHits: bad sum";
      return;
    }
    double aveTime = sigsumt / sigsum;
    // find the RMS
    sigsumt = 0.;
    for (unsigned short time = loTime; time < hiTime; ++time) {
      double dtime = time - aveTime;
      sigsumt += signal[time - loTime] * dtime * dtime;
    }
    const float RMS = std::sqrt(sigsumt / sigsum);
    // find the amplitude from the integrated charge and the RMS
    const float amplitude = chgsum * chgNorm / (2.507 * RMS);
    // modify the hit "in place" (actually completely overwrite it...)
    // TODO a lot of these quantities need revamp!!
    fHits[theHit] = recob::Hit(hit.Channel(),
                               hit.StartTick(),
                               hit.EndTick(),
                               aveTime, // peak_time
                               hit.SigmaPeakTime(),
                               RMS,       // rms
                               amplitude, // peak_amplitude
                               hit.SigmaPeakAmplitude(),
                               hit.SummedADC(),
                               chgsum, // hit_integral
                               hit.SigmaIntegral(),
                               NewMultiplicity, // multiplicity
                               0,               // local index
                               6666,            // GoodnessOfFit (flag for merged hit)
                               hit.DegreesOfFreedom(),
                               hit.View(),
                               hit.SignalType(),
                               hit.WireID());
    FixMultipletLocalIndices(MultipletRange.first, MultipletRange.second);
    didMerge = true;

  } // MergeHits()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::FixMultipletLocalIndices(size_t begin,
                                                   size_t end,
                                                   short int multiplicity /* = -1 */)
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
        if (inClus[iHit] < 0) continue;
        ++multiplicity;
      } // for
    }   // if no valid multiplicity is given

    // second pass: assign the correct multiplicity
    short int local_index = 0;
    for (size_t iHit = begin; iHit < end; ++iHit) {
      if (inClus[iHit] < 0) continue;

      // copy everything but overwrite the local index and multiplicity
      // TODO use a write wrapper!
      recob::Hit const& hit = fHits[iHit];
      fHits[iHit] = recob::Hit(hit.Channel(),
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
                               multiplicity, // multiplicity
                               local_index,  // local index
                               hit.GoodnessOfFit(),
                               hit.DegreesOfFreedom(),
                               hit.View(),
                               hit.SignalType(),
                               hit.WireID());

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
    if (tcl.size() < 2) return;

    // This code has some issues...
    return;

    vtxprt = (fDebugPlane == (int)plane && fDebugHit < 0);
    if (vtxprt) {
      mf::LogVerbatim("CC") << "FindStarVertices";
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
    for (unsigned short it1 = 0; it1 < tclsize - 1; ++it1) {
      // ignore obsolete clusters
      if (tcl[it1].ID < 0) continue;
      if (tcl[it1].CTP != clCTP) continue;
      // long dead-straight cluster?
      if (tcl[it1].tclhits.size() > 100) {
        dth = tcl[it1].BeginAng - tcl[it1].EndAng;
        if (std::abs(dth) < 0.1) continue;
      }
      es1 = tcl[it1].EndSlp;
      eth1 = tcl[it1].EndAng;
      ew1 = tcl[it1].EndWir;
      et1 = tcl[it1].EndTim;
      bw1 = tcl[it1].BeginWir;
      bt1 = tcl[it1].BeginTim;
      for (unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
        if (tcl[it2].ID < 0) continue;
        if (tcl[it2].CTP != clCTP) continue;
        // long dead-straight cluster?
        if (tcl[it2].tclhits.size() > 100) {
          dth = tcl[it2].BeginAng - tcl[it2].EndAng;
          if (std::abs(dth) < 0.05) continue;
        }
        es2 = tcl[it2].EndSlp;
        eth2 = tcl[it2].EndAng;
        ew2 = tcl[it2].EndWir;
        et2 = tcl[it2].EndTim;
        bw2 = tcl[it2].BeginWir;
        bt2 = tcl[it2].BeginTim;
        // require angle difference
        dth = std::abs(eth1 - eth2);
        if (dth < 0.3) continue;
        dsl = es2 - es1;
        fvw = (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
        // intersection within the wire boundaries
        if (fvw < ew1 || fvw > bw1) continue;
        if (fvw < ew2 || fvw > bw2) continue;
        if (vtxprt)
          mf::LogVerbatim("CC") << "Chk clusters " << tcl[it1].ID << " " << tcl[it2].ID
                                << " topo5 vtx wire " << fvw;
        // ensure that the intersection is close to the US end of the longer
        // cluster if it is more than 20 hits long
        if (tcl[it1].tclhits.size() > tcl[it2].tclhits.size()) {
          if (tcl[it1].tclhits.size() > 20) {
            nwirecut = 0.5 * tcl[it1].tclhits.size();
            if ((fvw - ew1) > nwirecut) continue;
          }
        }
        else {
          if (tcl[it2].tclhits.size() > 20) {
            nwirecut = 0.5 * tcl[it2].tclhits.size();
            if ((fvw - ew2) > nwirecut) continue;
          }
        }
        fvt = et1 + (fvw - ew1) * es1;
        // and time boundaries
        lotime = 9999;
        if (et1 < lotime) lotime = et1;
        if (bt1 < lotime) lotime = bt1;
        if (et2 < lotime) lotime = et2;
        if (bt2 < lotime) lotime = bt2;
        hitime = 0;
        if (et1 > hitime) hitime = et1;
        if (bt1 > hitime) hitime = bt1;
        if (et2 > hitime) hitime = et2;
        if (bt2 > hitime) hitime = bt2;
        if (fvt < lotime || fvt > hitime) continue;
        if (vtxprt)
          mf::LogVerbatim("CC") << " vertex time " << fvt << " lotime " << lotime << " hitime "
                                << hitime;
        unsigned int vw = (0.5 + fvw);
        // ensure that the vertex is near a hit on both clusters
        unsigned int pos1 = 0;
        for (unsigned short ii = 0; ii < tcl[it1].tclhits.size(); ++ii) {
          if (fHits[tcl[it1].tclhits[ii]].WireID().Wire <= vw) {
            if (std::abs(int(fHits[tcl[it1].tclhits[ii]].PeakTime() - fvt)) < 10) pos1 = ii;
            break;
          }
        } // ii
        // vertex is not near a hit on cluster 1
        if (pos1 == 0) continue;
        unsigned short pos2 = 0;
        for (unsigned short ii = 0; ii < tcl[it2].tclhits.size(); ++ii) {
          if (fHits[tcl[it2].tclhits[ii]].WireID().Wire <= vw) {
            if (std::abs(int(fHits[tcl[it2].tclhits[ii]].PeakTime() - fvt)) < 10) pos2 = ii;
            break;
          }
        } // ii
        // vertex is not near a hit on cluster 2
        if (pos2 == 0) continue;
        // ensure we aren't near an existing vertex
        for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if (std::abs(fvw - vtx[iv].Wire) < 2 && std::abs(fvt - vtx[iv].Time) < 10) continue;
        }
        // make a new vertex
        VtxStore newvx;
        newvx.Wire = fvw;
        newvx.WireErr = 1;
        newvx.Time = fvt;
        newvx.TimeErr = 1;
        newvx.Topo = 5;
        newvx.CTP = clCTP;
        newvx.Fixed = false;
        vtx.push_back(newvx);
        unsigned short ivx = vtx.size() - 1;
        if (vtxprt)
          mf::LogVerbatim("CC") << " new vertex " << ivx << " cluster " << tcl[it1].ID
                                << " split pos " << pos1;
        if (!SplitCluster(it1, pos1, ivx)) continue;
        tcl[tcl.size() - 1].ProcCode += 1000;
        tcl[tcl.size() - 2].ProcCode += 1000;
        if (vtxprt)
          mf::LogVerbatim("CC") << " new vertex " << ivx << " cluster " << tcl[it2].ID
                                << " split pos " << pos2;
        if (!SplitCluster(it2, pos2, ivx)) continue;
        tcl[tcl.size() - 1].ProcCode += 1000;
        tcl[tcl.size() - 2].ProcCode += 1000;
        FitVtx(ivx);
        break;
      } // it2
    }   // it1

    if (vtx.size() > vtxSizeIn) {
      // try to split other clusters
      VtxClusterSplit();
      // try to attach other clusters to it
      VertexCluster(vtx.size() - 1);
      FitAllVtx(clCTP);
    } // new vertex

    if (vtxprt) {
      mf::LogVerbatim("CC") << "Vertices " << vtx.size();
      for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if (vtx[iv].CTP != clCTP) continue;
        mf::LogVerbatim("CC") << "vtx " << iv << " wire " << vtx[iv].Wire << " time "
                              << (int)vtx[iv].Time << " NClusters " << vtx[iv].NClusters << " topo "
                              << vtx[iv].Topo;
      }
      PrintClusters();
    }

  } // FindStarVertices()

  //////////////////////////////////////////
  void ClusterCrawlerAlg::VertexCluster(unsigned short iv)
  {
    // try to attach clusters to the specified vertex
    if (vtx[iv].NClusters == 0) return;

    short dwb, dwe, dtb, dte;
    bool sigOK;

    for (unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != vtx[iv].CTP) continue;

      dwb = vtx[iv].Wire - tcl[icl].BeginWir;
      dtb = vtx[iv].Time - tcl[icl].BeginTim;
      dwe = vtx[iv].Wire - tcl[icl].EndWir;
      dte = vtx[iv].Time - tcl[icl].EndTim;

      float drb = dwb * dwb + dtb * dtb;
      float dre = dwe * dwe + dte * dte;

      bool bCloser = (drb < dre);

      // ignore clusters in showers
      if (bCloser) {
        if (tcl[icl].BeginChgNear > fChgNearCut) continue;
      }
      else {
        if (tcl[icl].EndChgNear > fChgNearCut) continue;
      }

      if (vtxprt)
        mf::LogVerbatim("CC") << "VertexCluster: Try icl ID " << tcl[icl].ID << " w vtx " << iv
                              << " dwb " << dwb << " dwe " << dwe << " drb " << drb << " dre "
                              << dre << " Begin closer? " << bCloser;

      if (tcl[icl].BeginVtx < 0 && bCloser && dwb > -3 && dwb < 3 && tcl[icl].EndVtx != iv) {
        sigOK = ChkSignal(tcl[icl].BeginWir, tcl[icl].BeginTim, vtx[iv].Wire, vtx[iv].Time);
        if (vtxprt)
          mf::LogVerbatim("CC") << " Attach cluster Begin to vtx? " << iv << " sigOK " << sigOK;
        if (sigOK) {
          if (vtxprt)
            mf::LogVerbatim("CC") << " check ClusterVertexChi " << ClusterVertexChi(icl, 0, iv);
          if (ClusterVertexChi(icl, 0, iv) < fVertex2DCut) {
            // do a fit and check the vertex error
            tcl[icl].BeginVtx = iv;
            FitVtx(iv);
            if (vtx[iv].ChiDOF > fVertex2DCut || vtx[iv].WireErr > fVertex2DWireErrCut) {
              tcl[icl].BeginVtx = -99;
              FitVtx(iv);
            }
          } // good DoCA
        }   // sigOK
      }     // check BEGIN

      if (tcl[icl].EndVtx < 0 && !bCloser && dwe > -3 && dwe < 3 && tcl[icl].BeginVtx != iv) {
        sigOK = ChkSignal(tcl[icl].EndWir, tcl[icl].EndTim, vtx[iv].Wire, vtx[iv].Time);
        if (vtxprt)
          mf::LogVerbatim("CC") << " Attach cluster End to vtx? " << iv << " sigOK " << sigOK;
        if (sigOK) {
          if (vtxprt)
            mf::LogVerbatim("CC") << " check ClusterVertexChi " << ClusterVertexChi(icl, 1, iv);
          if (ClusterVertexChi(icl, 1, iv) < 3) {
            // do a fit and check the vertex error
            tcl[icl].EndVtx = iv;
            FitVtx(iv);
            if (vtx[iv].ChiDOF > fVertex2DCut || vtx[iv].WireErr > fVertex2DWireErrCut) {
              tcl[icl].EndVtx = -99;
              FitVtx(iv);
            }
          } // good DoCA
        }   // sigOK
      }     // check END
    }       // icl
  }         // VertexCluster

  //////////////////////////////////////////
  void ClusterCrawlerAlg::FitAllVtx(CTP_t inCTP)
  {

    for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if (vtx[iv].CTP != inCTP) continue;
      FitVtx(iv);
    }

  } // FitAllVtx

  /////////////////////////////////////////
  void ClusterCrawlerAlg::FindVertices()
  {
    // try to make 2D vertices

    if (tcl.size() < 2) return;

    // sort clusters starting with the longest
    std::vector<CluLen> clulens;
    CluLen clulen;
    for (unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != clCTP) continue;
      if (tcl[icl].BeginVtx >= 0 && tcl[icl].EndVtx >= 0) continue;
      clulen.index = icl;
      clulen.length = tcl[icl].tclhits.size();
      clulens.push_back(clulen);
    }
    if (empty(clulens)) return;
    std::sort(clulens.begin(), clulens.end(), greaterThan);

    float nwires = fNumWires;
    float maxtime = fMaxTime;

    unsigned short vtxSizeIn = vtx.size();

    vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
    if (vtxprt) {
      mf::LogVerbatim("CC") << "FindVertices plane " << plane << " pass " << pass;
      PrintClusters();
    }

    float es1 = 0, eth1 = 0, ew1 = 0, et1 = 0;
    float bs1 = 0, bth1 = 0, bw1 = 0, bt1 = 0;
    float es2 = 0, eth2 = 0, ew2 = 0, et2 = 0;
    float bs2 = 0, bth2 = 0, bw2 = 0, bt2 = 0;
    float dth = 0, dsl = 0, fvw = 0, fvt = 0;
    float angcut = 0;
    unsigned int vw = 0, pass1, pass2, ii1, it1, ii2, it2;
    bool SigOK = false;
    for (ii1 = 0; ii1 < clulens.size() - 1; ++ii1) {
      it1 = clulens[ii1].index;
      es1 = tcl[it1].EndSlp;
      eth1 = tcl[it1].EndAng;
      ew1 = tcl[it1].EndWir;
      et1 = tcl[it1].EndTim;
      bs1 = tcl[it1].BeginSlp;
      bth1 = tcl[it1].BeginAng;
      bw1 = tcl[it1].BeginWir;
      bt1 = tcl[it1].BeginTim;
      pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
      for (ii2 = ii1 + 1; ii2 < clulens.size(); ++ii2) {
        it2 = clulens[ii2].index;
        // try to attach cluster to existing vertices at either end
        ClusterVertex(it2);
        // ignore if both clusters are short
        if (tcl[it1].tclhits.size() < 5 && tcl[it2].tclhits.size() < 5) continue;
        es2 = tcl[it2].EndSlp;
        eth2 = tcl[it2].EndAng;
        ew2 = tcl[it2].EndWir;
        et2 = tcl[it2].EndTim;
        bs2 = tcl[it2].BeginSlp;
        bth2 = tcl[it2].BeginAng;
        bw2 = tcl[it2].BeginWir;
        bt2 = tcl[it2].BeginTim;
        pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
        if (pass1 < pass2) { angcut = fKinkAngCut[pass2]; }
        else {
          angcut = fKinkAngCut[pass1];
        }
        // topo 1: check for vtx US of the ends of both clusters
        dth = fabs(eth1 - eth2);
        if (tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0 && dth > 0.1) {
          dsl = es2 - es1;
          // find vertex wire and vertex time in float precision (fvw, fvt)
          fvw = (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
          // vertex wire in the detector?
          if (fvw > 0. && fvw < nwires) {
            // require vtx in the range of wires with hits AND
            vw = (0.5 + fvw);
            // vtx US of both clusters AND
            // vtx not too far US of both clusters
            if (vw >= fFirstWire - 1 && fvw <= ew1 + 3 && fvw <= ew2 + 3 && fvw > (ew1 - 10) &&
                fvw > (ew2 - 10)) {
              fvt = et1 + (fvw - ew1) * es1;
              if (vtxprt)
                mf::LogVerbatim("CC")
                  << "Chk clusters topo1 " << tcl[it1].ID << " " << tcl[it2].ID << "  vtx wire "
                  << vw << " time " << (int)fvt << " dth " << dth;
              if (fvt > 0 && fvt < maxtime) {
                // Vertex wire US of cluster ends and time in the detector.
                // Check for signal at the vertex position and adjust the vertex by 1 wire
                // if necessary
                SigOK = ChkSignal(vw, fvt, vw, fvt);
                if (!SigOK) {
                  fvw += 1.;
                  vw = (0.5 + fvw);
                  SigOK = ChkSignal(vw, fvt, vw, fvt);
                }
                // Check this against existing vertices and update
                if (SigOK) ChkVertex(fvw, fvt, it1, it2, 1);
              } // fvt in detector
            }   // vw topo 1 check
          }     // fvw in detector
        }       // topo 1
                // topo 2: check for vtx US of it1 and DS of it2
        dth = std::abs(eth1 - bth2);
        if (tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0 && dth > angcut) {
          dsl = bs2 - es1;
          if (fabs(ew1 - bw2) < 3 && fabs(et1 - bt2) < 500) { fvw = 0.5 * (ew1 + bw2); }
          else {
            fvw = (et1 - ew1 * es1 - bt2 + bw2 * bs2) / dsl;
          }
          if (fvw > 0 && fvw < nwires) {
            // vertex wire in the detector
            vw = (0.5 + fvw);
            // require vtx US of cluster 1 End AND
            //         vtx DS of cluster 2 Begin
            if (fvw <= ew1 + 2 && fvw >= bw2 - 2) {
              fvt = et1 + (vw - ew1) * es1;
              fvt += bt2 + (vw - bw2) * bs2;
              fvt /= 2;
              if (vtxprt)
                mf::LogVerbatim("CC")
                  << "Chk clusters topo2 " << tcl[it1].ID << " " << tcl[it2].ID << " vtx wire "
                  << vw << " time " << (int)fvt << " dth " << dth;
              if (fvt > 0. && fvt < maxtime) {
                ChkVertex(fvw, fvt, it1, it2, 2);
              } // fvt in detector
            }   // vw topo 2 check
          }     // fvw in detector
        }       // topo 2
                // topo 3: check for vtx DS of it1 and US of it2
        dth = std::abs(bth1 - eth2);
        if (tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0 && dth > angcut) {
          dsl = bs1 - es2;
          if (fabs(bw1 - ew2) < 3 && fabs(bt1 - et2) < 500) { fvw = 0.5 * (bw1 + ew2); }
          else {
            fvw = (et2 - ew2 * es2 - bt1 + bw1 * bs1) / dsl;
          }
          if (fvw > 0 && fvw < nwires) {
            vw = (0.5 + fvw);
            // require vtx US of cluster 2 Begin AND
            //         vtx DS of cluster 1 End
            if (fvw <= ew2 + 2 && fvw >= bw1 - 2) {
              fvt = et2 + (fvw - ew2) * es2;
              fvt += bt1 + (fvw - bw1) * es1;
              fvt /= 2;
              if (vtxprt)
                mf::LogVerbatim("CC")
                  << "Chk clusters topo3 " << tcl[it1].ID << " " << tcl[it2].ID << " vtx wire "
                  << vw << " time " << (int)fvt << " dth " << dth;
              if (fvt > 0. && fvt < maxtime) {
                ChkVertex(fvw, fvt, it1, it2, 3);
              } // fvt in detector
            }   // vw topo 3 check
          }     // fvw in detector
        }       // topo 3
                // topo 4: check for vtx DS of it1 and DS of it2
        dth = std::abs(bth1 - bth2);
        if (tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0 && dth > 0.1) {
          dsl = bs2 - bs1;
          // find vertex wire and vertex time in float precision (fvw, fvt)
          // convert to integer if within the detector (vw, vt)
          fvw = (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl;
          if (vtxprt)
            mf::LogVerbatim("CC") << "Chk clusters topo4 " << bw1 << ":" << (int)bt1 << " " << bw2
                                  << ":" << (int)bt2 << " fvw " << fvw << " nwires " << nwires;
          if (fvw > 0 && fvw < nwires) {
            // vertex wire in the detector
            vw = (0.5 + fvw);
            // require vtx in the range of wires with hits AND vtx DS of both clusters AND
            // vtx not too far DS of both clusters
            float dwcut = 10;
            if (tcl[it1].tclhits.size() < 2 * dwcut) dwcut = tcl[it1].tclhits.size() / 2;
            if (tcl[it2].tclhits.size() < 2 * dwcut) dwcut = tcl[it2].tclhits.size() / 2;
            if (fvw <= fLastWire + 1 && fvw >= bw1 - dwcut && fvw <= bw1 + dwcut &&
                fvw >= bw2 - dwcut && fvw <= bw1 + dwcut) {
              fvt = bt1 + (fvw - bw1) * bs1;
              if (vtxprt)
                mf::LogVerbatim("CC")
                  << " vtx wire " << vw << " time " << fvt << " dth " << dth << " dwcut " << dwcut;
              if (fvt > 0. && fvt < maxtime) {
                // vertex wire US of cluster ends and time in the detector
                // Check for signal at the vertex position and adjust the vertex by 1 wire
                // if necessary
                SigOK = ChkSignal(vw, fvt, vw, fvt);
                if (!SigOK) {
                  fvw -= 1.;
                  vw = (0.5 + fvw);
                  SigOK = ChkSignal(vw, fvt, vw, fvt);
                }
                // Check this against existing vertices and update
                if (SigOK) ChkVertex(fvw, fvt, it1, it2, 4);
              } // fvt in detector
            }   // vw topo 4 check
          }     // fvw in detector
        }       // topo4
      }         // it2
    }           // it1

    // Fix vertices where both ends of a cluster are assigned to the
    // same vertex. This can happen with short clusters
    for (unsigned short it = 0; it < tcl.size(); ++it) {
      if (tcl[it].ID < 0) continue;
      if (tcl[it].CTP != clCTP) continue;
      if (tcl[it].BeginVtx > -1 && tcl[it].BeginVtx == tcl[it].EndVtx) {
        unsigned short iv = tcl[it].BeginVtx;
        float dwir = tcl[it].BeginWir - vtx[iv].Wire;
        float dtim = tcl[it].BeginTim - vtx[iv].Time;
        float rbeg = dwir * dwir + dtim * dtim;
        dwir = tcl[it].EndWir - vtx[iv].Wire;
        dtim = tcl[it].EndTim - vtx[iv].Time;
        float rend = dwir * dwir + dtim * dtim;
        if (rend < rbeg) { tcl[it].EndVtx = -99; }
        else {
          tcl[it].BeginVtx = -99;
        }
      } // tcl[it].BeginVtx == tcl[it].EndVtx
    }   // it

    if (vtx.size() > vtxSizeIn) FitAllVtx(clCTP);

    // "delete" any vertices that have only one cluster
    for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if (vtx[ivx].CTP != clCTP) continue;
      if (vtx[ivx].NClusters == 1) {
        vtx[ivx].NClusters = 0;
        for (unsigned short icl = 0; icl < tcl.size(); ++icl) {
          if (tcl[icl].CTP != clCTP) continue;
          if (tcl[icl].BeginVtx == ivx) {
            tcl[icl].BeginVtx = -99;
            break;
          }
          if (tcl[icl].EndVtx == ivx) {
            tcl[icl].EndVtx = -99;
            break;
          }
        } // icl
      }   // vtx[ivx].NClusters == 1
    }     // ivx

  } // FindVertices()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::ClusterVertex(unsigned short it)
  {
    // try to attach cluster it to an existing vertex

    if (vtx.size() == 0) return;

    unsigned short iv, jv;
    short dwib, dwjb, dwie, dwje;
    bool matchEnd, matchBegin;

    for (iv = 0; iv < vtx.size(); ++iv) {
      // ignore vertices in the wrong cryostat/TPC/Plane
      if (vtx[iv].CTP != clCTP) continue;
      // ignore deleted vertices
      if (vtx[iv].NClusters == 0) continue;
      // determine which end to match - begin or end. Handle short tracks
      matchEnd = false;
      matchBegin = false;
      if (tcl[it].tclhits.size() < 6) {
        // See which end is closer to this vertex vs other vertices
        dwib = std::abs(vtx[iv].Wire - tcl[it].BeginWir);
        if (dwib > 2) dwib = 2;
        dwie = std::abs(vtx[iv].Wire - tcl[it].EndWir);
        if (dwie > 2) dwie = 2;
        dwjb = 999;
        dwje = 999;
        for (jv = 0; jv < vtx.size(); ++jv) {
          if (iv == jv) continue;
          if (std::abs(vtx[jv].Time - tcl[it].BeginTim) < 50) {
            if (std::abs(vtx[jv].Wire - tcl[it].BeginWir) < dwjb)
              dwjb = std::abs(vtx[jv].Wire - tcl[it].BeginWir);
          } // std::abs(vtx[jv].Time - tcl[it].BeginTim) < 50
          if (std::abs(vtx[jv].Time - tcl[it].EndTim) < 50) {
            if (std::abs(vtx[jv].Wire - tcl[it].EndWir) < dwje)
              dwje = std::abs(vtx[jv].Wire - tcl[it].EndWir);
          } // std::abs(vtx[jv].Time - tcl[it].EndTim) < 50
        }   // jv
        matchEnd = tcl[it].EndVtx != iv && dwie < 3 && dwie < dwje && dwie < dwib;
        matchBegin = tcl[it].BeginVtx != iv && dwib < 3 && dwib < dwjb && dwib < dwie;
      }
      else {
        matchEnd = tcl[it].EndVtx < 0 && vtx[iv].Wire <= tcl[it].EndWir + 2;
        matchBegin = tcl[it].BeginVtx < 0 && vtx[iv].Wire >= tcl[it].BeginWir - 2;
      }
      if (matchEnd) {
        if (vtxprt)
          mf::LogVerbatim("CC") << " Match End chi " << ClusterVertexChi(it, 1, iv) << " to vtx "
                                << iv << " signalOK "
                                << ChkSignal(
                                     vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim);
        if (ClusterVertexChi(it, 1, iv) < fVertex2DCut &&
            ChkSignal(vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim)) {
          // good match
          tcl[it].EndVtx = iv;
          // re-fit it
          FitVtx(iv);
          if (vtxprt)
            mf::LogVerbatim("CC") << " Add End " << tcl[it].ID << " to vtx " << iv << " NClusters "
                                  << vtx[iv].NClusters;
          if (vtx[iv].ChiDOF < fVertex2DCut && vtx[iv].WireErr < fVertex2DWireErrCut) return;
          if (vtxprt)
            mf::LogVerbatim("CC") << " Bad fit. ChiDOF " << vtx[iv].ChiDOF << " WireErr "
                                  << vtx[iv].WireErr << " Undo it";
          tcl[it].EndVtx = -99;
          FitVtx(iv);
        } // tChi < 3
      }   // matchEnd

      if (matchBegin) {
        if (vtxprt)
          mf::LogVerbatim("CC") << " Match Begin chi " << ClusterVertexChi(it, 0, iv) << " to vtx "
                                << iv << " signalOK "
                                << ChkSignal(vtx[iv].Wire,
                                             vtx[iv].Time,
                                             tcl[it].BeginWir,
                                             tcl[it].BeginTim);
        if (ClusterVertexChi(it, 0, iv) < fVertex2DCut &&
            ChkSignal(vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim)) {
          // good match
          tcl[it].BeginVtx = iv;
          // re-fit it
          FitVtx(iv);
          if (vtxprt)
            mf::LogVerbatim("CC") << " Add Begin " << tcl[it].ID << " to vtx " << iv
                                  << " NClusters " << vtx[iv].NClusters;
          if (vtx[iv].ChiDOF < fVertex2DCut && vtx[iv].WireErr < fVertex2DWireErrCut) return;
          if (vtxprt)
            mf::LogVerbatim("CC") << " Bad fit. ChiDOF " << vtx[iv].ChiDOF << " WireErr "
                                  << vtx[iv].WireErr << " Undo it";
          tcl[it].BeginVtx = -99;
          FitVtx(iv);
        } // tChi < 3
      }   // matchBegin
    }     // iv
  }       // ClusterVertex()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::ChkVertex(float fvw,
                                    float fvt,
                                    unsigned short it1,
                                    unsigned short it2,
                                    short topo)
  {
    // Checks the vertex (vw, fvt) against the existing set of vertices.
    // If there a match, clusters it1 and/or it2 are associated with it
    // if there is signal between the existing vertex and the start of
    // the cluster. The topo flag indicates the type of vertex that was
    // found: 1 = US of it1 && US of it2, 2 = US of it1 && DS of it2,
    // 3 = DS of it1 && US of it2, 4 = DS of it1 and DS of it2.
    // didit is set true if a cluster is attached to a (new or old) vertex

    if (vtxprt)
      mf::LogVerbatim("CC") << " ChkVertex " << tcl[it1].EndWir << ":" << (int)tcl[it1].EndTim
                            << " - " << tcl[it1].BeginWir << ":" << (int)tcl[it1].BeginTim
                            << " and " << tcl[it2].EndWir << ":" << (int)tcl[it2].EndTim << " - "
                            << tcl[it2].BeginWir << ":" << (int)tcl[it2].BeginTim << " topo "
                            << topo << " fvw " << fvw << " fvt " << fvt;

    unsigned int vw = (unsigned int)(0.5 + fvw);
    unsigned int iht;

    // don't make vertices using short tracks that are close to
    // long straight clusters. These are most likely due to numerous
    // short delta ray clusters
    if (tcl[it1].tclhits.size() < 10 && tcl[it2].tclhits.size() < 10) {
      for (unsigned short it = 0; it < tcl.size(); ++it) {
        if (it == it1 || it == it2) continue;
        if (tcl[it].ID < 0) continue;
        if (tcl[it].CTP != clCTP) continue;
        if (tcl[it].tclhits.size() < 100) continue;
        if (std::abs(tcl[it].BeginAng - tcl[it].EndAng) > 0.1) continue;
        // don't reject because it is near the end of a long cluster
        if (vw < tcl[it].EndWir + 5) continue;
        if (vw > tcl[it].BeginWir - 5) continue;
        for (unsigned short ii = 0; ii < tcl[it].tclhits.size(); ++ii) {
          iht = tcl[it].tclhits[ii];
          if (fHits[iht].WireID().Wire <= vw) {
            if (std::abs(fHits[iht].PeakTime() - fvt) < 60) return;
          } // fHits[iht].WireID().Wire <= vWire
        }   // ii
      }     // it
    }

    // check vertex and clusters for proximity to existing vertices
    unsigned short nFitOK = 0;
    for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if (vtx[iv].CTP != clCTP) continue;
      // make this a really loose cut since the errors on the prospective vertex aren't known yet
      if (PointVertexChi(fvw, fvt, iv) > 300) continue;
      if (vtxprt)
        mf::LogVerbatim("CC") << " vtx " << iv << " PointVertexChi "
                              << PointVertexChi(fvw, fvt, iv);
      // got a match. Check the appropriate cluster end and attach
      if ((topo == 1 || topo == 2) && tcl[it1].EndVtx < 0) {
        if (ChkSignal(vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim)) {
          // try to fit
          tcl[it1].EndVtx = iv;
          FitVtx(iv);
          if (vtxprt)
            mf::LogVerbatim("CC") << " FitVtx " << iv << " WireErr " << vtx[iv].WireErr
                                  << " ChiDOF " << vtx[iv].ChiDOF;
          if (vtx[iv].WireErr < fVertex2DWireErrCut && vtx[iv].ChiDOF < 5) { ++nFitOK; }
          else {
            // bad fit
            tcl[it1].EndVtx = -99;
            FitVtx(iv);
          } // check fit
        }   // ChkSignal
      }
      else if ((topo == 3 || topo == 4) && tcl[it1].BeginVtx < 0) {
        if (ChkSignal(vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim)) {
          tcl[it1].BeginVtx = iv;
          FitVtx(iv);
          if (vtxprt)
            mf::LogVerbatim("CC") << " FitVtx " << iv << " WireErr " << vtx[iv].WireErr
                                  << " ChiDOF " << vtx[iv].ChiDOF;
          if (vtx[iv].WireErr < fVertex2DWireErrCut && vtx[iv].ChiDOF < 5) { ++nFitOK; }
          else {
            // bad fit
            tcl[it1].BeginVtx = -99;
            FitVtx(iv);
          } // bad fit
        }   // ChkSignal
      }     // cluster it2
      if ((topo == 1 || topo == 3) && tcl[it2].EndVtx < 0) {
        if (ChkSignal(vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim)) {
          tcl[it2].EndVtx = iv;
          FitVtx(iv);
          if (vtxprt)
            mf::LogVerbatim("CC") << " FitVtx " << iv << " WireErr " << vtx[iv].WireErr
                                  << " ChiDOF " << vtx[iv].ChiDOF;
          if (vtx[iv].WireErr < fVertex2DWireErrCut && vtx[iv].ChiDOF < 5) { ++nFitOK; }
          else {
            // bad fit
            tcl[it2].EndVtx = -99;
            FitVtx(iv);
          } // bad fit
        }   // ChkSignal
      }
      else if ((topo == 2 || topo == 4) && tcl[it2].BeginVtx < 0) {
        if (ChkSignal(vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim)) {
          tcl[it2].BeginVtx = iv;
          FitVtx(iv);
          if (vtxprt)
            mf::LogVerbatim("CC") << " FitVtx " << iv << " WireErr " << vtx[iv].WireErr
                                  << " ChiDOF " << vtx[iv].ChiDOF;
          if (vtx[iv].WireErr < fVertex2DWireErrCut && vtx[iv].ChiDOF < 5) { ++nFitOK; }
          else {
            // bad fit
            tcl[it2].BeginVtx = -99;
            FitVtx(iv);
          } // bad fit
        }   // ChkSignal
      }     // cluster it2
      if (nFitOK > 0) {
        if (vtxprt) mf::LogVerbatim("CC") << " Attached " << nFitOK << " clusters to vertex " << iv;
        return;
      }
    } // iv

    // no match to existing vertices. Ensure that there is a wire signal between
    // the vertex and the appropriate ends of the clusters
    bool SigOK = false;
    if (topo == 1 || topo == 2) { SigOK = ChkSignal(vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim); }
    else {
      SigOK = ChkSignal(vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim);
    }
    if (!SigOK) return;

    if (topo == 1 || topo == 3) { SigOK = ChkSignal(vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim); }
    else {
      SigOK = ChkSignal(vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim);
    }
    if (!SigOK) return;

    VtxStore newvx;
    newvx.Wire = vw;
    newvx.Time = fvt;
    newvx.Topo = topo;
    newvx.CTP = clCTP;
    newvx.Fixed = false;
    vtx.push_back(newvx);
    unsigned short iv = vtx.size() - 1;
    if (topo == 1 || topo == 2) {
      if (tcl[it1].EndVtx >= 0) {
        vtx.pop_back();
        return;
      }
      tcl[it1].EndVtx = iv;
    }
    else {
      if (tcl[it1].BeginVtx >= 0) {
        vtx.pop_back();
        return;
      }
      tcl[it1].BeginVtx = iv;
    }
    if (topo == 1 || topo == 3) {
      if (tcl[it2].EndVtx >= 0) {
        vtx.pop_back();
        return;
      }
      tcl[it2].EndVtx = iv;
    }
    else {
      if (tcl[it2].BeginVtx >= 0) {
        vtx.pop_back();
        return;
      }
      tcl[it2].BeginVtx = iv;
    }
    // fit it
    FitVtx(iv);
    // reject it if the fit is bad
    if (vtx[iv].ChiDOF < 5 && vtx[iv].WireErr < fVertex2DWireErrCut && vtx[iv].TimeErr < 20) {
      if (vtxprt)
        mf::LogVerbatim("CC") << " New vtx " << iv << " in plane " << plane << " topo " << topo
                              << " cls " << tcl[it1].ID << " " << tcl[it2].ID << " W:T " << (int)vw
                              << ":" << (int)fvt << " NClusters " << vtx[iv].NClusters;
      // Try to refine the cluster hits and vertex position
      if (fRefineVertexClusters) RefineVertexClusters(iv);
    }
    else {
      if (vtxprt)
        mf::LogVerbatim("CC") << " Bad vtx fit " << vtx[iv].ChiDOF << " wire err "
                              << vtx[iv].WireErr << " TimeErr " << vtx[iv].TimeErr;
      // clobber the vertex and references to it
      vtx.pop_back();
      if (tcl[it1].BeginVtx == iv) tcl[it1].BeginVtx = -99;
      if (tcl[it1].EndVtx == iv) tcl[it1].EndVtx = -99;
      if (tcl[it2].BeginVtx == iv) tcl[it2].BeginVtx = -99;
      if (tcl[it2].EndVtx == iv) tcl[it2].EndVtx = -99;
    } // bad fit

  } // ChkVertex()

  /////////////////////////////////////////
  bool ClusterCrawlerAlg::ChkSignal(unsigned int iht, unsigned int jht)
  {
    if (iht > fHits.size() - 1) return false;
    if (jht > fHits.size() - 1) return false;
    unsigned int wire1 = fHits[iht].WireID().Wire;
    float time1 = fHits[iht].PeakTime();
    unsigned int wire2 = fHits[jht].WireID().Wire;
    float time2 = fHits[jht].PeakTime();
    return ChkSignal(wire1, time1, wire2, time2);
  } // ChkSignal

  /////////////////////////////////////////
  bool ClusterCrawlerAlg::ChkSignal(unsigned int wire1,
                                    float time1,
                                    unsigned int wire2,
                                    float time2)
  {
    // returns  true if there is a signal on the line between
    // (wire1, time1) and (wire2, time2).
    // Be sure to set time1 < time2 if you are checking for signals on a single wire

    // Gaussian amplitude in bins of size 0.15
    const float gausAmp[20] = {1,    0.99, 0.96, 0.90, 0.84, 0.75, 0.67, 0.58, 0.49, 0.40,
                               0.32, 0.26, 0.20, 0.15, 0.11, 0.08, 0.06, 0.04, 0.03, 0.02};

    // return true if fMinAmp is set ignore wire signal checking in this plane
    if (fMinAmp[plane] <= 0) return true;

    // get the begin and end right
    unsigned int wireb = wire1;
    float timeb = time1;
    unsigned int wiree = wire2;
    float timee = time2;
    // swap them?
    if (wiree > wireb) {
      wireb = wire2;
      timeb = time2;
      wiree = wire1;
      timee = time1;
    }
    if (wiree < fFirstWire || wiree > fLastWire) return false;
    if (wireb < fFirstWire || wireb > fLastWire) return false;

    int wire0 = wiree;
    // checking a single wire?
    float slope = 0;
    bool oneWire = false;
    float prTime, prTimeLo = 0, prTimeHi = 0;
    if (wireb == wiree) {
      oneWire = true;
      if (time1 < time2) {
        prTimeLo = time1;
        prTimeHi = time2;
      }
      else {
        prTimeLo = time2;
        prTimeHi = time1;
      }
    }
    else {
      slope = (timeb - timee) / (wireb - wiree);
    }

    int bin;
    unsigned short nmissed = 0;
    for (unsigned int wire = wiree; wire < wireb + 1; ++wire) {
      if (oneWire) { prTime = (prTimeLo + prTimeHi) / 2; }
      else {
        prTime = timee + (wire - wire0) * slope;
      }
      // skip dead wires
      if (WireHitRange[wire].first == -1) continue;
      // no hits on this wire
      if (WireHitRange[wire].first == -2) ++nmissed;
      unsigned int firsthit = WireHitRange[wire].first;
      unsigned int lasthit = WireHitRange[wire].second;
      float amp = 0;
      for (unsigned int khit = firsthit; khit < lasthit; ++khit) {
        if (oneWire) {
          // TODO: This sometimes doesn't work with overlapping hits
          // A not totally satisfactory solution
          if (prTime < fHits[khit].StartTick()) continue;
          if (prTime > fHits[khit].EndTick()) continue;
          return true;
        }
        else {
          // skip checking if we are far away from prTime on the positive side
          if (fHits[khit].PeakTime() - prTime > 500) continue;
          bin = std::abs(fHits[khit].PeakTime() - prTime) / fHits[khit].RMS();
          bin /= 0.15;
          if (bin > 19) continue;
          if (bin < 0) continue;
          // add amplitude from all hits
          amp += fHits[khit].PeakAmplitude() * gausAmp[bin];
        }
      } // khit
      if (amp < fMinAmp[plane]) ++nmissed;
    } // wire
    if (nmissed > fAllowNoHitWire) return false;
    return true;

  } // ChkSignal()

  /////////////////////////////////////////
  bool ClusterCrawlerAlg::SplitCluster(unsigned short icl, unsigned short pos, unsigned short ivx)
  {
    // split cluster icl into two clusters

    if (tcl[icl].ID < 0) return false;

    // declare icl obsolete
    MakeClusterObsolete(icl);

    unsigned short ii, iclnew;
    unsigned int iht;

    if (pos > 2) {
      // Have enough hits to make a cluster at the Begin end
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
      for (ii = 0; ii < pos; ++ii) {
        iht = tcl[icl].tclhits[ii];
        fcl2hits.push_back(iht);
      }
      // determine the pass in which this cluster was created
      pass = tcl[icl].ProcCode - 10 * (tcl[icl].ProcCode / 10);
      if (pass > fNumPass - 1) pass = fNumPass - 1;
      // fit the end hits
      FitCluster();
      clEndSlp = clpar[1];
      clEndSlpErr = clparerr[1];
      clEndAng = std::atan(fScaleF * clEndSlp);
      // find the charge at the end
      FitClusterChg();
      clEndChg = fAveChg;
      if (!TmpStore()) {
        RestoreObsoleteCluster(icl);
        return false;
      }
      // associate the End with the supplied vertex
      iclnew = tcl.size() - 1;
      tcl[iclnew].EndVtx = ivx;
      tcl[iclnew].BeginVtx = tcl[icl].BeginVtx;
      if (vtxprt)
        mf::LogVerbatim("CC") << "SplitCluster made cluster " << iclnew << " attached to Begin vtx "
                              << ivx;
    } // pos > 2

    if (pos < tcl[icl].tclhits.size() - 3) {
      // have enough hits to make a cluster at the End
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
      for (ii = pos; ii < tcl[icl].tclhits.size(); ++ii) {
        iht = tcl[icl].tclhits[ii];
        if (inClus[iht] != 0) {
          RestoreObsoleteCluster(icl);
          return false;
        }
        fcl2hits.push_back(iht);
        // define the Begin parameters
        if (fcl2hits.size() == fMaxHitsFit[pass] || fcl2hits.size() == fMinHits[pass]) {
          FitCluster();
          clBeginSlp = clpar[1];
          clBeginAng = std::atan(fScaleF * clBeginSlp);
          clBeginSlpErr = clparerr[1];
          didFit = true;
        }
        if ((unsigned short)fcl2hits.size() == fNHitsAve[pass] + 1) {
          FitClusterChg();
          clBeginChg = fAveChg;
          didFit = true;
        }
      } // ii
      // do a fit using all hits if one wasn't done
      if (!didFit) {
        FitCluster();
        FitClusterChg();
        clBeginChg = fAveChg;
      }
      if (!TmpStore()) {
        // clobber the previously stored cluster
        MakeClusterObsolete(tcl.size() - 1);
        RestoreObsoleteCluster(icl);
        return false;
      }
      // associate the End with the supplied vertex
      iclnew = tcl.size() - 1;
      tcl[iclnew].BeginVtx = ivx;
      tcl[iclnew].EndVtx = tcl[icl].EndVtx;
      if (vtxprt)
        mf::LogVerbatim("CC") << "SplitCluster made cluster " << iclnew << " attached to End vtx "
                              << ivx;
    }

    return true;
  } // SplitCluster()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::ChkMerge()
  {
    // Try to merge clusters. Clusters that have been subsumed in other
    // clusters, i.e. no longer valid, have ID < 0

    if (tcl.size() < 2) return;
    // The size of the ClusterStore vector will increase while merging
    // is on-going so the upper limit on it1 is fixed tcl.size() - 1
    // before merging starts

    prt = (fDebugPlane == (short)plane && fDebugWire < 0);
    if (prt) mf::LogVerbatim("CC") << "ChkMerge on pass " << pass;

    unsigned short it1, it2, nh1, pass1, pass2;
    float bs1, bth1, bt1, bc1, es1, eth1, et1, ec1;
    float bs2, bth2, bt2, bc2, es2, eth2, et2, ec2;
    int bw1, ew1, bw2, ew2, ndead;
    float dth, dw, angcut, chgrat, chgcut, dtim, timecut, bigslp;
    bool bothLong, NoVtx;

    int skipcut, tclsize = tcl.size();
    int maxOverlap = 3;

    for (it1 = 0; it1 < tclsize - 1; ++it1) {
      // ignore already merged clusters
      if (tcl[it1].ID < 0) continue;
      if (tcl[it1].CTP != clCTP) continue;
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
      if (pass1 > fNumPass) pass1 = fNumPass;
      for (it2 = it1 + 1; it2 < tclsize; ++it2) {
        // ignore already merged clusters
        if (tcl[it1].ID < 0) continue;
        if (tcl[it2].ID < 0) continue;
        // only merge if they are in the right cryostat/TPC/plane
        if (tcl[it2].CTP != clCTP) continue;
        // Don't bother if these clusters, if merged, would fail the
        // cluster hit fraction cut
        if (!ChkMergedClusterHitFrac(it1, it2)) { continue; }
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
        if (pass2 > fNumPass) pass2 = fNumPass;
        // use the more promiscuous pass for cuts
        angcut = fKinkAngCut[pass1];
        if (fKinkAngCut[pass2] > angcut) angcut = fKinkAngCut[pass2];
        skipcut = fMaxWirSkip[pass1];
        if (fMaxWirSkip[pass2] > skipcut) skipcut = fMaxWirSkip[pass2];
        chgcut = fMergeChgCut[pass1];
        if (fMergeChgCut[pass2] > chgcut) chgcut = fMergeChgCut[pass2];
        // tweak the cuts for long straight tracks
        bothLong = (nh1 > 100 && tcl[it2].tclhits.size() > 100);

        // look for US and DS broken clusters w similar angle.
        // US cluster 2 merge with DS cluster 1?
        // This is the most likely occurrence given the order in which
        // clusters are created so put it first.
        dth = std::abs(bth2 - eth1);
        ndead = DeadWireCount(bw2, ew1);
        dw = ew1 - bw2 - ndead;
        // require no vertex between
        NoVtx = (tcl[it1].EndVtx < 0) && (tcl[it2].BeginVtx < 0);
        if (prt && bw2 < (ew1 + maxOverlap))
          mf::LogVerbatim("CC") << "Chk1 ID1-2 " << tcl[it1].ID << "-" << tcl[it2].ID << " " << ew1
                                << ":" << (int)et1 << " " << bw2 << ":" << (int)bt2 << " dw " << dw
                                << " ndead " << ndead << " skipcut " << skipcut << " dth " << dth
                                << " angcut " << angcut;
        if (NoVtx && bw2 < (ew1 + maxOverlap) && dw < skipcut && dth < angcut) {
          chgrat = 2 * fabs(bc2 - ec1) / (bc2 + ec1);
          // ignore the charge cut for long tracks with small dth
          if (bothLong && dth < 0.05) chgrat = 0.;
          // project bw2,bt2 to ew1
          dtim = fabs(bt2 + (ew1 - bw2) * bs2 - et1);
          bigslp = std::abs(bs2);
          if (std::abs(es1) > bigslp) bigslp = std::abs(es1);
          timecut = fTimeDelta[pass2] * AngleFactor(bigslp);
          if (prt)
            mf::LogVerbatim("CC") << " dtim " << dtim << " timecut " << (int)timecut << " ec1 "
                                  << ec1 << " bc2 " << bc2 << " chgrat " << chgrat << " chgcut "
                                  << chgcut << " es1 " << es1 << " ChkSignal "
                                  << ChkSignal(ew1, et1, bw2, bt2);
          if (chgrat < chgcut && dtim < timecut) {
            // ensure there is a signal between cluster ends
            if (ChkSignal(ew1, et1, bw2, bt2)) {
              DoMerge(it2, it1, 10);
              tclsize = tcl.size();
              break;
            }
          } // chgrat < chgcut ...
        }   // US cluster 2 merge with DS cluster 1?

        // look for US and DS broken clusters w similar angle
        // US cluster 1 merge with DS cluster 2?
        dth = fabs(bth1 - eth2);
        ndead = DeadWireCount(bw1, ew2);
        dw = ew2 - bw1 - ndead;
        if (prt && bw1 < (ew2 + maxOverlap) && dw < skipcut)
          mf::LogVerbatim("CC") << "Chk2 ID1-2 " << tcl[it1].ID << "-" << tcl[it2].ID << " " << bw1
                                << ":" << (int)bt1 << " " << bw2 << ":" << (int)et2 << " dw " << dw
                                << " ndead " << ndead << " skipcut " << skipcut << " dth " << dth
                                << " angcut " << angcut;
        // require no vertex between
        NoVtx = (tcl[it2].EndVtx < 0) && (tcl[it1].BeginVtx < 0);
        if (NoVtx && bw1 < (ew2 + maxOverlap) && dw < skipcut && dth < angcut) {
          chgrat = 2 * fabs((bc1 - ec2) / (bc1 + ec2));
          // ignore the charge cut for long tracks with small dth
          if (bothLong && dth < 0.05) chgrat = 0.;
          // project sw1,st1 to ew2
          dtim = std::abs(bt1 + (ew2 - bw1) * bs1 - et2);
          bigslp = std::abs(bs1);
          if (std::abs(bs2) > bigslp) bigslp = std::abs(bs2);
          timecut = fTimeDelta[pass2] * AngleFactor(bigslp);
          if (prt)
            mf::LogVerbatim("CC") << " dtim " << dtim << " err " << dtim << " timecut "
                                  << (int)timecut << " chgrat " << chgrat << " chgcut " << chgcut
                                  << " ChkSignal " << ChkSignal(bw1, bt1, ew2, et2);
          // TODO: we should be checking for a signal here like we did above
          if (chgrat < chgcut && dtim < timecut) {
            if (ChkSignal(bw1, bt1, ew2, et2)) {
              DoMerge(it1, it2, 10);
              tclsize = tcl.size();
              break;
            }
          } // chgrat < chgcut ...
        }   // US cluster 1 merge with DS cluster 2

        if (bw2 < bw1 && ew2 > ew1) {
          // look for small cl2 within the wire boundary of cl1
          // with similar times and slopes for both clusters
          dth = fabs(eth2 - eth1);
          dtim = fabs(et1 + (ew2 - ew1 - 1) * es1 - et2);
          bigslp = std::abs(es1);
          if (std::abs(es1) > bigslp) bigslp = std::abs(es1);
          timecut = fTimeDelta[pass2] * AngleFactor(bigslp);
          // count the number of wires with no hits on cluster 1
          short nmiss1 = bw1 - ew1 + 1 - tcl[it1].tclhits.size();
          // compare with the number of hits in cluster 2
          short nin2 = tcl[it2].tclhits.size();
          if (prt)
            mf::LogVerbatim("CC") << "cl2: " << ew2 << ":" << (int)et2 << " - " << bw2 << ":"
                                  << (int)bt2 << " within cl1 " << ew1 << ":" << (int)et1 << " - "
                                  << bw1 << ":" << (int)bt1 << " ? dth " << dth << " dtim " << dtim
                                  << " nmissed " << nmiss1 << " timecut " << timecut
                                  << " FIX THIS CODE";
          // make rough cuts before calling ChkMerge12
          // this may not work well for long wandering clusters
          // TODO fix this code
          bool didit = false;
          if (dth < 1 && dtim < timecut && nmiss1 >= nin2) ChkMerge12(it1, it2, didit);
          if (didit) {
            tclsize = tcl.size();
            break;
          } //didit
        }   // small cl2 within the wire boundary of cl1

        if (bw1 < bw2 && ew1 > ew2) {
          // look for small cl1 within the wire boundary of cl2
          // with similar times and slopes for both clusters
          dth = std::abs(eth2 - eth1);
          dtim = std::abs(et2 + (ew1 - ew2 - 1) * es2 - et1);
          bigslp = std::abs(es1);
          if (std::abs(es1) > bigslp) bigslp = std::abs(es1);
          timecut = fTimeDelta[pass2] * AngleFactor(bigslp);
          // count the number of wires with no hits on cluster 2
          short nmiss2 = bw2 - ew2 + 1 - tcl[it2].tclhits.size();
          // compare with the number of hits in cluster 1
          short nin1 = tcl[it1].tclhits.size();
          if (prt)
            mf::LogVerbatim("CC") << "cl1: " << ew1 << ":" << (int)et1 << " - " << bw1 << ":"
                                  << (int)bt1 << " within cl2 " << ew2 << ":" << (int)et2 << " - "
                                  << bw2 << ":" << (int)bt2 << " ? dth " << dth << " dtim " << dtim
                                  << " nmissed " << nmiss2 << " timecut " << timecut
                                  << " FIX THIS CODE";
          // make rough cuts before calling ChkMerge12
          // this may not work well for long wandering clusters
          bool didit = false;
          if (dth < 1 && dtim < timecut && nmiss2 >= nin1) ChkMerge12(it2, it1, didit);
          if (didit) {
            tclsize = tcl.size();
            break;
          } // didit
        }   // small cl1 within the wire boundary of cl2

        if (tcl[it1].ID < 0) break;
      } // cluster 2
      if (tcl[it1].ID < 0) continue;
    } // cluster 1
  }

  /////////////////////////////////////////
  void ClusterCrawlerAlg::ChkMerge12(unsigned short it1, unsigned short it2, bool& didit)
  {
    // Calling routine has done a rough check that cluster it2 is a candidate
    // for merging with cluster it1. The wire range spanned by it2 lies
    // within the wire range of it1 and the clusters are reasonably close
    // together in time.

    // assume that no merging was done
    didit = false;

    if (prt) mf::LogVerbatim("CC") << "ChkMerge12 " << tcl[it1].ID << " " << tcl[it2].ID;

    ClusterStore& cl1 = tcl[it1];
    // fill a vector spanning the length of cluster 1 and filled with the hit
    // time
    int ew1 = tcl[it1].EndWir;
    int bw1 = tcl[it1].BeginWir;
    unsigned int iht, hit;
    int wire;
    std::vector<unsigned int> cl1hits(bw1 + 1 - ew1);
    // put in the hit IDs
    for (iht = 0; iht < cl1.tclhits.size(); ++iht) {
      hit = cl1.tclhits[iht];
      wire = fHits[hit].WireID().Wire;
      if (wire - ew1 < 0 || wire - ew1 > (int)cl1hits.size()) { return; }
      cl1hits[wire - ew1] = hit;
    }
    unsigned int ew2 = tcl[it2].EndWir;
    float et2 = tcl[it2].EndTim;
    // look for the closest wire with a hit on cluster 1
    unsigned int wiron1 = 0;
    // count the number of missing hits
    short nmiss = 0;
    for (wire = ew2 - 1; wire > ew1; --wire) {
      if (cl1hits[wire - ew1] > 0) {
        wiron1 = wire;
        break;
      }
      ++nmiss;
    } // wire
    if (prt) mf::LogVerbatim("CC") << "chk next US wire " << wiron1 << " missed " << nmiss;
    if (wiron1 == 0) return;
    if (nmiss > fMaxWirSkip[pass]) return;

    // compare the wires with hits on cluster 2 with the gap in cluster 1
    // the number of hit wires that fit in the gap
    unsigned int hiton2;
    int wiron2;
    unsigned short nfit = 0;
    for (iht = 0; iht < tcl[it2].tclhits.size(); ++iht) {
      hiton2 = tcl[it2].tclhits[iht];
      wiron2 = fHits[hiton2].WireID().Wire;
      if (wiron2 < ew1 || wiron2 > bw1) return;
      if (cl1hits[wiron2 - ew1] == 0) ++nfit;
    }
    // require complete filling of the gap
    if (nfit < tcl[it2].tclhits.size()) return;

    // decode the pass for both clusters and select the matching cuts
    unsigned short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
    unsigned short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
    unsigned short cpass = pass1;
    // use the tighter cuts
    if (pass2 < pass1) cpass = pass2;

    // ***** Check End of Cluster 2 matching with middle of cluster 1
    if ((int)wiron1 - ew1 < 0) return;
    unsigned int hiton1 = cl1hits[wiron1 - ew1];
    if (hiton1 > fHits.size() - 1) { return; }
    // check the time difference
    float timon1 = fHits[hiton1].PeakTime();
    float dtim = std::abs(et2 + (wiron1 - ew2) * tcl[it2].EndSlp - timon1);
    if (dtim > fTimeDelta[cpass]) return;
    // check the slope difference. First do a local fit on cluster 1 near
    // the matching point
    FitClusterMid(it1, hiton1, 3);
    if (clChisq > 20.) return;
    // fit parameters are now in clpar.
    // check for angle consistency
    float dth = std::abs(std::atan(fScaleF * clpar[1]) - std::atan(fScaleF * tcl[it2].EndSlp));
    if (prt) mf::LogVerbatim("CC") << "US dtheta " << dth << " cut " << fKinkAngCut[cpass];
    if (dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut. fAveChg was calculated in FitClusterMid
    float chgrat = 2 * std::abs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
    if (prt) mf::LogVerbatim("CC") << "US chgrat " << chgrat << " cut " << fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    bool SigOK;
    SigOK = ChkSignal(wiron1, timon1, ew2, et2);
    if (prt) mf::LogVerbatim("CC") << "US SigOK? " << SigOK;
    if (!SigOK) return;

    // ***** Check Begin of Cluster 2 matching with middle of cluster 1
    unsigned int bw2 = tcl[it2].BeginWir;
    float bt2 = tcl[it2].BeginTim;
    nmiss = 0;
    wiron1 = 0;
    for (wire = bw2 + 1; wire < bw1; ++wire) {
      if (cl1hits[wire - ew1] > 0) {
        wiron1 = wire;
        break;
      }
      ++nmiss;
    }
    if (wiron1 == 0) return;
    if (nmiss > fMaxWirSkip[pass]) return;
    // fit this section of cluster 1 with 4 hits starting at the hit on the
    // closest wire and moving DS
    hiton1 = cl1hits[wiron1 - ew1];
    if (hiton1 > fHits.size() - 1) { return; }
    timon1 = fHits[hiton1].PeakTime();
    dtim = std::abs(bt2 - (wiron1 - bw2) * tcl[it2].BeginSlp - timon1);
    if (dtim > fTimeDelta[cpass]) return;
    FitClusterMid(it1, hiton1, -3);
    if (clChisq > 20.) return;
    // check for angle consistency
    dth = std::abs(std::atan(fScaleF * clpar[1]) - std::atan(fScaleF * tcl[it2].BeginSlp));
    if (prt) mf::LogVerbatim("CC") << "DS dtheta " << dth << " cut " << fKinkAngCut[cpass];
    if (dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * std::abs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
    if (prt) mf::LogVerbatim("CC") << "DS chgrat " << chgrat << " cut " << fMergeChgCut[pass];
    // ensure that there is a signal on any missing wires at the US end of 1
    SigOK = ChkSignal(wiron1, timon1, bw2, bt2);
    if (prt) mf::LogVerbatim("CC") << "DS SigOK? " << SigOK;
    if (!SigOK) return;

    if (prt) mf::LogVerbatim("CC") << "Merge em";
    // success. Merge them
    DoMerge(it1, it2, 100);
    didit = true;
  } // ChkMerge12()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::DoMerge(unsigned short it1, unsigned short it2, short inProcCode)
  {
    // Merge clusters.

    ClusterStore& cl1 = tcl[it1];
    ClusterStore& cl2 = tcl[it2];

    if (cl1.tclhits.size() < 3) return;
    if (cl2.tclhits.size() < 3) return;

    unsigned int lowire, hiwire, whsize, ii, iht, indx;
    // do a fit across the boundary between cl1 and cl2 to
    // ensure that they truly should be merged
    unsigned int fithit;
    // Find the low and high wire for both clusters.
    // Assume that cluster 1 is DS
    bool cl1DS = true;
    hiwire = cl1.BeginWir;
    fithit = cl1.tclhits[cl1.tclhits.size() - 2];
    if (cl2.BeginWir > hiwire) {
      hiwire = cl2.BeginWir;
      fithit = cl2.tclhits[cl2.tclhits.size() - 2];
      cl1DS = false;
    }
    lowire = cl1.EndWir;
    if (cl2.EndWir < lowire) lowire = cl2.EndWir;

    // make a vector of wire hits
    whsize = hiwire + 2 - lowire;
    std::vector<int> wirehit(whsize, -1);
    // put in the hit IDs for cluster 2
    for (ii = 0; ii < cl2.tclhits.size(); ++ii) {
      iht = cl2.tclhits[ii];
      indx = fHits[iht].WireID().Wire - lowire;
      wirehit[indx] = iht;
    } // iht
    // now cluster 1
    for (ii = 0; ii < cl1.tclhits.size(); ++ii) {
      iht = cl1.tclhits[ii];
      indx = fHits[iht].WireID().Wire - lowire;
      wirehit[indx] = iht;
    } // iht
    // make the new cluster
    fcl2hits.clear();
    for (std::vector<int>::reverse_iterator rit = wirehit.rbegin(); rit != wirehit.rend(); ++rit) {
      if (*rit < 0) continue;
      fcl2hits.push_back(*rit);
    } // rit

    // fit the 6 hits that are near the merging point
    short nhitfit = 6;
    FitClusterMid(fcl2hits, fithit, nhitfit);
    if (clChisq > 5) return;

    // mark cl1 and cl2 obsolete
    MakeClusterObsolete(it1);
    MakeClusterObsolete(it2);

    short endVtx = 0;
    short begVtx = 0;
    short del1Vtx = -99;
    short del2Vtx = -99;
    if (cl1DS) {
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
    }
    else {
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
    if (!TmpStore()) return;
    unsigned short itnew = tcl.size() - 1;
    if (prt)
      mf::LogVerbatim("CC") << "DoMerge " << cl1.ID << " " << cl2.ID << " -> " << tcl[itnew].ID;
    // stuff the processor code with the current pass
    tcl[itnew].ProcCode = inProcCode + pass;
    // transfer the vertex info
    // delete a vertex between these two?
    if (del1Vtx >= 0 && del1Vtx == del2Vtx) vtx[del1Vtx].NClusters = 0;
    // preserve the vertex assignments
    tcl[itnew].BeginVtx = begVtx;
    tcl[itnew].EndVtx = endVtx;
  } // DoMerge

  /////////////////////////////////////////
  void ClusterCrawlerAlg::PrintVertices()
  {

    mf::LogVerbatim myprt("CC");

    if (vtx3.size() > 0) {
      // print out 3D vertices
      myprt
        << "****** 3D vertices ******************************************__2DVtx_Indx__*******\n";
      myprt
        << "Vtx  Cstat  TPC Proc     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for (unsigned short iv = 0; iv < vtx3.size(); ++iv) {
        myprt << std::right << std::setw(3) << std::fixed << iv << std::setprecision(1);
        myprt << std::right << std::setw(7) << vtx3[iv].CStat;
        myprt << std::right << std::setw(5) << vtx3[iv].TPC;
        myprt << std::right << std::setw(5) << vtx3[iv].ProcCode;
        myprt << std::right << std::setw(8) << vtx3[iv].X;
        myprt << std::right << std::setw(8) << vtx3[iv].Y;
        myprt << std::right << std::setw(8) << vtx3[iv].Z;
        myprt << std::right << std::setw(5) << vtx3[iv].XErr;
        myprt << std::right << std::setw(5) << vtx3[iv].YErr;
        myprt << std::right << std::setw(5) << vtx3[iv].ZErr;
        myprt << std::right << std::setw(5) << vtx3[iv].Ptr2D[0];
        myprt << std::right << std::setw(5) << vtx3[iv].Ptr2D[1];
        myprt << std::right << std::setw(5) << vtx3[iv].Ptr2D[2];
        myprt << std::right << std::setw(5) << vtx3[iv].Wire;
        if (vtx3[iv].Wire < 0) { myprt << "    Matched in all planes"; }
        else {
          myprt << "    Incomplete";
        }
        myprt << "\n";
      }
    } // vtx3.size

    if (vtx.size() > 0) {
      // print out 2D vertices
      myprt << "************ 2D vertices ************\n";
      myprt << "Vtx   CTP    wire     error   tick     error  ChiDOF  NCl  topo  cluster IDs\n";
      for (unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if (fDebugPlane < 3 && fDebugPlane != (int)vtx[iv].CTP) continue;
        myprt << std::right << std::setw(3) << std::fixed << iv << std::setprecision(1);
        myprt << std::right << std::setw(6) << vtx[iv].CTP;
        myprt << std::right << std::setw(8) << vtx[iv].Wire << " +/- ";
        myprt << std::right << std::setw(4) << vtx[iv].WireErr;
        myprt << std::right << std::setw(8) << vtx[iv].Time << " +/- ";
        myprt << std::right << std::setw(4) << vtx[iv].TimeErr;
        myprt << std::right << std::setw(8) << vtx[iv].ChiDOF;
        myprt << std::right << std::setw(5) << vtx[iv].NClusters;
        myprt << std::right << std::setw(6) << vtx[iv].Topo;
        myprt << "    ";
        // display the cluster IDs
        for (unsigned short ii = 0; ii < tcl.size(); ++ii) {
          if (fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
          if (tcl[ii].ID < 0) continue;
          if (tcl[ii].BeginVtx == (short)iv) myprt << std::right << std::setw(4) << tcl[ii].ID;
          if (tcl[ii].EndVtx == (short)iv) myprt << std::right << std::setw(4) << tcl[ii].ID;
        }
        myprt << "\n";
      } // iv
    }   // vtx.size

  } // PrintVertices

  /////////////////////////////////////////
  void ClusterCrawlerAlg::PrintClusters()
  {

    // prints clusters to the screen for code development
    mf::LogVerbatim myprt("CC");

    PrintVertices();

    float aveRMS, aveRes;
    myprt << "*************************************** Clusters "
             "*********************************************************************\n";
    myprt << "  ID CTP nht Stop  Proc  beg_W:T    bAng  bSlp  Err  bChg end_W:T      eAng  eSlp  "
             "Err  eChg  bVx  eVx aveRMS Qual cnt\n";
    for (unsigned short ii = 0; ii < tcl.size(); ++ii) {
      // print clusters in all planes (fDebugPlane = 3) or in a selected plane
      if (fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
      myprt << std::right << std::setw(4) << tcl[ii].ID;
      myprt << std::right << std::setw(3) << tcl[ii].CTP;
      myprt << std::right << std::setw(5) << tcl[ii].tclhits.size();
      myprt << std::right << std::setw(4) << tcl[ii].StopCode;
      myprt << std::right << std::setw(6) << tcl[ii].ProcCode;
      unsigned int iTime = tcl[ii].BeginTim;
      myprt << std::right << std::setw(6) << tcl[ii].BeginWir << ":" << iTime;
      if (iTime < 10) { myprt << "   "; }
      else if (iTime < 100) {
        myprt << "  ";
      }
      else if (iTime < 1000)
        myprt << " ";
      myprt << std::right << std::setw(7) << std::fixed << std::setprecision(2) << tcl[ii].BeginAng;
      if (std::abs(tcl[ii].BeginSlp) < 100) {
        myprt << std::right << std::setw(6) << std::fixed << std::setprecision(2)
              << tcl[ii].BeginSlp;
      }
      else {
        myprt << std::right << std::setw(6) << (int)tcl[ii].BeginSlp;
      }
      myprt << std::right << std::setw(6) << std::fixed << std::setprecision(2)
            << tcl[ii].BeginSlpErr;
      myprt << std::right << std::setw(5) << (int)tcl[ii].BeginChg;
      iTime = tcl[ii].EndTim;
      myprt << std::right << std::setw(6) << tcl[ii].EndWir << ":" << iTime;
      if (iTime < 10) { myprt << "   "; }
      else if (iTime < 100) {
        myprt << "  ";
      }
      else if (iTime < 1000)
        myprt << " ";
      myprt << std::right << std::setw(7) << std::fixed << std::setprecision(2) << tcl[ii].EndAng;
      if (std::abs(tcl[ii].EndSlp) < 100) {
        myprt << std::right << std::setw(6) << std::fixed << std::setprecision(2) << tcl[ii].EndSlp;
      }
      else {
        myprt << std::right << std::setw(6) << (int)tcl[ii].EndSlp;
      }
      myprt << std::right << std::setw(6) << std::fixed << std::setprecision(2)
            << tcl[ii].EndSlpErr;
      myprt << std::right << std::setw(5) << (int)tcl[ii].EndChg;
      myprt << std::right << std::setw(5) << tcl[ii].BeginVtx;
      myprt << std::right << std::setw(5) << tcl[ii].EndVtx;
      aveRMS = 0;
      unsigned int iht = 0;
      for (unsigned short jj = 0; jj < tcl[ii].tclhits.size(); ++jj) {
        iht = tcl[ii].tclhits[jj];
        aveRMS += fHits[iht].RMS();
      }
      aveRMS /= (float)tcl[ii].tclhits.size();
      myprt << std::right << std::setw(5) << std::fixed << std::setprecision(1) << aveRMS;
      aveRes = 0;
      // find cluster tracking resolution
      unsigned int hit0, hit1, hit2, cnt = 0;
      float arg;
      for (unsigned short iht = 1; iht < tcl[ii].tclhits.size() - 1; ++iht) {
        hit1 = tcl[ii].tclhits[iht];
        hit0 = tcl[ii].tclhits[iht - 1];
        hit2 = tcl[ii].tclhits[iht + 1];
        // require hits on adjacent wires
        if (fHits[hit1].WireID().Wire + 1 != fHits[hit0].WireID().Wire) continue;
        if (fHits[hit2].WireID().Wire + 1 != fHits[hit1].WireID().Wire) continue;
        arg = (fHits[hit0].PeakTime() + fHits[hit2].PeakTime()) / 2 - fHits[hit1].PeakTime();
        aveRes += arg * arg;
        ++cnt;
      }
      if (cnt > 1) {
        aveRes /= (float)cnt;
        aveRes = sqrt(aveRes);
        // convert to a quality factor
        aveRes /= (aveRMS * fHitErrFac);
        myprt << std::right << std::setw(6) << std::fixed << std::setprecision(1) << aveRes;
        myprt << std::right << std::setw(5) << std::fixed << cnt;
      }
      else {
        myprt << "    NA";
        myprt << std::right << std::setw(5) << std::fixed << cnt;
      }
      myprt << "\n";
    } // ii

  } // PrintClusters()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::TmpGet(unsigned short it1)
  {
    // copies temp cluster it1 into the fcl2hits vector, etc. This is
    // effectively the inverse of cl2TmpStore

    if (it1 > tcl.size()) return;

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

    if (fcl2hits.size() < 2) return false;
    if (fcl2hits.size() > fHits.size()) return false;

    if (NClusters == SHRT_MAX) return false;

    ++NClusters;

    unsigned int hit0 = fcl2hits[0];
    unsigned int hit;
    unsigned int tCST = fHits[hit0].WireID().Cryostat;
    unsigned int tTPC = fHits[hit0].WireID().TPC;
    unsigned int tPlane = fHits[hit0].WireID().Plane;
    unsigned int lastWire = 0;

    for (unsigned short it = 0; it < fcl2hits.size(); ++it) {
      hit = fcl2hits[it];
      if (inClus[hit] != 0) {
        --NClusters;
        return false;
      }
      // check for WireID() consistency
      if (fHits[hit].WireID().Cryostat != tCST || fHits[hit].WireID().TPC != tTPC ||
          fHits[hit].WireID().Plane != tPlane) {
        --NClusters;
        return false;
      }
      // check for decreasing wire number
      if (clProcCode < 900 && it > 0 && fHits[hit].WireID().Wire >= lastWire) {
        --NClusters;
        return false;
      }
      lastWire = fHits[hit].WireID().Wire;
      inClus[hit] = NClusters;
    }

    // ensure that the cluster begin/end info is correct

    // define the begin/end charge if it wasn't done already
    if (clEndChg <= 0) {
      // use the next to the last two hits. The End hit may have low charge
      unsigned int ih0 = fcl2hits.size() - 2;
      hit = fcl2hits[ih0];
      clEndChg = fHits[hit].Integral();
      hit = fcl2hits[ih0 - 1];
      clEndChg += fHits[hit].Integral();
      clEndChg = clEndChg / 2.;
    }
    if (clBeginChg <= 0) {
      // use the 2nd and third hit. The Begin hit may have low charge
      hit = fcl2hits[1];
      clBeginChg = fHits[hit].Integral();
      hit = fcl2hits[2];
      clBeginChg += fHits[hit].Integral();
      clBeginChg = clBeginChg / 2.;
    }

    hit0 = fcl2hits[0];
    hit = fcl2hits[fcl2hits.size() - 1];

    // store the cluster in the temporary ClusterStore struct
    ClusterStore clstr;

    clstr.ID = NClusters;
    clstr.BeginSlp = clBeginSlp;
    clstr.BeginSlpErr = clBeginSlpErr;
    clstr.BeginAng = std::atan(fScaleF * clBeginSlp);
    clstr.BeginWir = fHits[hit0].WireID().Wire;
    clstr.BeginTim = fHits[hit0].PeakTime();
    clstr.BeginChg = clBeginChg;
    clstr.BeginChgNear = clBeginChgNear;
    clstr.EndSlp = clEndSlp;
    clstr.EndSlpErr = clEndSlpErr;
    clstr.EndAng = std::atan(fScaleF * clEndSlp);
    clstr.EndWir = fHits[hit].WireID().Wire;
    clstr.EndTim = fHits[hit].PeakTime();
    clstr.EndChg = clEndChg;
    clstr.EndChgNear = clEndChgNear;
    clstr.StopCode = clStopCode;
    clstr.ProcCode = clProcCode;
    clstr.BeginVtx = -99;
    clstr.EndVtx = -99;
    clstr.CTP = EncodeCTP(tCST, tTPC, tPlane);
    clstr.tclhits = fcl2hits;
    tcl.push_back(clstr);

    return true;
  } // TmpStore()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::LACrawlUS()
  {
    // Crawl a large angle cluster upstream. Similar to CrawlUS but require
    // that a hit be added on each wire

    unsigned int dhit = fcl2hits[0];
    short dwir = fHits[dhit].WireID().Wire;
    clLA = true;
    prt = false;
    if (fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
      prt = std::abs(fHits[dhit].PeakTime() - fDebugHit) < 40;

    if (prt) {
      mf::LogVerbatim myprt("CC");
      myprt << "******************* LACrawlUS PASS " << pass << " Hits ";
      for (unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
        unsigned int iht = fcl2hits[fcl2hits.size() - 1 - ii];
        myprt << fHits[iht].WireID().Wire << ":" << (int)fHits[iht].PeakTime() << " ";
      }
    }

    bool SigOK = true;
    bool HitOK = true;
    short nmissed = 0;
    // count the number of kinks encountered. Hits US of the kink are removed
    // and crawling continues unless another kink is encountered
    unsigned short kinkOnWire = 0;
    unsigned int it = fcl2hits.size() - 1;
    unsigned int lasthit = fcl2hits[it];
    unsigned int lastwire = fHits[lasthit].WireID().Wire;
    unsigned int prevHit, prevWire;
    bool ChkCharge = false;
    for (unsigned int nextwire = lastwire - 1; nextwire >= fFirstWire; --nextwire) {
      if (prt)
        mf::LogVerbatim("CC") << "LACrawlUS: next wire " << nextwire << " HitRange "
                              << WireHitRange[nextwire].first;
      // stop crawling if there is a nearby vertex
      if (CrawlVtxChk(nextwire)) {
        if (prt) mf::LogVerbatim("CC") << "LACrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // AddLAHit will merge the hit on nextwire if necessary
      AddLAHit(nextwire, ChkCharge, HitOK, SigOK);
      if (prt)
        mf::LogVerbatim("CC") << "LACrawlUS: HitOK " << HitOK << " SigOK " << SigOK
                              << " nmissed SigOK " << nmissed << " cut " << fAllowNoHitWire;
      if (SigOK) nmissed = 0;
      if (!SigOK) {
        ++nmissed;
        if (nmissed > fAllowNoHitWire) {
          clStopCode = 1;
          break;
        }
        continue;
      }
      // If a hit was added after a gap, check to see if there is indeed
      // a wire signal in the gap
      if (HitOK) {
        prevHit = fcl2hits[fcl2hits.size() - 2];
        prevWire = fHits[prevHit].WireID().Wire;
        if (prevWire > nextwire + 2) {
          if (!ChkSignal(fcl2hits[fcl2hits.size() - 1], fcl2hits[fcl2hits.size() - 2])) {
            // no hit so trim the last hit and quit
            FclTrimUS(1);
            break;
          } // no signal
        }   // prevWire > nextwire + 2
      }     // HitOK
      // Merge all of the hit multiplets in the fcl2hits array into single
      // hits when enough hits have been found to call this a credible large
      // angle cluster. The last hit was already merged in AddHit
      if (fcl2hits.size() == 4) {
        bool didMerge;
        for (unsigned short kk = 0; kk < fcl2hits.size() - 1; ++kk) {
          unsigned int hit = fcl2hits[kk];
          if (mergeAvailable[hit]) MergeHits(hit, didMerge);
        }
        // update the fit
        FitCluster();
        clBeginSlp = clpar[1];
        // start checking the charge ratio when adding new hits
        ChkCharge = true;
        continue;
      } // fcl2hits.size() == 4
      unsigned short chsiz = chifits.size() - 1;
      // chsiz is fcl2hits.size() - 1...
      if (chsiz < 6) continue;
      if (fKinkChiRat[pass] <= 0) continue;
      if (chifits.size() != fcl2hits.size()) {
        mf::LogError("CC") << "LACrawlUS: chifits size error " << chifits.size() << " "
                           << fcl2hits.size();
        return;
      }
      if (prt)
        mf::LogVerbatim("CC") << "Kink chk " << chifits[chsiz] << " " << chifits[chsiz - 1] << " "
                              << chifits[chsiz - 2] << " " << chifits[chsiz - 3];
      if (chifits[chsiz - 1] > fKinkChiRat[pass] * chifits[chsiz - 2] &&
          chifits[chsiz] > fKinkChiRat[pass] * chifits[chsiz - 1]) {
        // find the kink angle (crudely) from the 0th and 2nd hit
        unsigned int ih0 = fcl2hits.size() - 1;
        unsigned int hit0 = fcl2hits[ih0];
        unsigned int ih2 = ih0 - 2;
        unsigned int hit2 = fcl2hits[ih2];
        float dt02 = fHits[hit2].PeakTime() - fHits[hit0].PeakTime();
        float dw02 = fHits[hit2].WireID().Wire - fHits[hit0].WireID().Wire;
        float th02 = std::atan(fScaleF * dt02 / dw02);
        // and the 3rd and 5th hit
        unsigned int ih3 = ih2 - 1;
        unsigned int hit3 = fcl2hits[ih3];
        unsigned int ih5 = ih3 - 2;
        unsigned int hit5 = fcl2hits[ih5];
        float dt35 = fHits[hit5].PeakTime() - fHits[hit3].PeakTime();
        float dw35 = fHits[hit5].WireID().Wire - fHits[hit3].WireID().Wire;
        float th35 = std::atan(fScaleF * dt35 / dw35);
        float dth = std::abs(th02 - th35);
        if (prt)
          mf::LogVerbatim("CC") << " Kink angle " << std::setprecision(3) << dth << " cut "
                                << fKinkAngCut[pass];
        if (dth > fKinkAngCut[pass]) {
          // hit a kink. Lop of the first 3 hits, refit and keep crawling?
          FclTrimUS(3);
          FitCluster();
          // See if this is a second kink and it is close to the first
          // kink (which had hits removed).
          if (kinkOnWire > 0) {
            if (kinkOnWire - nextwire < 4) {
              if (prt)
                mf::LogVerbatim("CC")
                  << "Hit a second kink. kinkOnWire = " << kinkOnWire << " Stopping";
              // set the kink stop code
              clStopCode = 3;
              break;
            }
          }
          kinkOnWire = nextwire;
          if (prt) mf::LogVerbatim("CC") << "Removed kink hits";
        } // kinkang check
      }   // chifits test
    }     // nextwire

    CheckClusterHitFrac(prt);

    clProcCode += 300;
    if (prt) mf::LogVerbatim("CC") << "LACrawlUS done. Nhits = " << fcl2hits.size();
    prt = false;
  } // LACrawlUS

  /////////////////////////////////////////
  void ClusterCrawlerAlg::CrawlUS()
  {
    // Crawl along a trail of hits moving upstream

    if (fcl2hits.size() < 2) return;

    unsigned int dhit = fcl2hits[0];
    int dwir = fHits[dhit].WireID().Wire;
    clLA = false;

    prt = false;
    if (fDebugPlane == (short)plane && dwir == fDebugWire && fDebugHit > 0)
      prt = std::abs(fHits[dhit].PeakTime() - fDebugHit) < 20;

    if (prt) {
      mf::LogVerbatim myprt("CC");
      myprt << "******************* Start CrawlUS on pass " << pass << " Hits: ";
      for (unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
        unsigned int iht = fcl2hits[fcl2hits.size() - 1 - ii];
        myprt << fHits[iht].WireID().Wire << ":" << (int)fHits[iht].PeakTime() << " ";
      }
      myprt << "\n";
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
    if (lasthit > fHits.size() - 1) {
      mf::LogError("CC") << "CrawlUS bad lasthit " << lasthit;
      return;
    }
    unsigned int lastwire = fHits[lasthit].WireID().Wire;
    if (prt) mf::LogVerbatim("CC") << "CrawlUS: last wire " << lastwire << " hit " << lasthit;

    unsigned int lastWireWithSignal = lastwire;
    unsigned short nDroppedHit = 0;

    for (unsigned int nextwire = lastwire - 1; nextwire > fFirstWire; --nextwire) {
      if (prt)
        mf::LogVerbatim("CC") << "CrawlUS: next wire " << nextwire << " HitRange "
                              << WireHitRange[nextwire].first;
      // stop crawling if there is a nearby vertex
      if (CrawlVtxChk(nextwire)) {
        if (prt) mf::LogVerbatim("CC") << "CrawlUS: stop at vertex";
        clStopCode = 6;
        break;
      }
      // Switch to large angle crawling?
      if (std::abs(clpar[1]) > fLAClusSlopeCut) {
        if (prt) mf::LogVerbatim("CC") << ">>>>> CrawlUS: Switching to LACrawlUS";
        LACrawlUS();
        return;
      }
      // add hits and check for PH and width consistency
      AddHit(nextwire, HitOK, SigOK);
      if (prt)
        mf::LogVerbatim("CC") << "CrawlUS: HitOK " << HitOK << " SigOK " << SigOK << " nmissed "
                              << nmissed;
      if (SigOK) lastWireWithSignal = nextwire;
      if (!HitOK) {
        // no hit on this wire. Was there a signal or dead wire?
        if (SigOK) {
          // no hit on the wire but there is a signal
          ++nmissed;
          // see if we are in the PostSkip phase and missed more than 1 wire
          if (PostSkip && nmissed > fMinWirAfterSkip[pass]) {
            // cluster is really short
            if ((int)(fcl2hits.size() - nHitAfterSkip) < 4) {
              fcl2hits.clear();
              return;
            }
            if (prt) mf::LogVerbatim("CC") << " PostSkip && nmissed = " << nmissed;
            clStopCode = 2;
            FclTrimUS(nHitAfterSkip);
            FitCluster();
            if (clChisq > 90.) {
              fcl2hits.clear();
              return;
            }
            FitCluster();
            return;
          } // PostSkip && nmissed >
          if (nmissed > 1) {
            DidaSkip = true;
            PostSkip = false;
          }
        } // SigOK
        else {
          // SigOK is false
          clStopCode = 0;
          lasthit = fcl2hits[fcl2hits.size() - 1];
          if ((lastWireWithSignal - nextwire) > fAllowNoHitWire) {
            if (prt)
              mf::LogVerbatim("CC")
                << "No hit or signal on wire " << nextwire << " last wire with signal "
                << lastWireWithSignal << " exceeding fAllowNoHitWire " << fAllowNoHitWire
                << " Break!";
            break;
          }
        } // else SigOK false
      }   // !HitOK
      else {
        if (prt)
          mf::LogVerbatim("CC") << " CrawlUS check clChisq " << clChisq << " cut " << fChiCut[pass];
        if (clChisq > fChiCut[pass]) {
          if (fcl2hits.size() < 3) {
            fcl2hits.clear();
            return;
          }
          // a fit error occurred. Lop off the leading hit and refit
          if (prt) mf::LogVerbatim("CC") << "ADD- Bad clChisq " << clChisq << " dropping hit";
          FclTrimUS(1);
          FitCluster();
          ++nDroppedHit;
          if (nDroppedHit > 4) {
            if (prt)
              mf::LogVerbatim("CC")
                << " Too many dropped hits: " << nDroppedHit << " Stop crawling";
            break;
          } // too many dropped hits
          if (clChisq > fChiCut[pass]) {
            if (prt)
              mf::LogVerbatim("CC")
                << " Bad clChisq " << clChisq << " after dropping hit. Stop crawling";
            break;
          }
          FitClusterChg();
          continue;
        } // clChisq > fChiCut[0]
        // monitor the onset of a kink. Look for a progressive increase
        // in chisq for the previous 0 - 2 hits.
        if (chifits.size() > 5 && fKinkChiRat[pass] > 0) {
          if (chifits.size() != fcl2hits.size()) {
            mf::LogError("CC") << "CrawlUS: chifits size error " << chifits.size() << " "
                               << fcl2hits.size();
            return;
          }
          unsigned short chsiz = chifits.size() - 1;
          if (prt)
            mf::LogVerbatim("CC") << "Kink chk " << chifits[chsiz] << " " << chifits[chsiz - 1]
                                  << " " << chifits[chsiz - 2] << " " << chifits[chsiz - 3];
          if (chifits[chsiz - 1] > fKinkChiRat[pass] * chifits[chsiz - 2] &&
              chifits[chsiz] > fKinkChiRat[pass] * chifits[chsiz - 1]) {
            if (fcl2hits.size() != chifits.size()) {
              mf::LogError("CC") << "bad kink check size " << chifits.size() << " "
                                 << fcl2hits.size() << " plane " << plane << " cluster " << dwir
                                 << ":" << dhit;
              continue;
            }
            if (EndKinkAngle() > fKinkAngCut[pass]) {
              if (prt)
                mf::LogVerbatim("CC")
                  << "******************* Stopped tracking - kink angle " << EndKinkAngle() << " > "
                  << fKinkAngCut[pass] << " Removing 3 hits";
              // kill the last 3 hits and refit
              FclTrimUS(3);
              FitCluster();
              FitClusterChg();
            } // kinkang check
          }   // chifits check
        }     // chifits.size() > 5
        // done with kink check
        // update the cluster Begin information?
        if (fcl2hits.size() == fMaxHitsFit[pass] || fcl2hits.size() == fMinHits[pass]) {
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
        }
        // define the begin cluster charge if it's not defined yet
        if (clBeginChg <= 0 && fAveChg > 0) {
          clBeginChg = fAveChg;
          if (prt) mf::LogVerbatim("CC") << " Set clBeginChg " << clBeginChg;
        }
        // reset nmissed
        nmissed = 0;
        // start counting hits added after skipping
        if (DidaSkip) {
          // start PostSkip phase
          PostSkip = true;
          DidaSkip = false;
          nHitAfterSkip = 0;
        } // DidaSkip
        // check for PostSkip phase
        if (PostSkip) {
          // end the PostSkip phase if there are enough hits
          ++nHitAfterSkip;
          if (nHitAfterSkip == fMinWirAfterSkip[pass]) PostSkip = false;
        }
        // check for bad chisq
        if (clChisq > fChiCut[pass]) {
          if (prt) mf::LogVerbatim("CC") << "<<ADD- Bad chisq " << clChisq;
          // remove the last few hits if there is a systematic increase in chisq and re-fit
          float chirat;
          unsigned short lopped = 0;
          for (unsigned short nlop = 0; nlop < 4; ++nlop) {
            unsigned short cfsize = chifits.size() - 1;
            chirat = chifits[cfsize] / chifits[cfsize - 1];
            if (prt)
              mf::LogVerbatim("CC")
                << "chirat " << chirat << " last hit " << fcl2hits[fcl2hits.size() - 1];
            if (chirat < 1.2) break;
            if (prt) mf::LogVerbatim("CC") << "<<ADD- Bad chisq. Bad chirat " << chirat;
            FclTrimUS(1);
            ++lopped;
            if (fcl2hits.size() < 6) break;
            if (chifits.size() < 6) break;
          } // nlop
          if (fcl2hits.size() < 6) {
            clStopCode = 4;
            if (prt) mf::LogVerbatim("CC") << "Bad fit chisq - short cluster. Break";
            break;
          }
          if (lopped == 0 && fcl2hits.size() > 5) {
            if (prt) mf::LogVerbatim("CC") << "<<ADD- Bad chisq. remove 1 hit";
            FclTrimUS(1);
            ++lopped;
          }
          FitCluster();
          FitClusterChg();
          if (prt)
            mf::LogVerbatim("CC") << "Bad fit chisq - removed " << lopped << " hits. Continue";
        } // clChisq > fChiCut[pass]
      }   // !HitOK check
    }     // nextwire
    if (prt)
      mf::LogVerbatim("CC") << "******************* CrawlUS normal stop. size " << fcl2hits.size();

    bool reFit = false;
    // end kink angle check
    if (fcl2hits.size() > 5) {
      // check for a kink at the US end
      if (prt)
        mf::LogVerbatim("CC") << "EndKinkAngle check " << EndKinkAngle() << " cut "
                              << fKinkAngCut[pass];
      if (EndKinkAngle() > fKinkAngCut[pass]) {
        if (prt) mf::LogVerbatim("CC") << "EndKinkAngle removes 3 hits ";
        FclTrimUS(3);
        reFit = true;
      }
    } // fcl2hits.size() > 5

    // count the number of hits on adjacent wires at the leading edge and
    // ensure that the count is consistent with fMinWirAfterSkip
    if ((unsigned short)fcl2hits.size() > fMinWirAfterSkip[pass] + 3) {
      unsigned int ih0 = fcl2hits.size() - 1;
      unsigned int hit0 = fcl2hits[ih0];
      unsigned int uswir = fHits[hit0].WireID().Wire;
      unsigned int nxtwir;
      unsigned short nAdjHit = 0;
      for (unsigned short ii = ih0 - 1; ii > 0; --ii) {
        nxtwir = fHits[fcl2hits[ii]].WireID().Wire;
        ++nAdjHit;
        if (nxtwir != uswir + 1) break;
        // break if there are enough hits
        if (nAdjHit == fMinWirAfterSkip[pass]) break;
        uswir = nxtwir;
      } // ii
      // lop off hits?
      if (nAdjHit < fMinWirAfterSkip[pass]) {
        if (prt) mf::LogVerbatim("CC") << "fMinWirAfterSkip removes " << nAdjHit << " hits ";
        FclTrimUS(nAdjHit);
        reFit = true;
      }
    } // fcl2hits.size() > fMinWirAfterSkip[pass] + 3

    // check for a bad hit on the end
    if (!reFit && fcl2hits.size() > 3) {
      float chirat = chifits[chifits.size() - 1] / chifits[chifits.size() - 2];
      if (prt)
        mf::LogVerbatim("CC") << "Last hit chirat " << chirat << " cut " << fKinkChiRat[pass];
      if (prt)
        mf::LogVerbatim("CC") << "Check " << clChisq << " " << chifits[chifits.size() - 1] << " "
                              << chifits[chifits.size() - 2];
      if (chirat > fKinkChiRat[pass]) {
        if (prt) mf::LogVerbatim("CC") << "<<ADD- at end";
        FclTrimUS(1);
        reFit = true;
      }
    } // !reFit

    if (reFit) {
      FitCluster();
      FitClusterChg();
    }
    CheckClusterHitFrac(prt);
    if (prt)
      mf::LogVerbatim("CC") << "******************* CrawlUS done. Size " << fcl2hits.size()
                            << " min length for this pass " << fMinHits[pass];

    prt = false;
  } // CrawlUS()

  /////////////////////////////////////////
  float ClusterCrawlerAlg::EndKinkAngle()
  {
    // find the kink angle (crudely) from the 0th and 2nd hit on the cluster under construction

    unsigned int ih0 = fcl2hits.size() - 1;
    unsigned int hit0 = fcl2hits[ih0];
    unsigned int ih2 = ih0 - 2;
    unsigned int hit2 = fcl2hits[ih2];
    float dt02 = fHits[hit2].PeakTime() - fHits[hit0].PeakTime();
    float dw02 = fHits[hit2].WireID().Wire - fHits[hit0].WireID().Wire;
    float th02 = std::atan(fScaleF * dt02 / dw02);
    // and the 3rd and 5th hit
    unsigned int ih3 = ih2 - 1;
    unsigned int hit3 = fcl2hits[ih3];
    unsigned int ih5 = ih3 - 2;
    unsigned int hit5 = fcl2hits[ih5];
    float dt35 = fHits[hit5].PeakTime() - fHits[hit3].PeakTime();
    float dw35 = fHits[hit5].WireID().Wire - fHits[hit3].WireID().Wire;
    float th35 = std::atan(fScaleF * dt35 / dw35);
    return std::abs(th02 - th35);
  }

  /////////////////////////////////////////
  bool ClusterCrawlerAlg::ChkMergedClusterHitFrac(unsigned short it1, unsigned short it2)
  {
    // This routine is called before two tcl clusters, it1 and it2, are slated to be
    // merged to see if they will satisfy the minimum hit fraction criterion
    if (it1 > tcl.size() - 1) return false;
    if (it2 > tcl.size() - 1) return false;
    unsigned int eWire = 99999;
    unsigned int bWire = 0, wire;
    if (tcl[it1].BeginWir > bWire) bWire = tcl[it1].BeginWir;
    if (tcl[it2].BeginWir > bWire) bWire = tcl[it2].BeginWir;
    if (tcl[it1].EndWir < eWire) eWire = tcl[it1].EndWir;
    if (tcl[it2].EndWir < eWire) eWire = tcl[it2].EndWir;
    unsigned int mergedClusterLength = bWire - eWire + 1;
    // define a vector of size = length of the wire range
    std::vector<bool> cHits(mergedClusterLength, false);
    // set the elements true if there is a hit
    unsigned int ii, iht, indx;
    for (ii = 0; ii < tcl[it1].tclhits.size(); ++ii) {
      iht = tcl[it1].tclhits[ii];
      if (iht > fHits.size() - 1) {
        mf::LogError("CC") << "ChkMergedClusterHitFrac bad iht ";
        return false;
      }
      indx = fHits[iht].WireID().Wire - eWire;
      cHits[indx] = true;
    } // ii
    for (ii = 0; ii < tcl[it2].tclhits.size(); ++ii) {
      iht = tcl[it2].tclhits[ii];
      if (iht > fHits.size() - 1) {
        mf::LogError("CC") << "ChkMergedClusterHitFrac bad iht ";
        return false;
      }
      indx = fHits[iht].WireID().Wire - eWire;
      cHits[indx] = true;
    } // ii
    // set cHits true if the wire is dead
    for (ii = 0; ii < cHits.size(); ++ii) {
      wire = eWire + ii;
      if (WireHitRange[wire].first == -1) cHits[ii] = true;
    }
    // count the number of wires with hits
    float nhits = std::count(cHits.begin(), cHits.end(), true);
    float hitFrac = nhits / (float)cHits.size();
    return (hitFrac > fMinHitFrac);
  } // ChkMergedClusterHitFrac

  /////////////////////////////////////////
  void ClusterCrawlerAlg::CheckClusterHitFrac(bool prt)
  {

    // Find the fraction of the wires on the cluster that have
    // hits.
    unsigned int iht = fcl2hits[fcl2hits.size() - 1];
    clEndWir = fHits[iht].WireID().Wire;
    clBeginWir = fHits[fcl2hits[0]].WireID().Wire;
    float hitFrac = (float)(fcl2hits.size() + DeadWireCount()) / (float)(clBeginWir - clEndWir + 1);

    if (hitFrac < fMinHitFrac) {
      if (prt)
        mf::LogVerbatim("CC") << "CheckClusterHitFrac: Poor hit fraction " << hitFrac
                              << " clBeginWir " << clBeginWir << " clEndWir " << clEndWir
                              << " size " << fcl2hits.size() << " DeadWireCount "
                              << DeadWireCount();
      fcl2hits.clear();
      return;
    } // hitFrac

    /* TODO: Does this make sense?
    // lop off the last hit if it is part of a hit multiplet
    if(fHits[iht].Multiplicity() > 1) {
      fcl2hits.resize(fcl2hits.size() - 1);
    }
*/
    // check for short track ghosts
    if (fcl2hits.size() < 5) {
      unsigned short nsing = 0;
      for (unsigned short iht = 0; iht < fcl2hits.size(); ++iht)
        if (fHits[fcl2hits[iht]].Multiplicity() == 1) ++nsing;
      hitFrac = (float)nsing / (float)fcl2hits.size();
      if (hitFrac < fMinHitFrac) {
        fcl2hits.clear();
        if (prt)
          mf::LogVerbatim("CC") << "CheckClusterHitFrac: Poor short track hit fraction " << hitFrac;
        return;
      } // hitFrac
    }   // short ghost track check

    // analyze the pattern of nearby charge
    // will need the cluster charge so calculate it here if it isn't defined yet
    if (clBeginChg <= 0) {
      unsigned int iht, nht = 0;
      for (unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
        iht = fcl2hits[ii];
        clBeginChg += fHits[iht].Integral();
        ++nht;
        if (nht == 8) break;
      }
      clBeginChg /= (float)nht;
    } // clBeginChg == 0
    // handle short vs long clusters
    unsigned short cnt = chgNear.size() / 2;
    // get the average charge from <= 30 hits at each end
    if (chgNear.size() > 60) cnt = 30;
    clBeginChgNear = 0;
    clEndChgNear = 0;
    for (unsigned short ids = 0; ids < cnt; ++ids) {
      clBeginChgNear += chgNear[ids];
      clEndChgNear += chgNear[chgNear.size() - 1 - ids];
    }
    clBeginChgNear /= (float)cnt;
    clEndChgNear /= (float)cnt;

  } // CheckClusterHitFrac()

  /////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterMid(unsigned short it1, unsigned int ihtin, short nhit)
  {
    FitClusterMid(tcl[it1].tclhits, ihtin, nhit);
  } // FitClusterMid

  /////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterMid(std::vector<unsigned int>& hitVec,
                                        unsigned int ihtin,
                                        short nhit)
  {
    // Fits hits on temp cluster it1 to a line starting at hit ihtin and including
    // nhit hits incrementing towards the hit vector End when nhit > 0 and
    // decrementing towards the hit vector Begin when nhit < 0.
    // The fit params are stashed in the clpar and clparerr arrays.
    // fAveChg is re-calculated as well.

    // set chisq bad in case something doesn't work out
    clChisq = 99.;

    if (hitVec.size() < 3) return;

    std::vector<float> xwir;
    std::vector<float> ytim;
    std::vector<float> ytimerr2;

    unsigned short ii, hitcnt = 0, nht = 0, usnhit;
    float wire0 = 0;
    unsigned int iht;
    bool UseEm = false;
    fAveChg = 0.;
    fChgSlp = 0.;

    if (nhit > 0) {
      usnhit = nhit;
      // find the first desired hit and move towards the End
      for (ii = 0; ii < hitVec.size(); ++ii) {
        iht = hitVec[ii];
        if (iht > fHits.size() - 1) {
          mf::LogError("CC") << "FitClusterMid bad iht " << iht;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if (iht == ihtin) {
          UseEm = true;
          wire0 = fHits[iht].WireID().Wire;
        }
        // use hits after finding the first desired hit
        if (UseEm) {
          xwir.push_back((float)fHits[iht].WireID().Wire - wire0);
          ytim.push_back(fHits[iht].PeakTime());
          // pass the error^2 to the fitter
          float terr = fHitErrFac * fHits[iht].RMS();
          ytimerr2.push_back(terr * terr);
          fAveChg += fHits[iht].Integral();
          ++hitcnt;
          if (hitcnt == usnhit) break;
        }
      }
      nht = hitcnt;
    }
    else {
      usnhit = -nhit;
      // find the first desired hit and move towards the Begin
      for (auto ii = hitVec.crbegin(); ii != hitVec.crend(); ++ii) {
        iht = *ii;
        if (iht > fHits.size() - 1) {
          mf::LogVerbatim("CC") << "FitClusterMid bad iht " << iht;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if (iht == ihtin) {
          UseEm = true;
          wire0 = fHits[iht].WireID().Wire;
        }
        // use hits after finding the first desired hit
        if (UseEm) {
          xwir.push_back((float)fHits[iht].WireID().Wire - wire0);
          ytim.push_back(fHits[iht].PeakTime());
          float terr = fHitErrFac * fHits[iht].RMS();
          ytimerr2.push_back(terr * terr);
          fAveChg += fHits[iht].Integral();
          ++hitcnt;
          if (hitcnt == usnhit) break;
        }
      }
      nht = hitcnt;
    }

    if (nht < 2) return;
    fAveChg = fAveChg / (float)nht;
    fChgSlp = 0.;

    float intcpt = 0.;
    float slope = 0.;
    float intcpterr = 0.;
    float slopeerr = 0.;
    float chidof = 0.;
    fLinFitAlg.LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    clChisq = chidof;
    if (clChisq > fChiCut[0]) return;
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

    if (pass > fNumPass - 1) {
      mf::LogError("CC") << "FitCluster called on invalid pass " << pass;
      return;
    }

    unsigned short ii, nht = 0;
    // fit all hits or truncate?
    nht = fcl2hits.size();
    if (clLA) {
      // Fit large angle cluster
      if (nht > fLAClusMaxHitsFit) nht = fLAClusMaxHitsFit;
    }
    else {
      if (nht > fMaxHitsFit[pass]) nht = fMaxHitsFit[pass];
    }
    if (nht < 2) return;

    std::vector<float> xwir;
    std::vector<float> ytim;
    std::vector<float> ytimerr2;
    // apply an angle dependent scale factor.
    float angfactor = AngleFactor(clpar[1]);

    unsigned int wire;
    unsigned int wire0 = fHits[fcl2hits[fcl2hits.size() - 1]].WireID().Wire;
    unsigned int ihit;
    float terr, qave = 0, qwt;
    for (ii = 0; ii < nht; ++ii) {
      ihit = fcl2hits[fcl2hits.size() - 1 - ii];
      qave += fHits[ihit].Integral();
    } // ii
    qave /= (float)nht;
    for (ii = 0; ii < nht; ++ii) {
      ihit = fcl2hits[fcl2hits.size() - 1 - ii];
      wire = fHits[ihit].WireID().Wire;
      xwir.push_back(wire - wire0);
      ytim.push_back(fHits[ihit].PeakTime());
      // Scale error by hit multiplicity to account for bias in hit
      // multiplet fitting
      terr = fHitErrFac * fHits[ihit].RMS() * fHits[ihit].Multiplicity();
      if (fAveChg > 0) {
        // increase the error for large charge hits
        qwt = (fHits[ihit].Integral() / qave) - 1;
        if (qwt < 1) qwt = 1;
        terr *= qwt;
      }
      if (terr < 0.01) terr = 0.01;
      ytimerr2.push_back(angfactor * terr * terr);
    }
    CalculateAveHitWidth();
    if (prt) {
      mf::LogVerbatim myprt("CC");
      myprt << "FitCluster W:T ";
      unsigned short cnt = 0;
      for (std::vector<unsigned int>::reverse_iterator it = fcl2hits.rbegin();
           it != fcl2hits.rend();
           ++it) {
        unsigned int ihit = *it;
        unsigned short wire = fHits[ihit].WireID().Wire;
        myprt << wire << ":" << (short)fHits[ihit].PeakTime() << " ";
        ++cnt;
        if (cnt == 8) {
          myprt << " .... ";
          break;
        }
      }
    } // prt

    if (xwir.size() < 2) return;

    float intcpt = 0.;
    float slope = 0.;
    float intcpterr = 0.;
    float slopeerr = 0.;
    float chidof = 0.;
    fLinFitAlg.LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    clChisq = chidof;
    if (chidof > fChiCut[0]) return;
    clpar[0] = intcpt;
    clpar[1] = slope;
    clpar[2] = wire0;
    clparerr[0] = intcpterr;
    clparerr[1] = slopeerr;

    if (prt)
      mf::LogVerbatim("CC") << "nht " << nht << " fitpar " << (int)clpar[0] << "+/-"
                            << (int)intcpterr << " " << clpar[1] << "+/-" << slopeerr << " clChisq "
                            << clChisq;
  }
  /////////////////////////////////////////
  float ClusterCrawlerAlg::AngleFactor(float slope)
  {
    // returns an angle dependent cluster projection error factor for fitting
    // and hit finding

    float slp = std::abs(slope);
    if (slp > 15) slp = 15;
    // return a value between 1 and 4
    float angfac = 1 + 0.03 * slp * slp;
    return angfac;
  }

  /////////////////////////////////////////
  void ClusterCrawlerAlg::CalculateAveHitWidth()
  {
    fAveHitWidth = 0;
    for (unsigned short ii = 0; ii < fcl2hits.size(); ++ii)
      fAveHitWidth += fHits[fcl2hits[ii]].EndTick() - fHits[fcl2hits[ii]].StartTick();
    fAveHitWidth /= (float)fcl2hits.size();
  } // CalculateAveHitWidth

  /////////////////////////////////////////
  void ClusterCrawlerAlg::FitClusterChg()
  {
    // Fits the charge of hits on the fcl2hits vector to a line, or simply
    // uses the average of 1 or 2 hits as determined by NHitsAve

    if (fcl2hits.size() == 0) return;
    unsigned int ih0 = fcl2hits.size() - 1;

    if (pass >= fNumPass) {
      mf::LogError("CC") << "FitClusterChg bad pass " << pass;
      return;
    }

    // Handle Large Angle clusters
    if (clLA) {
      // simple average of the charge (and the hit width
      // while we are here)
      unsigned short nhave = fLAClusMaxHitsFit;
      if (nhave > fcl2hits.size()) nhave = fcl2hits.size();
      fAveChg = 0;
      fChgSlp = 0;
      fAveHitWidth = 0;
      unsigned int iht;
      for (unsigned short ii = 0; ii < nhave; ++ii) {
        iht = fcl2hits[fcl2hits.size() - 1 - ii];
        fAveChg += fHits[iht].Integral();
        fAveHitWidth += (fHits[iht].EndTick() - fHits[iht].StartTick());
      } // ii
      fAveChg /= (float)fcl2hits.size();
      fAveHitWidth /= (float)fcl2hits.size();
      return;
    } // clLA

    // number of hits at the leading edge that we will fit
    unsigned short fitLen = fNHitsAve[pass];
    // start fitting charge when there are at least 6 hits if we are tracking
    // long clusters
    if (fitLen > 5 &&             // Fit 6 hits when tracking long clusters AND
        fcl2hits.size() > 5 &&    // there are at least 6 hits AND
        fcl2hits.size() < fitLen) // there are less than fNHitsAve[pass]
      fitLen = 5;

    // don't find the average charge --> no charge cut is made
    if (fNHitsAve[pass] < 1) return;

    if (fNHitsAve[pass] == 1) {
      // simply use the charge and width the last hit
      fAveChg = fHits[fcl2hits[ih0]].Integral();
      fChgSlp = 0.;
    }
    else if (fNHitsAve[pass] == 2) {
      // average the last two points if requested
      fAveChg = (fHits[fcl2hits[ih0]].Integral() + fHits[fcl2hits[ih0 - 1]].Integral()) / 2.;
      fChgSlp = 0.;
    }
    else if ((unsigned short)fcl2hits.size() > fitLen) {
      // do a real fit
      std::vector<float> xwir;
      std::vector<float> ychg;
      std::vector<float> ychgerr2;
      // origin of the fit
      unsigned int wire0 = fHits[fcl2hits[fcl2hits.size() - 1]].WireID().Wire;
      // find the mean and rms of the charge
      unsigned short npt = 0;
      unsigned short imlast = 0;
      float ave = 0.;
      float rms = 0.;
      // this loop intentionally ignores the Begin hit
      for (unsigned int ii = fcl2hits.size() - 1; ii > 0; --ii) {
        ++npt;
        float chg = fHits[fcl2hits[ii]].Integral();
        ave += chg;
        rms += chg * chg;
        if (npt == fitLen) {
          imlast = ii;
          break;
        }
      }
      float fnpt = npt;
      ave /= fnpt;
      rms = std::sqrt((rms - fnpt * ave * ave) / (fnpt - 1));
      float chgcut = ave + rms;
      float chg;
      unsigned int wire;
      for (unsigned short ii = fcl2hits.size() - 1; ii > imlast; --ii) {
        wire = fHits[fcl2hits[ii]].WireID().Wire;
        chg = fHits[fcl2hits[ii]].Integral();
        if (chg > chgcut) continue;
        xwir.push_back((float)(wire - wire0));
        ychg.push_back(chg);
        ychgerr2.push_back(chg);
      }
      if (ychg.size() < 3) return;
      float intcpt;
      float slope;
      float intcpterr;
      float slopeerr;
      float chidof;
      fLinFitAlg.LinFit(xwir, ychg, ychgerr2, intcpt, slope, intcpterr, slopeerr, chidof);
      if (prt)
        mf::LogVerbatim("CC") << "FitClusterChg wire " << wire0 << " chidof " << (int)chidof
                              << " npt " << xwir.size() << " charge = " << (int)intcpt
                              << " slope = " << (int)slope << " first ave " << (int)ave << " rms "
                              << (int)rms;
      if (chidof > 100.) return;
      // fit must have gone wrong if the fStepCrawlChgRatCut average is greater than
      // the average using all points
      if (intcpt > ave) return;
      // ensure that change does not exceed 30%
      if (fAveChg > 0) {
        ave = intcpt / fAveChg;
        if (ave > 1.3) return;
        if (ave < 0.77) return;
      }
      fAveChg = intcpt;
      fChgSlp = slope;
    }
  } // fitchg

  /////////////////////////////////////////
  void ClusterCrawlerAlg::AddLAHit(unsigned int kwire, bool& ChkCharge, bool& HitOK, bool& SigOK)
  {
    // A variant of AddHit for large angle clusters

    SigOK = false;
    HitOK = false;

    // not in the range of wires with hits
    if (kwire < fFirstWire || kwire > fLastWire) return;

    if (fcl2hits.size() == 0) return;

    // skip bad wire and assume the track was there
    if (WireHitRange[kwire].first == -1) {
      SigOK = true;
      return;
    }
    // return SigOK false if no hit on a good wire
    if (WireHitRange[kwire].first == -2) return;

    unsigned int firsthit = WireHitRange[kwire].first;
    unsigned int lasthit = WireHitRange[kwire].second;

    // max allowable time difference between projected cluster and a hit
    float timeDiff = 40 * AngleFactor(clpar[1]);
    float dtime;

    // the last hit added to the cluster
    unsigned int lastClHit = UINT_MAX;
    if (fcl2hits.size() > 0) {
      lastClHit = fcl2hits[fcl2hits.size() - 1];
      if (lastClHit == 0) {
        fAveChg = fHits[lastClHit].Integral();
        fAveHitWidth = fHits[lastClHit].EndTick() - fHits[lastClHit].StartTick();
      }
    } // fcl2hits.size() > 0

    // the projected time of the cluster on this wire
    float prtime = clpar[0] + ((float)kwire - clpar[2]) * clpar[1];
    float chgWinLo = prtime - fChgNearWindow;
    float chgWinHi = prtime + fChgNearWindow;
    float chgrat, hitWidth;
    float hitWidthCut = 0.5 * fAveHitWidth;
    float cnear = 0;
    //    float fom;
    if (prt)
      mf::LogVerbatim("CC") << "AddLAHit: wire " << kwire << " prtime " << prtime
                            << " max timeDiff " << timeDiff << " fAveChg " << fAveChg;
    unsigned int imbest = 0;
    unsigned int khit;
    for (khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if (inClus[khit] < 0) continue;
      dtime = std::abs(fHits[khit].PeakTime() - prtime);
      hitWidth = fHits[khit].EndTick() - fHits[khit].StartTick();
      chgrat = 1;
      if (ChkCharge && fAveChg > 0) { chgrat = fHits[khit].Integral() / fAveChg; }
      if (prt)
        mf::LogVerbatim("CC") << " Chk W:T " << kwire << ":" << (short)fHits[khit].PeakTime()
                              << " dT " << std::fixed << std::setprecision(1) << dtime << " InClus "
                              << inClus[khit] << " mult " << fHits[khit].Multiplicity() << " width "
                              << (int)hitWidth << " MergeAvail " << mergeAvailable[khit] << " Chi2 "
                              << std::fixed << std::setprecision(2) << fHits[khit].GoodnessOfFit()
                              << " Charge " << (int)fHits[khit].Integral() << " chgrat "
                              << std::fixed << std::setprecision(1) << chgrat << " index " << khit;
      // count charge in the window
      if (fHits[khit].PeakTime() > chgWinLo && fHits[khit].PeakTime() < chgWinHi)
        cnear += fHits[khit].Integral();
      // projected time outside the Signal time window?
      if (prtime < fHits[khit].StartTick() - timeDiff) continue;
      if (prtime > fHits[khit].EndTick() + timeDiff) continue;
      SigOK = true;
      // hit used?
      if (inClus[khit] > 0) continue;
      // ignore narrow width hits
      if (hitWidth < hitWidthCut) continue;
      // ignore very low charge hits
      if (chgrat < 0.1) continue;
      if (dtime < timeDiff) {
        HitOK = true;
        imbest = khit;
        timeDiff = dtime;
      }
    } // khit

    if (prt && !HitOK) mf::LogVerbatim("CC") << " no hit found ";

    if (!HitOK) return;

    if (prt)
      mf::LogVerbatim("CC") << " Pick hit time " << (int)fHits[imbest].PeakTime() << " hit index "
                            << imbest;

    // merge hits in a multiplet?
    short hnear = 0;
    if (lastClHit != UINT_MAX && fHits[imbest].Multiplicity() > 1) {
      bool doMerge = true;
      // Standard code
      // don't merge if we are close to a vertex
      for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if (vtx[ivx].CTP != clCTP) continue;
        if (prt)
          mf::LogVerbatim("CC") << " close vtx chk W:T " << vtx[ivx].Wire << ":"
                                << (int)vtx[ivx].Time;
        if (std::abs(kwire - vtx[ivx].Wire) < 5 &&
            std::abs(int(fHits[imbest].PeakTime() - vtx[ivx].Time)) < 20) {
          if (prt) mf::LogVerbatim("CC") << " Close to a vertex. Don't merge hits";
          doMerge = false;
        }
      } // ivx
      // Decide which hits in the multiplet to merge. Hits that are well
      // separated from each other should not be merged
      if (doMerge) {
        unsigned short nused = 0;
        // the total charge of the hit multiplet
        float multipletChg = 0.;
        float chicut = AngleFactor(clpar[1]) * fHitMergeChiCut * fHits[lastClHit].RMS();
        // look for a big separation between adjacent hits
        std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(imbest);
        for (size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
          if (inClus[jht] < 0) continue;
          if (inClus[jht] == 0)
            multipletChg += fHits[jht].Integral();
          else
            ++nused;
          // check the neighbor hit separation
          if (jht > MultipletRange.first) {
            // pick the larger RMS of the two hits
            float hitRMS = fHits[jht].RMS();
            if (fHits[jht - 1].RMS() > hitRMS) hitRMS = fHits[jht - 1].RMS();
            const float tdiff =
              std::abs(fHits[jht].PeakTime() - fHits[jht - 1].PeakTime()) / hitRMS;
            if (prt) mf::LogVerbatim("CC") << " Hit RMS chisq " << tdiff << " chicut " << chicut;
            if (tdiff > chicut) doMerge = false;
          } // jht > 0
        }   // jht
        if (prt) {
          if (!doMerge) mf::LogVerbatim("CC") << " Hits are well separated. Don't merge them ";
        }
        if (doMerge && nused == 0) {
          // compare the charge with the last hit added?
          if (ChkCharge) {
            // there is a nearby hit
            hnear = 1;
            float chgrat = multipletChg / fHits[lastClHit].Integral();
            if (prt)
              mf::LogVerbatim("CC") << " merge hits charge check " << (int)multipletChg
                                    << " Previous hit charge " << (int)fHits[lastClHit].Integral();
            if (chgrat > 2) doMerge = false;
          }
        } // doMerge && nused == 0
      }   // doMerge true
      if (doMerge) {
        // there is a nearby hit and it will be merged
        hnear = -1;
        bool didMerge;
        MergeHits(imbest, didMerge);
      } // doMerge
    }   // Hits[imbest].Multiplicity() > 1

    // attach to the cluster and fit
    fcl2hits.push_back(imbest);
    FitCluster();
    FitClusterChg();
    chifits.push_back(clChisq);
    hitNear.push_back(hnear);
    // remove the charge of the just added hit
    cnear -= fHits[imbest].Integral();
    if (cnear < 0) cnear = 0;
    // divide by the just added hit charge
    cnear /= fHits[imbest].Integral();
    chgNear.push_back(cnear);
    if (prt) {
      hitWidth = fHits[imbest].EndTick() - fHits[imbest].StartTick();
      mf::LogVerbatim("CC") << " >>LADD" << pass << " W:T " << PrintHit(imbest) << " dT "
                            << timeDiff << " clChisq " << clChisq << " Chg "
                            << (int)fHits[imbest].Integral() << " AveChg " << (int)fAveChg
                            << " width " << (int)hitWidth << " AveWidth " << (int)fAveHitWidth
                            << " fcl2hits size " << fcl2hits.size();
    } // prt
    // decide what to do with a bad fit
    if (clChisq > fChiCut[pass]) {
      FclTrimUS(1);
      FitCluster();
      HitOK = false;
      SigOK = false;
      if (prt) mf::LogVerbatim("CC") << "   LADD- Removed bad hit. Stopped tracking";
    }
  } // AddLAHit()

  /////////////////////////////////////////
  bool ClusterCrawlerAlg::ClusterHitsOK(short nHitChk)
  {
    // Check StartTick and EndTick of hits on adjacent wires overlap as illustrated below.
    // >>>>>> This is OK
    // Wire StartTick    EndTick
    // n    |--------------|
    // n+1              |--------------|
    // n+2                        |--------------|
    // >>>>>> This is NOT OK
    // n    |------|
    // n+1              |-----|
    // n+2                        |------|

    if (fcl2hits.size() == 0) return true;

    unsigned short nHitToChk = fcl2hits.size();
    if (nHitChk > 0) nHitToChk = nHitChk + 1;
    unsigned short ii, indx;

    // require that they overlap
    // add a tolerance to the StartTick - EndTick overlap
    raw::TDCtick_t tol = 30;
    // expand the tolerance for induction planes
    if (plane < geom->TPC(geo::TPCID(cstat, tpc)).Nplanes() - 1) tol = 40;

    bool posSlope =
      (fHits[fcl2hits[0]].PeakTime() > fHits[fcl2hits[fcl2hits.size() - 1]].PeakTime());

    if (prt) {
      for (ii = 0; ii < nHitToChk; ++ii) {
        indx = fcl2hits.size() - 1 - ii;
        mf::LogVerbatim("CC") << "ClusterHitsOK chk " << fHits[fcl2hits[indx]].WireID().Wire
                              << " start " << fHits[fcl2hits[indx]].StartTick() << " peak "
                              << fHits[fcl2hits[indx]].PeakTime() << " end "
                              << fHits[fcl2hits[indx]].EndTick() << " posSlope " << posSlope;
      }
    }

    raw::TDCtick_t hiStartTick, loEndTick;
    for (ii = 0; ii < nHitToChk - 1; ++ii) {
      indx = fcl2hits.size() - 1 - ii;
      // ignore if not on adjacent wires
      if (lar::util::absDiff(fHits[fcl2hits[indx]].WireID().Wire,
                             fHits[fcl2hits[indx - 1]].WireID().Wire) > 1)
        continue;
      hiStartTick =
        std::max(fHits[fcl2hits[indx]].StartTick(), fHits[fcl2hits[indx - 1]].StartTick());
      loEndTick = std::min(fHits[fcl2hits[indx]].EndTick(), fHits[fcl2hits[indx - 1]].EndTick());
      if (posSlope) {
        if (loEndTick + tol < hiStartTick) { return false; }
      }
      else {
        if (loEndTick + tol < hiStartTick) { return false; }
      }
    } // ii
    return true;
  } // ClusterHitsOK

  /////////////////////////////////////////
  void ClusterCrawlerAlg::AddHit(unsigned int kwire, bool& HitOK, bool& SigOK)
  {
    // Add a hit to the cluster if it meets several criteria:
    // similar pulse height to the cluster (if fAveChg is defined)
    // closest hit to the project cluster position.
    // Return SigOK if there is a nearby hit that was missed due to the cuts

    SigOK = false;
    HitOK = false;

    // not in the range of wires with hits
    if (kwire < fFirstWire || kwire > fLastWire) return;

    unsigned int lastClHit = UINT_MAX;
    if (fcl2hits.size() > 0) lastClHit = fcl2hits[fcl2hits.size() - 1];

    // the last hit added to the cluster
    unsigned int wire0 = clpar[2];

    // return if no signal and no hit
    if (fAllowNoHitWire == 0) {
      if (WireHitRange[kwire].first == -2) return;
    }
    else {
      // allow a number of wires with no hits
      if (WireHitRange[kwire].first == -2 && (wire0 - kwire) > fAllowNoHitWire) {
        SigOK = true;
        return;
      }
    }
    // skip bad wire, but assume the track was there
    if (WireHitRange[kwire].first == -1) {
      SigOK = true;
      return;
    }

    unsigned int firsthit = WireHitRange[kwire].first;
    unsigned int lasthit = WireHitRange[kwire].second;

    // the projected time of the cluster on this wire
    float dw = (float)kwire - (float)wire0;
    float prtime = clpar[0] + dw * clpar[1];
    if (prtime < 0 || (unsigned int)prtime > fMaxTime) return;
    // Find the projected time error including the projection error and the
    // error from the last hit added
    float prtimerr2 = std::abs(dw) * clparerr[1] * clparerr[1];

    // apply an angle dependent scale factor to the hit error. Default is very large error
    float hiterr = 10;
    if (lastClHit != UINT_MAX) hiterr = 3 * fHitErrFac * fHits[lastClHit].RMS();
    float err = std::sqrt(prtimerr2 + hiterr * hiterr);
    // Time window for accepting a hit.
    float hitWin = fClProjErrFac * err;

    float prtimeLo = prtime - hitWin;
    float prtimeHi = prtime + hitWin;
    float chgWinLo = prtime - fChgNearWindow;
    float chgWinHi = prtime + fChgNearWindow;
    if (prt) {
      mf::LogVerbatim("CC") << "AddHit: wire " << kwire << " prtime Lo " << (int)prtimeLo
                            << " prtime " << (int)prtime << " Hi " << (int)prtimeHi << " prtimerr "
                            << sqrt(prtimerr2) << " hiterr " << hiterr << " fAveChg "
                            << (int)fAveChg << " fAveHitWidth " << std::setprecision(3)
                            << fAveHitWidth;
    }
    // loop through the hits
    unsigned int imbest = INT_MAX;
    float best = 9999., dtime;
    float cnear = 0;
    float hitTime, hitChg, hitStartTick, hitEndTick;
    for (unsigned int khit = firsthit; khit < lasthit; ++khit) {
      // obsolete hit?
      if (inClus[khit] < 0) continue;
      hitTime = fHits[khit].PeakTime();
      dtime = std::abs(hitTime - prtime);
      if (dtime > 1000) continue;
      hitStartTick = fHits[khit].StartTick();
      hitEndTick = fHits[khit].EndTick();
      // weight by the charge difference
      if (fAveChg > 0) dtime *= std::abs(fHits[khit].Integral() - fAveChg) / fAveChg;
      if (prt && std::abs(dtime) < 100) {
        mf::LogVerbatim("CC") << " Chk W:T " << PrintHit(khit) << " dT " << std::fixed
                              << std::setprecision(1) << (hitTime - prtime) << " InClus "
                              << inClus[khit] << " mult " << fHits[khit].Multiplicity() << " RMS "
                              << std::fixed << std::setprecision(1) << fHits[khit].RMS() << " Chi2 "
                              << std::fixed << std::setprecision(1) << fHits[khit].GoodnessOfFit()
                              << " Charge " << (int)fHits[khit].Integral() << " Peak " << std::fixed
                              << std::setprecision(1) << fHits[khit].PeakAmplitude() << " LoT "
                              << (int)hitStartTick << " HiT " << (int)hitEndTick << " index "
                              << khit;
      }
      // count charge in the window
      if (fHits[khit].StartTick() > chgWinLo && fHits[khit].EndTick() < chgWinHi)
        cnear += fHits[khit].Integral();
      // check for signal
      if (prtimeHi < hitStartTick) continue;
      if (prtimeLo > hitEndTick) continue;
      SigOK = true;
      // check for good hit
      if (hitTime < prtimeLo) continue;
      if (hitTime > prtimeHi) continue;
      // hit used?
      if (inClus[khit] > 0) continue;
      if (dtime < best) {
        best = dtime;
        imbest = khit;
      }
    } // khit

    if (!SigOK) {
      if (fAllowNoHitWire == 0) return;
      if (prt)
        mf::LogVerbatim("CC") << " wire0 " << wire0 << " kwire " << kwire << " max "
                              << fAllowNoHitWire << " imbest " << imbest;
      if ((wire0 - kwire) > fAllowNoHitWire) return;
      SigOK = true;
    }

    if (imbest == INT_MAX) return;

    recob::Hit const& hit = fHits[imbest];
    hitChg = hit.Integral();

    if (prt) mf::LogVerbatim("CC") << " Best hit time " << (int)hit.PeakTime();

    short hnear = 0;
    // merge hits in a doublet?
    bool didMerge = false;
    if (lastClHit != UINT_MAX && fAveHitWidth > 0 && fHitMergeChiCut > 0 &&
        hit.Multiplicity() == 2) {
      bool doMerge = true;
      for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if (std::abs(kwire - vtx[ivx].Wire) < 10 &&
            std::abs(int(hit.PeakTime() - vtx[ivx].Time)) < 20) {
          doMerge = false;
          break;
        }
      } // ivx
      // quit if localindex does not make sense.
      if (hit.LocalIndex() != 0 && imbest == 0) doMerge = false;
      if (doMerge) {
        // find the neighbor hit
        unsigned int oht;
        if (hit.LocalIndex() == 0) { oht = imbest + 1; }
        else {
          oht = imbest - 1;
        } // hit.LocalIndex() == 0
        // check the hit time separation
        recob::Hit const& other_hit = fHits[oht];
        float hitSep = std::abs(hit.PeakTime() - other_hit.PeakTime());
        hitSep /= hit.RMS();
        // check the the charge similarity
        float totChg = hitChg + other_hit.Integral();
        float lastHitChg = fAveChg;
        if (lastHitChg < 0) lastHitChg = fHits[lastClHit].Integral();
        hnear = 1;
        if (prt)
          mf::LogVerbatim("CC") << " Chk hit merge hitsep " << hitSep << " dChg "
                                << std::abs(totChg - lastHitChg) << " Cut "
                                << std::abs(hit.Integral() - lastHitChg);
        if (inClus[oht] == 0 && hitSep < fHitMergeChiCut) {
          if (prt) mf::LogVerbatim("CC") << "  Merging hit doublet " << imbest;
          MergeHits(imbest, didMerge);
          if (prt && !didMerge) mf::LogVerbatim("CC") << "  Hit merge failed ";
        } // not in a cluster, hitSep OK, total charge OK
      }   // doMerge
    }     // fHitMergeChiCut > 0 && hit.Multiplicity() == 2

    // Make a charge similarity cut if the average charge is defined
    bool fitChg = true;
    if (fAveChg > 0.) {

      float chgrat = (hitChg - fAveChg) / fAveChg;
      if (prt) mf::LogVerbatim("CC") << " Chgrat " << std::setprecision(2) << chgrat;

      // charge is way too high?
      if (chgrat > 3 * fChgCut[pass]) {
        if (prt)
          mf::LogVerbatim("CC") << " fails 3 x high charge cut " << fChgCut[pass] << " on pass "
                                << pass;
        return;
      }

      // Determine if the last hit added was a large (low) charge hit
      // This will be used to prevent adding large (low) charge hits on two
      // consecutive fits. This cut is only applied to hits on adjacent wires
      float bigchgcut = 1.5 * fChgCut[pass];
      bool lasthitbig = false;
      bool lasthitlow = false;
      if (lastClHit != UINT_MAX && lar::util::absDiff(wire0, kwire) == 1) {
        float lastchgrat = (fHits[lastClHit].Integral() - fAveChg) / fAveChg;
        lasthitbig = (lastchgrat > bigchgcut);
        lasthitlow = (lastchgrat < -fChgCut[pass]);
      }

      // the last hit added was low charge and this one is as well
      if (lasthitlow && chgrat < -fChgCut[pass]) {
        if (prt) mf::LogVerbatim("CC") << " fails low charge cut. Stop crawling.";
        SigOK = false;
        return;
      } // lasthitlow

      // the last hit was high charge and this one is also
      if (lasthitbig && chgrat > fChgCut[pass]) {
        if (prt)
          mf::LogVerbatim("CC") << " fails 2nd hit high charge cut. Last hit was high also. ";
        return;
      } // lasthitbig

      // require that large charge hits have a very good projection error
      if (chgrat > fChgCut[pass]) {
        if (best > 2 * err) {
          if (prt) mf::LogVerbatim("CC") << " high charge && bad dT= " << best << " err= " << err;
          return;
        }
      } // chgrat > fChgCut[pass]

      // decide whether to fit the charge
      fitChg = (chgrat < std::abs(fChgCut[pass]));
    } // fAveChg > 0

    // we now have a hit that meets all the criteria. Fit it
    fcl2hits.push_back(imbest);
    // This is strictly only necessary when calling AddHit for seed clusters
    if (fcl2hits.size() == 3) std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
    FitCluster();
    chifits.push_back(clChisq);
    hitNear.push_back(hnear);
    // remove the charge of the just added hit
    cnear -= fHits[imbest].Integral();
    if (cnear < 0) cnear = 0;
    // divide by the just added hit charge
    cnear /= fHits[imbest].Integral();
    chgNear.push_back(cnear);
    // nearby hit check
    //    ChkClusterNearbyHits(prt);
    HitOK = true;

    if (chgNear.size() != fcl2hits.size()) {
      mf::LogError("CC") << "AddHit: Bad length";
      return;
    }

    if (prt)
      mf::LogVerbatim("CC") << " >>ADD" << pass << " W:T " << PrintHit(imbest) << " dT " << best
                            << " clChisq " << clChisq << " Chg " << (int)fHits[imbest].Integral()
                            << " AveChg " << (int)fAveChg << " fcl2hits size " << fcl2hits.size();

    if (!fitChg) return;
    if (prt) mf::LogVerbatim("CC") << " Fit charge ";
    FitClusterChg();
  } // AddHit()

  //////////////////////////////////////
  void ClusterCrawlerAlg::ChkClusterNearbyHits(bool prt)
  {
    // analyze the hitnear vector
    //  0 = no nearby hit exists
    //  1 = a nearby hit exists but was not merged
    // -1 = a nearby hit was merged

    if (fHitMergeChiCut <= 0) return;

    if (hitNear.size() != fcl2hits.size()) {
      mf::LogWarning("CC") << "Coding error: hitNear size != fcl2hits";
      return;
    }

    // Analyze the last 6 hits added but don't consider the first few hits
    if (hitNear.size() < 12) return;

    // TODO move into loops
    unsigned short ii, indx;
    unsigned short merged = 0;
    unsigned short notmerged = 0;
    for (ii = 0; ii < 6; ++ii) {
      indx = hitNear.size() - 1 - ii;
      if (hitNear[indx] > 0) ++notmerged;
      if (hitNear[indx] < 0) ++merged;
    }

    if (prt)
      mf::LogVerbatim("CC") << "ChkClusterNearbyHits: nearby hits merged " << merged
                            << " not merged " << notmerged;

    if (notmerged < 2) return;

    // a number of nearby hits were not merged while crawling, so the
    // average charge is probably wrong. Look at the last 6 hits added
    // and merge them if they are close
    bool didMerge;
    for (ii = 0; ii < 6; ++ii) {
      indx = fcl2hits.size() - 1 - ii;
      const unsigned int iht = fcl2hits[indx];
      recob::Hit const& hit = fHits[iht];
      if (hit.Multiplicity() == 2) {
        // quit if localindex does not make sense.
        if (hit.LocalIndex() != 0 && iht == 0) continue;
        // hit doublet. Get the index of the other hit
        unsigned int oht;
        if (hit.LocalIndex() == 0) { oht = iht + 1; }
        else {
          oht = iht - 1;
        } // hit.LocalIndex() == 0
        recob::Hit const& other_hit = fHits[oht];
        // TODO use Hit::TimeDistanceAsRMS()
        float hitSep = std::abs(hit.PeakTime() - other_hit.PeakTime());
        hitSep /= hit.RMS();
        if (hitSep < fHitMergeChiCut && inClus[oht] == 0) {
          if (prt)
            mf::LogVerbatim("CC") << "Merging hit doublet " << iht << " W:T "
                                  << fHits[iht].WireID().Wire << ":" << fHits[iht].PeakTime();
          MergeHits(iht, didMerge);
          if (didMerge) hitNear[indx] = -1;
        } // hitSep OK and not in a cluster
      }   // hit doublet
    }     // ii

    // now re-fit
    FitCluster();
    FitClusterChg();

    if (prt) mf::LogVerbatim("CC") << "ChkClusterNearbyHits refit cluster. fAveChg= " << fAveChg;

  } // ChkClusterHitNear()

  //////////////////////////////////////
  void ClusterCrawlerAlg::FitVtx(unsigned short iv)
  {
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> ey2;
    float arg;

    // don't fit fixed vertices
    if (vtx[iv].Fixed) return;

    // Set this large in case something bad happens
    vtx[iv].ChiDOF = 99;

    // make a list of clusters
    unsigned short icl;
    std::vector<unsigned short> vcl;
    for (icl = 0; icl < tcl.size(); ++icl) {
      if (tcl[icl].ID < 0) continue;
      if (tcl[icl].CTP != vtx[iv].CTP) continue;
      if (tcl[icl].EndVtx == iv) vcl.push_back(icl);
      if (tcl[icl].BeginVtx == iv) vcl.push_back(icl);
    }

    vtx[iv].NClusters = vcl.size();

    if (vcl.size() == 0) return;

    // don't let the time error be less than the expected
    // time error of hits on a cluster. This is crude by
    // probably good enough
    icl = vcl[0];
    unsigned short indx = tcl[icl].tclhits.size() - 1;
    unsigned int hit = tcl[icl].tclhits[indx];
    float minTimeErr = fHitErrFac * fHits[hit].RMS() * fHits[hit].Multiplicity();

    if (vcl.size() == 1) {
      icl = vcl[0];
      // Put the vertex at the appropriate end of the cluster
      if (tcl[icl].EndVtx == iv) {
        vtx[iv].Wire = tcl[icl].EndWir;
        vtx[iv].WireErr = 1;
        vtx[iv].Time = tcl[icl].EndTim;
        // set the vertex time error to the hit error used for fitting
        indx = tcl[icl].tclhits.size() - 1;
        hit = tcl[icl].tclhits[indx];
        vtx[iv].TimeErr = fHitErrFac * fHits[hit].RMS() * fHits[hit].Multiplicity();
        vtx[iv].ChiDOF = 0;
      }
      if (tcl[icl].BeginVtx == iv) {
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

    std::vector<double> slps;
    std::vector<double> slperrs;
    for (unsigned short ii = 0; ii < vcl.size(); ++ii) {
      icl = vcl[ii];
      if (tcl[icl].EndVtx == iv) {
        x.push_back(tcl[icl].EndSlp);
        slps.push_back(tcl[icl].EndSlp);
        slperrs.push_back(tcl[icl].EndSlpErr);
        arg = tcl[icl].EndSlp * tcl[icl].EndWir - tcl[icl].EndTim;
        y.push_back(arg);
        if (tcl[icl].EndSlpErr > 0.) { arg = tcl[icl].EndSlpErr * tcl[icl].EndWir; }
        else {
          arg = .1 * tcl[icl].EndWir;
        }
        ey2.push_back(arg * arg);
      }
      else if (tcl[icl].BeginVtx == iv) {
        x.push_back(tcl[icl].BeginSlp);
        slps.push_back(tcl[icl].BeginSlp);
        slperrs.push_back(tcl[icl].BeginSlpErr);
        arg = tcl[icl].BeginSlp * tcl[icl].BeginWir - tcl[icl].BeginTim;
        y.push_back(arg);
        if (tcl[icl].BeginSlpErr > 0.) { arg = tcl[icl].BeginSlpErr * tcl[icl].BeginWir; }
        else {
          arg = .1 * tcl[icl].BeginWir;
        }
        ey2.push_back(arg * arg);
      }
    } // ii
    if (x.size() < 2) return;

    // calculate error
    double sumerr = 0, cnt = 0;
    for (unsigned short ii = 0; ii < slps.size() - 1; ++ii) {
      for (unsigned short jj = ii + 1; jj < slps.size(); ++jj) {
        arg = std::min(slperrs[ii], slperrs[jj]);
        arg /= (slps[ii] - slps[jj]);
        sumerr += arg * arg;
        ++cnt;
      } // jj
    }   // ii
    sumerr /= cnt;

    float vTime = 0.;
    float vTimeErr = 0.;
    float vWire = 0.;
    float vWireErr = 0.;
    float chiDOF;
    fLinFitAlg.LinFit(x, y, ey2, vTime, vWire, vTimeErr, vWireErr, chiDOF);
    if (chiDOF > 900) return;
    vTime = -vTime;
    // a crazy time from the fit?
    if (vTime < 0 || vTime > fMaxTime) return;
    // a crazy wire from the fit?
    geo::PlaneID iplID = DecodeCTP(vtx[iv].CTP);
    if (vWire < 0 || vWire > geom->Nwires(iplID)) return;
    vtx[iv].ChiDOF = chiDOF;
    vtx[iv].Wire = vWire;
    vtx[iv].Time = vTime;
    vtx[iv].WireErr = vWire * sqrt(sumerr);
    vtx[iv].TimeErr = vTime * fabs(sumerr);
    // set minimum wire error
    if (vtx[iv].WireErr < 1) vtx[iv].WireErr = 2;
    // set minimum time error
    if (vtx[iv].TimeErr < minTimeErr) vtx[iv].TimeErr = minTimeErr;

  } // FitVtx

  //////////////////////////////////////
  void ClusterCrawlerAlg::Vtx3ClusterMatch(detinfo::DetectorPropertiesData const& det_prop,
                                           geo::TPCID const& tpcid)
  {
    // Look for clusters that end/begin near the expected wire/time
    // for incomplete 3D vertices
    if (empty(vtx3)) return;

    const unsigned int cstat = tpcid.Cryostat;
    const unsigned int tpc = tpcid.TPC;

    unsigned int thePlane, theWire;
    float theTime;
    int dwb, dwe;

    for (unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
      // A complete 3D vertex with matching 2D vertices in all planes?
      if (vtx3[ivx].Wire < 0) continue;
      if (vtx3[ivx].CStat != cstat || vtx3[ivx].TPC != tpc) continue;
      // Find the plane that is missing a 2D vertex
      thePlane = 3;
      theWire = vtx3[ivx].Wire;
      for (plane = 0; plane < 3; ++plane) {
        if (vtx3[ivx].Ptr2D[plane] >= 0) continue;
        thePlane = plane;
        break;
      } // plane
      if (thePlane > 2) continue;
      theTime = det_prop.ConvertXToTicks(vtx3[ivx].X, thePlane, tpc, cstat);
      clCTP = EncodeCTP(cstat, tpc, thePlane);
      // Create a new 2D vertex and see how many clusters we can attach to it
      VtxStore vnew;
      vnew.Wire = theWire;
      vnew.Time = theTime;
      vnew.CTP = clCTP;
      vnew.Topo = 7;
      vnew.Fixed = false;
      vtx.push_back(vnew);
      unsigned short ivnew = vtx.size() - 1;
      std::vector<short> vclIndex;
      for (unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if (tcl[icl].ID < 0) continue;
        if (tcl[icl].CTP != clCTP) continue;
        dwb = lar::util::absDiff(theWire, tcl[icl].BeginWir);
        dwe = lar::util::absDiff(theWire, tcl[icl].EndWir);
        // rough cut to start
        if (dwb > 10 && dwe > 10) continue;
        if (dwb < dwe && dwb < 10 && tcl[icl].BeginVtx < 0) {
          // cluster begin is closer
          if (theWire < tcl[icl].BeginWir + 5) continue;
          if (ClusterVertexChi(icl, 0, ivnew) > fVertex3DCut) continue;
          tcl[icl].BeginVtx = ivnew;
          vclIndex.push_back(icl);
        }
        else if (dwe < 10 && tcl[icl].EndVtx < 0) {
          // cluster end is closer
          if (theWire > tcl[icl].EndWir - 5) continue;
          if (ClusterVertexChi(icl, 1, ivnew) > fVertex3DCut) continue;
          tcl[icl].EndVtx = ivnew;
          vclIndex.push_back(icl);
        } // dwb/dwe check
      }   // icl
      bool goodVtx = false;
      if (vclIndex.size() > 0) {
        FitVtx(ivnew);
        goodVtx = (vtx[ivnew].ChiDOF < fVertex3DCut);
        vtx3[ivx].Ptr2D[thePlane] = ivnew;
      }
      if (goodVtx) {
        vtx3[ivx].Ptr2D[thePlane] = ivnew;
        vtx3[ivx].Wire = -1;
      }
      else {
        // clobber the vertex
        vtx.pop_back();
        for (unsigned short ii = 0; ii < vclIndex.size(); ++ii) {
          unsigned short icl = vclIndex[ii];
          if (tcl[icl].BeginVtx == ivnew) tcl[icl].BeginVtx = -99;
          if (tcl[icl].EndVtx == ivnew) tcl[icl].EndVtx = -99;
        } // ii
      }
    } // ivx
  }   // Vtx3ClusterMatch

  //////////////////////////////////////
  void ClusterCrawlerAlg::Vtx3ClusterSplit(detinfo::DetectorPropertiesData const& det_prop,
                                           geo::TPCID const& tpcid)
  {
    // Try to split clusters in a view in which there is no 2D vertex
    // assigned to a 3D vertex
    if (empty(vtx3)) return;

    const unsigned int cstat = tpcid.Cryostat;
    const unsigned int tpc = tpcid.TPC;

    vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);

    unsigned int lastplane = 5, kcl, kclID;
    float theTime;
    unsigned int thePlane, theWire, plane;
    unsigned int loWire, hiWire;

    for (unsigned short ivx = 0; ivx < vtx3.size(); ++ivx) {
      if (vtx3[ivx].CStat != cstat || vtx3[ivx].TPC != tpc) continue;
      // Complete 3D vertex with matching 2D vertices in all planes?
      if (vtxprt)
        mf::LogVerbatim("CC") << "Vtx3ClusterSplit ivx " << ivx << " Ptr2D " << vtx3[ivx].Ptr2D[0]
                              << " " << vtx3[ivx].Ptr2D[1] << " " << vtx3[ivx].Ptr2D[2] << " wire "
                              << vtx3[ivx].Wire;
      if (vtx3[ivx].Wire < 0) continue;
      // find the plane that needs to be studied
      thePlane = 3;
      theWire = vtx3[ivx].Wire;
      for (plane = 0; plane < 3; ++plane) {
        if (vtx3[ivx].Ptr2D[plane] >= 0) continue;
        thePlane = plane;
        break;
      } // plane
      if (thePlane > 2) continue;
      theTime = det_prop.ConvertXToTicks(
        (double)vtx3[ivx].X, (int)thePlane, (int)tpcid.TPC, (int)tpcid.Cryostat);
      // get the hit range if necessary
      if (thePlane != lastplane) {
        clCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, thePlane);
        GetHitRange(clCTP);
        lastplane = thePlane;
      }
      // make a list of clusters that have hits near this point on nearby wires
      std::vector<short> clIDs;
      if (theWire > fFirstWire + 5) { loWire = theWire - 5; }
      else {
        loWire = fFirstWire;
      }
      if (theWire < fLastWire - 5) { hiWire = theWire + 5; }
      else {
        hiWire = fLastWire;
      }
      if (vtxprt)
        mf::LogVerbatim("CC") << "3DVtx " << ivx << " look for cluster hits near P:W:T " << thePlane
                              << ":" << theWire << ":" << (int)theTime << " Wire range " << loWire
                              << " to " << hiWire;
      for (unsigned int wire = loWire; wire < hiWire; ++wire) {
        // ignore dead wires or wires with no hits
        if (WireHitRange[wire].first < 0) continue;
        unsigned int firsthit = WireHitRange[wire].first;
        unsigned int lasthit = WireHitRange[wire].second;
        for (unsigned int khit = firsthit; khit < lasthit; ++khit) {
          // ignore obsolete and un-assigned hits
          if (inClus[khit] <= 0) continue;
          if ((unsigned int)inClus[khit] > tcl.size() + 1) {
            mf::LogError("CC") << "Invalid hit InClus. " << khit << " " << inClus[khit];
            continue;
          }
          // check an expanded time range
          if (theTime < fHits[khit].StartTick() - 20) continue;
          if (theTime > fHits[khit].EndTick() + 20) continue;
          kclID = inClus[khit];
          kcl = kclID - 1;
          // ignore obsolete clusters
          if (tcl[kcl].ID < 0) continue;
          // ignore short clusters
          if (tcl[kcl].tclhits.size() < 6) continue;
          // ignore long straight clusters
          if (tcl[kcl].tclhits.size() > 100 && std::abs(tcl[kcl].BeginAng - tcl[kcl].EndAng) < 0.1)
            continue;
          // put the cluster in the list if it's not there already
          if (vtxprt)
            mf::LogVerbatim("CC") << "Bingo " << ivx << " plane " << thePlane << " wire " << wire
                                  << " hit " << fHits[khit].WireID().Wire << ":"
                                  << (int)fHits[khit].PeakTime() << " inClus " << inClus[khit]
                                  << " P:W:T " << thePlane << ":" << tcl[kcl].BeginWir << ":"
                                  << (int)tcl[kcl].BeginTim;
          if (std::find(clIDs.begin(), clIDs.end(), kclID) == clIDs.end()) {
            clIDs.push_back(kclID);
          } // std::find
        }   // khit
      }     // wire
      if (clIDs.size() == 0) continue;
      if (vtxprt)
        for (unsigned int ii = 0; ii < clIDs.size(); ++ii)
          mf::LogVerbatim("CC") << " clIDs " << clIDs[ii];

      unsigned short ii, icl, jj;
      unsigned int iht;
      short nhitfit;
      bool didit;
      // find a reasonable time error using the 2D vertices that comprise this
      // incomplete 3D vertex
      float tErr = 1;
      unsigned short i2Dvx = 0;
      for (ii = 0; ii < 3; ++ii) {
        if (ii == thePlane) continue;
        i2Dvx = vtx3[ivx].Ptr2D[ii];
        if (i2Dvx > vtx.size() - 1) {
          mf::LogError("CC") << "Vtx3ClusterSplit: Coding error";
          return;
        }
        if (vtx[i2Dvx].TimeErr > tErr) tErr = vtx[i2Dvx].TimeErr;
      } // ii -> plane

      // do a local fit near the crossing point and make a tighter cut
      for (ii = 0; ii < clIDs.size(); ++ii) {
        icl = clIDs[ii] - 1;
        didit = false;
        for (jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
          iht = tcl[icl].tclhits[jj];
          if (fHits[iht].WireID().Wire < theWire) {
            nhitfit = 3;
            if (jj > 3) nhitfit = -3;
            FitClusterMid(icl, iht, nhitfit);
            float doca = DoCA(-1, 1, theWire, theTime);
            if (vtxprt)
              mf::LogVerbatim("CC") << " cls " << icl << " DoCA " << doca << " tErr " << tErr;
            if ((doca / tErr) > 2) clIDs[ii] = -1;
            didit = true;
            break;
          } // fHits[iht].WireID().Wire < theWire
          if (didit) break;
        } // jj
        if (didit) break;
      } // ii
      if (vtxprt) {
        mf::LogVerbatim("CC") << "clIDs after fit " << clIDs.size();
        for (ii = 0; ii < clIDs.size(); ++ii)
          mf::LogVerbatim("CC") << " clIDs " << clIDs[ii];
      }

      // see if any candidates remain
      unsigned short nok = 0;
      for (ii = 0; ii < clIDs.size(); ++ii)
        if (clIDs[ii] >= 0) ++nok;
      if (nok == 0) continue;

      // make a new 2D vertex
      VtxStore vnew;
      vnew.Wire = theWire;
      vnew.WireErr = 1;
      vnew.Time = theTime;
      vnew.TimeErr = 1;
      vnew.Topo = 8;
      vnew.CTP = clCTP;
      vnew.Fixed = false;
      vtx.push_back(vnew);
      // update the 2D -> 3D vertex pointer
      unsigned short ivnew = vtx.size() - 1;
      if (vtxprt)
        mf::LogVerbatim("CC") << "Make new 2D vtx " << ivnew << " in plane " << thePlane
                              << " from 3D vtx " << ivx;
      vtx3[ivx].Ptr2D[thePlane] = ivnew;
      // either split or attach clusters to this vertex
      for (ii = 0; ii < clIDs.size(); ++ii) {
        if (clIDs[ii] < 0) continue;
        icl = clIDs[ii] - 1;
        // find the split position
        // split the cluster. Find the split position
        unsigned short pos = 0;
        for (unsigned short jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
          iht = tcl[icl].tclhits[jj];
          if (fHits[iht].WireID().Wire < theWire) {
            pos = jj;
            break;
          }
        } // jj
        if (pos == 0) {
          // vertex is DS of the cluster Begin
          tcl[icl].BeginVtx = ivnew;
          if (vtxprt) mf::LogVerbatim("CC") << "Attach to Begin " << icl;
        }
        else if (pos > tcl[icl].tclhits.size()) {
          // vertex is US of the cluster Eend
          tcl[icl].EndVtx = ivnew;
          if (vtxprt) mf::LogVerbatim("CC") << "Attach to End " << icl;
        }
        else {
          // vertex is in the middle of the cluster
          if (vtxprt) mf::LogVerbatim("CC") << "Split cluster " << clIDs[ii] << " at pos " << pos;
          if (!SplitCluster(icl, pos, ivnew)) {
            if (vtxprt) mf::LogVerbatim("CC") << "SplitCluster failed";
            continue;
          }
          tcl[icl].ProcCode += 10000;
          tcl[tcl.size() - 1].ProcCode += 10000;
        } // pos check
      }   // ii
      // Fit the vertex position
      FitVtx(ivnew);
      if (vtx[ivnew].ChiDOF < 5 && vtx[ivnew].WireErr < fVertex2DWireErrCut) {
        // mark the 3D vertex as complete
        vtx3[ivx].Wire = -1;
      }
      else {
        if (vtxprt)
          mf::LogVerbatim("CC") << "Bad vtx fit " << ivnew << " ChiDOF " << vtx[ivnew].ChiDOF
                                << " WireErr " << vtx[ivnew].WireErr;
        // Recover (partially) from a bad fit. Leave the ProcCode as-is to trace this problem
        vtx.pop_back();
        vtx3[ivx].Ptr2D[thePlane] = -1;
        // find the cluster - vertex associations
        for (jj = 0; jj < tcl.size(); ++jj) {
          if (tcl[jj].BeginVtx == ivnew) tcl[jj].BeginVtx = -99;
          if (tcl[jj].EndVtx == ivnew) tcl[jj].EndVtx = -99;
        } // jj
      }
    } // ivx

  } // Vtx3ClusterSplit()

  //////////////////////////////////////
  void ClusterCrawlerAlg::FindHammerClusters(detinfo::DetectorPropertiesData const& det_prop)
  {
    // look for a long cluster that stops at a short cluster in two views. This can occur in a CCmu
    // interaction where two protons are emitted back-to-back and are therefore reconstructed as one cluster
    // This routine looks for this signature, and if found, splits the short clusters and creates a new 3D vertex.
    // This routine only considers the case where the long cluster intersects the short cluster at the US (End) end.

    unsigned int nPln = geom->TPC(geo::TPCID(cstat, tpc)).Nplanes();
    if (nPln != 3) return;

    float ew1, ew2, bw2, fvw;

    struct Hammer {
      bool Used;
      unsigned int Wire; // intersection point of the long cluster and the short cluster
      float Tick;        // intersection point of the long cluster and the short cluster
      float X;
      unsigned short longClIndex;
      unsigned short shortClIndex;
      unsigned short splitPos;
    };
    std::array<std::vector<Hammer>, 3> hamrVec;

    unsigned int ipl;
    bool useit = false;
    for (ipl = 0; ipl < 3; ++ipl) {
      clCTP = EncodeCTP(cstat, tpc, ipl);
      for (unsigned short ic1 = 0; ic1 < tcl.size(); ++ic1) {
        if (tcl[ic1].ID < 0) continue;
        // require a long cluster
        if (tcl[ic1].tclhits.size() < 20) continue;
        if (tcl[ic1].CTP != clCTP) continue;
        // ignore long clusters with an End vertex assignment
        if (tcl[ic1].EndVtx >= 0) continue;
        ew1 = tcl[ic1].EndWir;
        for (unsigned short ic2 = 0; ic2 < tcl.size(); ++ic2) {
          if (tcl[ic2].ID < 0) continue;
          // require a short cluster
          if (tcl[ic2].tclhits.size() > 20) continue;
          // but not too short cluster
          if (tcl[ic2].tclhits.size() < 6) continue;
          if (tcl[ic2].CTP != clCTP) continue;
          ew2 = tcl[ic2].EndWir;
          bw2 = tcl[ic2].BeginWir;
          // check the US end. The End of the long cluster must lie between the Begin and End wire of the
          // short cluster
          if (ew1 < ew2 || ew1 > bw2) continue;
          // look for intersection of the two clusters
          float best = 10;
          short ibst = -1;
          unsigned short spos = 0;
          for (unsigned short ii = 0; ii < tcl[ic2].tclhits.size(); ++ii) {
            unsigned int iht = tcl[ic2].tclhits[ii];
            float dw = fHits[iht].WireID().Wire - tcl[ic1].EndWir;
            float dt = fabs(fHits[iht].PeakTime() - tcl[ic1].EndTim - tcl[ic1].EndSlp * dw);
            if (dt < best) {
              best = dt;
              ibst = iht;
              spos = ii;
            }
          } // ii
          if (ibst < 0) continue;
          fvw = 0.5 + fHits[ibst].WireID().Wire;
          Hammer aHam;
          aHam.Used = false;
          aHam.Wire = (0.5 + fvw);
          aHam.Tick = fHits[ibst].PeakTime();
          aHam.X = det_prop.ConvertTicksToX((double)aHam.Tick, (int)ipl, (int)tpc, (int)cstat);
          aHam.longClIndex = ic1;
          aHam.shortClIndex = ic2;
          aHam.splitPos = spos;
          unsigned short indx = hamrVec[ipl].size();
          hamrVec[ipl].resize(indx + 1);
          hamrVec[ipl][indx] = aHam;
          useit = true;
        } // ic2
        if (useit) break;
      } // ic1
    }   // ipl

    unsigned short noham = 0;
    for (ipl = 0; ipl < 3; ++ipl)
      if (hamrVec[ipl].size() == 0) ++noham;
    if (noham > 1) return;

    // Y,Z limits of the detector

    geo::TPCID const tpcid(cstat, tpc);
    const geo::TPCGeo& thetpc = geom->TPC(tpcid);
    auto const world = thetpc.GetCenter();
    float YLo = world.Y() - thetpc.HalfHeight() + 1;
    float YHi = world.Y() + thetpc.HalfHeight() - 1;
    float ZLo = world.Z() - thetpc.Length() / 2 + 1;
    float ZHi = world.Z() + thetpc.Length() / 2 - 1;

    // Match in 3D
    float dX;
    double y, z;
    unsigned short icl, jpl, jcl, kpl, splitPos;
    for (ipl = 0; ipl < 3; ++ipl) {
      if (hamrVec[ipl].size() == 0) continue;
      jpl = (ipl + 1) % nPln;
      kpl = (jpl + 1) % nPln;
      for (unsigned short ii = 0; ii < hamrVec[ipl].size(); ++ii) {
        if (hamrVec[ipl][ii].Used) continue;
        for (unsigned short jj = 0; jj < hamrVec[jpl].size(); ++jj) {
          if (hamrVec[jpl][jj].Used) continue;
          dX = hamrVec[ipl][ii].X - hamrVec[jpl][jj].X;
          if (fabs(dX) > fVertex3DCut) continue;
          geo::PlaneID const plane_i{tpcid, ipl};
          geo::PlaneID const plane_j{tpcid, jpl};
          geom->IntersectionPoint(geo::WireID{plane_i, hamrVec[ipl][ii].Wire},
                                  geo::WireID{plane_j, hamrVec[jpl][jj].Wire},
                                  y,
                                  z);
          if (y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
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
          newVtx3.ZErr = 1; // TODO
          newVtx3.CStat = cstat;
          newVtx3.TPC = tpc;

          // make 2D vertex in ipl
          VtxStore newVtx2;
          newVtx2.Wire = hamrVec[ipl][ii].Wire;
          newVtx2.WireErr = 2;
          newVtx2.Time = hamrVec[ipl][ii].Tick;
          newVtx2.TimeErr = 5;
          newVtx2.Topo = 6;
          newVtx2.Fixed = false;
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
          if (!SplitCluster(icl, splitPos, ivnew)) return;

          // make 2D vertex in jpl
          newVtx2.Wire = hamrVec[jpl][jj].Wire;
          newVtx2.Time = hamrVec[jpl][jj].Tick;
          newVtx2.Topo = 6;
          jcl = hamrVec[jpl][jj].longClIndex;
          newVtx2.CTP = tcl[jcl].CTP;
          vtx.push_back(newVtx2);
          ivnew = vtx.size() - 1;
          // associate the new vertex with the long cluster
          tcl[jcl].EndVtx = ivnew;
          // stash the index in newVtx3
          newVtx3.Ptr2D[jpl] = (short)(vtx.size() - 1);
          // split the short cluster and associate the new clusters with the new
          // vtx
          jcl = hamrVec[jpl][jj].shortClIndex;
          splitPos = hamrVec[jpl][jj].splitPos;
          if (!SplitCluster(jcl, splitPos, vtx.size() - 1)) return;
          FitVtx(ivnew);
          // set the kpl 2D vertex index < 0. Let follow-on code find the 3rd
          // plane vertex
          newVtx3.Ptr2D[kpl] = -1;
          geo::Point_t const WPos{0, y, z};
          try {
            newVtx3.Wire = geom->NearestWireID(WPos, geo::PlaneID{cstat, tpc, kpl}).Wire;
          }
          catch (geo::InvalidWireError const& e) {
            newVtx3.Wire = e.suggestedWireID().Wire; // pick the closest valid wire
          }
          vtx3.push_back(newVtx3);
        } // jj
      }   // ii
    }

  } // FindHammerClusters

  //////////////////////////////////////
  void ClusterCrawlerAlg::VtxMatch(detinfo::DetectorPropertiesData const& det_prop,
                                   geo::TPCID const& tpcid)
  {
    // Create 3D vertices from 2D vertices. 3D vertices that are matched
    // in all three planes have Ptr2D >= 0 for all planes

    geo::TPCGeo const& TPC = geom->TPC(tpcid);

    // Y,Z limits of the detector
    auto const world = TPC.GetCenter();

    // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
    float YLo = world.Y() - TPC.HalfHeight() + 1;
    float YHi = world.Y() + TPC.HalfHeight() - 1;
    float ZLo = world.Z() - TPC.Length() / 2 + 1;
    float ZHi = world.Z() + TPC.Length() / 2 - 1;

    vtxprt = (fDebugPlane >= 0) && (fDebugHit == 6666);

    if (vtxprt) {
      mf::LogVerbatim("CC") << "Inside VtxMatch";
      PrintVertices();
    }

    // wire spacing in cm
    float wirePitch = geom->WirePitch(geo::PlaneID{tpcid, 0});

    // fill temp vectors of 2D vertex X and X errors
    std::vector<float> vX(vtx.size());
    std::vector<float> vXsigma(vtx.size());
    float vXp;
    for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if (vtx[ivx].NClusters == 0) continue;
      geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
      if (iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
      // Convert 2D vertex time error to X error
      vX[ivx] =
        det_prop.ConvertTicksToX((double)vtx[ivx].Time, (int)iplID.Plane, (int)tpc, (int)cstat);
      vXp = det_prop.ConvertTicksToX(
        (double)(vtx[ivx].Time + vtx[ivx].TimeErr), (int)iplID.Plane, (int)tpc, (int)cstat);
      vXsigma[ivx] = fabs(vXp - vX[ivx]);
    } // ivx

    // create a array/vector of 2D vertex indices in each plane
    std::array<std::vector<unsigned short>, 3> vIndex;
    unsigned short indx, ipl;
    for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if (vtx[ivx].NClusters == 0) continue;
      geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
      if (iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
      ipl = iplID.Plane;
      if (ipl > 2) continue;
      indx = vIndex[ipl].size();
      vIndex[ipl].resize(indx + 1);
      vIndex[ipl][indx] = ivx;
    }

    // vector of 2D vertices -> 3D vertices.
    std::vector<short> vPtr;
    for (unsigned short ii = 0; ii < vtx.size(); ++ii)
      vPtr.push_back(-1);

    // temp vector of all 2D vertex matches
    std::vector<Vtx3Store> v3temp;

    double y = 0, z = 0;
    geo::Point_t WPos{0, 0, 0};
    // i, j, k indicates 3 different wire planes
    unsigned short ii, jpl, jj, kpl, kk, ivx, jvx, kvx;
    unsigned int iWire, jWire;
    unsigned short v3dBest = 0;
    float xbest = 0, ybest = 0, zbest = 0;
    float kX, kWire;
    // compare vertices in each view
    bool gotit = false;
    for (ipl = 0; ipl < 2; ++ipl) {
      geo::PlaneID const plane_i{tpcid, ipl};
      for (ii = 0; ii < vIndex[ipl].size(); ++ii) {
        ivx = vIndex[ipl][ii];
        if (ivx > vtx.size() - 1) {
          mf::LogError("CC") << "VtxMatch: bad ivx " << ivx;
          return;
        }
        // vertex has been matched already
        if (vPtr[ivx] >= 0) continue;
        iWire = vtx[ivx].Wire;
        float best = fVertex3DCut;
        // temp array of 2D vertex indices in each plane
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        std::array<short, 3> t2dIndex = {{-1, -1, -1}};
        std::array<short, 3> tmpIndex = {{-1, -1, -1}};
        for (jpl = ipl + 1; jpl < 3; ++jpl) {
          geo::PlaneID const plane_j{tpcid, jpl};
          for (jj = 0; jj < vIndex[jpl].size(); ++jj) {
            jvx = vIndex[jpl][jj];
            if (jvx > vtx.size() - 1) {
              mf::LogError("CC") << "VtxMatch: bad jvx " << jvx;
              return;
            }
            // vertex has been matched already
            if (vPtr[jvx] >= 0) continue;
            jWire = vtx[jvx].Wire;
            // new stuff
            float dX = fabs(vX[ivx] - vX[jvx]);
            float dXSigma = sqrt(vXsigma[ivx] * vXsigma[ivx] + vXsigma[jvx] * vXsigma[jvx]);
            float dXChi = dX / dXSigma;

            if (vtxprt)
              mf::LogVerbatim("CC")
                << "VtxMatch: ipl " << ipl << " ivx " << ivx << " ivX " << vX[ivx] << " jpl " << jpl
                << " jvx " << jvx << " jvX " << vX[jvx] << " W:T " << (int)vtx[jvx].Wire << ":"
                << (int)vtx[jvx].Time << " dXChi " << dXChi << " fVertex3DCut " << fVertex3DCut;

            if (dXChi > fVertex3DCut) continue;
            geom->IntersectionPoint(geo::WireID{plane_i, iWire}, geo::WireID{plane_j, jWire}, y, z);
            if (y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
            WPos.SetY(y);
            WPos.SetZ(z);
            kpl = 3 - ipl - jpl;
            kX = 0.5 * (vX[ivx] + vX[jvx]);
            kWire = -1;
            if (TPC.Nplanes() > 2) {
              try {
                kWire = geom->NearestWireID(WPos, geo::PlaneID{cstat, tpc, kpl}).Wire;
              }
              catch (geo::InvalidWireError const& e) {
                kWire = e.suggestedWireID().Wire; // pick the closest valid wire
              }
            }
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
            float yzSigma = wirePitch * sqrt(vtx[ivx].WireErr * vtx[ivx].WireErr +
                                             vtx[jvx].WireErr * vtx[jvx].WireErr);
            v3d.YErr = yzSigma;
            v3d.Z = z;
            v3d.ZErr = yzSigma;
            v3d.Wire = kWire;
            v3d.CStat = cstat;
            v3d.TPC = tpc;
            v3temp.push_back(v3d);

            if (vtxprt)
              mf::LogVerbatim("CC")
                << "VtxMatch: 2 Plane match ivx " << ivx << " P:W:T " << ipl << ":"
                << (int)vtx[ivx].Wire << ":" << (int)vtx[ivx].Time << " jvx " << jvx << " P:W:T "
                << jpl << ":" << (int)vtx[jvx].Wire << ":" << (int)vtx[jvx].Time << " dXChi "
                << dXChi << " yzSigma " << yzSigma;

            if (TPC.Nplanes() == 2) continue;
            // look for a 3 plane match
            best = fVertex3DCut;
            for (kk = 0; kk < vIndex[kpl].size(); ++kk) {
              kvx = vIndex[kpl][kk];
              if (vPtr[kvx] >= 0) continue;
              // Wire difference error
              float dW = wirePitch * (vtx[kvx].Wire - kWire) / yzSigma;
              // X difference error
              float dX = (vX[kvx] - kX) / dXSigma;
              float kChi = 0.5 * sqrt(dW * dW + dX * dX);
              if (kChi < best) {
                best = kChi;
                xbest = (vX[kvx] + 2 * kX) / 3;
                ybest = y;
                zbest = z;
                t2dIndex[ipl] = ivx;
                t2dIndex[jpl] = jvx;
                t2dIndex[kpl] = kvx;
                v3dBest = v3temp.size() - 1;
              }

              if (vtxprt)
                mf::LogVerbatim("CC")
                  << " kvx " << kvx << " kpl " << kpl << " wire " << (int)vtx[kvx].Wire << " kTime "
                  << (int)vtx[kvx].Time << " kChi " << kChi << " best " << best << " dW "
                  << vtx[kvx].Wire - kWire;

            } // kk
            if (vtxprt)
              mf::LogVerbatim("CC") << " done best = " << best << " fVertex3DCut " << fVertex3DCut;
            if (TPC.Nplanes() > 2 && best < fVertex3DCut) {
              // create a real 3D vertex using the previously entered incomplete 3D vertex as a template
              if (v3dBest > v3temp.size() - 1) {
                mf::LogError("CC") << "VtxMatch: bad v3dBest " << v3dBest;
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
              for (unsigned short jj = 0; jj < 3; ++jj)
                if (t2dIndex[jj] >= 0) vPtr[t2dIndex[jj]] = vtx3.size() - 1;

              if (vtxprt)
                mf::LogVerbatim("CC")
                  << "New 3D vtx " << vtx3.size() << " X " << v3d.X << " Y " << v3d.Y << " Z "
                  << v3d.Z << " t2dIndex " << t2dIndex[0] << " " << t2dIndex[1] << " "
                  << t2dIndex[2] << " best Chi " << best;

            } // best < dRCut
            if (gotit) break;
          } // jj
          if (gotit) break;
        } // jpl
        if (gotit) break;
      } // ii
    }   // ipl

    // Store incomplete 3D vertices but ignore those that are part of a complete 3D vertex
    unsigned short vsize = vtx3.size();
    for (unsigned short it = 0; it < v3temp.size(); ++it) {
      bool keepit = true;
      for (unsigned short i3d = 0; i3d < vsize; ++i3d) {
        for (unsigned short plane = 0; plane < 3; ++plane) {
          if (v3temp[it].Ptr2D[plane] == vtx3[i3d].Ptr2D[plane]) {
            keepit = false;
            break;
          }
        } // plane
        if (!keepit) break;
      } // i3d

      if (keepit) vtx3.push_back(v3temp[it]);
    } // it

    // Modify Ptr2D for 2-plane detector
    if (TPC.Nplanes() == 2) {
      for (unsigned short iv3 = 0; iv3 < vtx3.size(); ++iv3) {
        vtx3[iv3].Ptr2D[2] = 666;
      } //iv3
    }   // 2 planes

    if (vtxprt) {
      for (unsigned short it = 0; it < vtx3.size(); ++it) {
        mf::LogVerbatim("CC") << "vtx3 " << it << " Ptr2D " << vtx3[it].Ptr2D[0] << " "
                              << vtx3[it].Ptr2D[1] << " " << vtx3[it].Ptr2D[2] << " wire "
                              << vtx3[it].Wire;
      }
    }

  } // VtxMatch

  //////////////////////////////////
  void ClusterCrawlerAlg::GetHitRange(CTP_t CTP)
  {
    // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
    // Hits must have been sorted by increasing wire number
    fFirstHit = 0;
    geo::PlaneID planeID = DecodeCTP(CTP);
    unsigned int nwires = geom->Nwires(planeID);
    WireHitRange.resize(nwires + 1);

    // These will be re-defined later
    fFirstWire = 0;
    fLastWire = 0;

    unsigned int wire, iht;
    unsigned int nHitInPlane;
    std::pair<int, int> flag;

    // Define the "no hits on wire" condition
    flag.first = -2;
    flag.second = -2;
    for (auto& apair : WireHitRange)
      apair = flag;

    nHitInPlane = 0;

    std::vector<bool> firsthit;
    firsthit.resize(nwires + 1, true);
    bool firstwire = true;
    for (iht = 0; iht < fHits.size(); ++iht) {
      if (fHits[iht].WireID().asPlaneID() != planeID) continue;
      wire = fHits[iht].WireID().Wire;
      // define the first hit start index in this TPC, Plane
      if (firsthit[wire]) {
        WireHitRange[wire].first = iht;
        firsthit[wire] = false;
      }
      if (firstwire) {
        fFirstWire = wire;
        firstwire = false;
      }
      WireHitRange[wire].second = iht + 1;
      fLastWire = wire + 1;
      ++nHitInPlane;
    }
    // overwrite with the "dead wires" condition
    lariov::ChannelStatusProvider const& channelStatus =
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    flag.first = -1;
    flag.second = -1;
    unsigned int nbad = 0;
    for (wire = 0; wire < nwires; ++wire) {
      raw::ChannelID_t chan = geom->PlaneWireToChannel(geo::WireID(planeID, wire));
      if (!channelStatus.IsGood(chan)) {
        WireHitRange[wire] = flag;
        ++nbad;
      }
    } // wire
    // define the MergeAvailable vector and check for errors
    if (mergeAvailable.size() < fHits.size())
      throw art::Exception(art::errors::LogicError)
        << "GetHitRange: Invalid mergeAvailable vector size " << mergeAvailable.size()
        << fHits.size();
    unsigned int firstHit, lastHit;
    unsigned int cnt;
    cnt = 0;
    float maxRMS, chiSep, peakCut;
    for (wire = 0; wire < nwires; ++wire) {
      // ignore dead wires and wires with no hits
      if (WireHitRange[wire].first < 0) continue;
      firstHit = WireHitRange[wire].first;
      lastHit = WireHitRange[wire].second;
      for (iht = firstHit; iht < lastHit; ++iht) {
        if (fHits[iht].WireID().Wire != wire)
          throw art::Exception(art::errors::LogicError)
            << "Bad WireHitRange on wire " << wire << "\n";
        ++cnt;
        if (fHits[iht].Multiplicity() > 1) {
          peakCut = 0.6 * fHits[iht].PeakAmplitude();
          std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(iht);
          for (size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
            if (jht == iht) continue;
            // require that the j hit be similar in magnitude to the i hit
            if (fHits[jht].PeakAmplitude() < peakCut) continue;
            maxRMS = std::max(fHits[iht].RMS(), fHits[jht].RMS());
            chiSep = std::abs(fHits[iht].PeakTime() - fHits[jht].PeakTime()) / maxRMS;
            if (chiSep < fHitMergeChiCut) {
              mergeAvailable[iht] = true;
              break;
            }
          } // jht
        }   // fHits[iht].Multiplicity() > 1
      }     // iht
    }       // wire
    if (cnt != nHitInPlane)
      mf::LogWarning("CC") << "Bad WireHitRange count " << cnt << " " << nHitInPlane << "\n";

    if (!fMergeAllHits) return;

    // Merge all of the hits
    bool didMerge;
    for (wire = 0; wire < nwires; ++wire) {
      if (WireHitRange[wire].first < 0) continue;
      firstHit = WireHitRange[wire].first;
      lastHit = WireHitRange[wire].second;
      for (iht = firstHit; iht < lastHit; ++iht) {
        if (!mergeAvailable[iht]) continue;
        // already merged?
        if (fHits[iht].GoodnessOfFit() == 6666) continue;
        MergeHits(iht, didMerge);
        mergeAvailable[iht] = false;
      } // iht
    }   // wire

  } // GetHitRange()

  //////////////////////////////////////////
  unsigned int ClusterCrawlerAlg::DeadWireCount(unsigned int inWire1, unsigned int inWire2)
  {
    if (inWire1 > inWire2) {
      // put in increasing order
      unsigned int tmp = inWire1;
      inWire1 = inWire2;
      inWire2 = tmp;
    } // inWire1 > inWire2
    ++inWire2;
    unsigned int wire, ndead = 0;
    for (wire = inWire1; wire < inWire2; ++wire)
      if (WireHitRange[wire].first == -1) ++ndead;
    return ndead;
  } // DeadWireCount

  //////////////////////////////////////////
  unsigned int ClusterCrawlerAlg::DeadWireCount()
  {
    // Counts the number of dead wires in the range spanned by fcl2hits
    if (fcl2hits.size() < 2) return 0;
    unsigned int wire, ndead = 0;
    unsigned int iht = fcl2hits[fcl2hits.size() - 1];
    unsigned int eWire = fHits[iht].WireID().Wire;
    iht = fcl2hits[0];
    unsigned int bWire = fHits[iht].WireID().Wire + 1;
    for (wire = eWire; wire < bWire; ++wire)
      if (WireHitRange[wire].first == -1) ++ndead;
    return ndead;
  } // DeadWireCount

  //////////////////////////////////////////
  bool ClusterCrawlerAlg::areInSameMultiplet(recob::Hit const& first_hit,
                                             recob::Hit const& second_hit)
  {
    return (first_hit.StartTick() == second_hit.StartTick()) &&
           (first_hit.WireID() == second_hit.WireID());
  } // ClusterCrawlerAlg::areInSameMultiplet()

  //////////////////////////////////////////
  std::pair<size_t, size_t> ClusterCrawlerAlg::FindHitMultiplet(size_t iHit) const
  {
    std::pair<size_t, size_t> range{iHit, iHit};

    range.first = iHit - fHits[iHit].LocalIndex();
    range.second = range.first + fHits[iHit].Multiplicity();

    return range;
  } // ClusterCrawlerAlg::FindHitMultiplet()

  //////////////////////////////////////////
  bool ClusterCrawlerAlg::CheckHitDuplicates(std::string location,
                                             std::string marker /* = "" */) const
  {
    // currently unused, only for debug
    unsigned int nDuplicates = 0;
    std::set<unsigned int> hits;
    for (unsigned int hit : fcl2hits) {
      if (hits.count(hit)) {
        ++nDuplicates;
        mf::LogProblem log("CC");
        log << "Hit #" << hit
            << " being included twice in the future cluster (ID=" << (tcl.size() + 1)
            << "?) at location: " << location;
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

    if (icl > (short)tcl.size()) return 9999;

    float cwire, cslp, ctick;
    // figure out which cluster to use
    if (icl < 0) {
      if (fcl2hits.size() == 0) return 9999;
      // cluster under construction
      if (end == 0) {
        cwire = clBeginWir;
        cslp = clBeginSlp;
        ctick = clBeginTim;
      }
      else {
        cwire = clpar[2];
        cslp = clpar[1];
        ctick = clpar[0];
      } // end
    }
    else {
      // tcl cluster
      if (end == 0) {
        cwire = tcl[icl].BeginWir;
        cslp = tcl[icl].BeginSlp;
        ctick = tcl[icl].BeginTim;
      }
      else {
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

    if (icl > (short)tcl.size()) return 9999;
    if (ivx > vtx.size()) return 9999;

    float cwire, cslp, cslpErr, ctick;
    // figure out which cluster to use
    if (icl < 0) {
      if (fcl2hits.size() == 0) return 9999;
      // cluster under construction
      if (end == 0) {
        cwire = clBeginWir;
        cslp = clBeginSlp;
        cslpErr = clBeginSlpErr;
        ctick = clBeginTim;
      }
      else {
        cwire = clpar[2];
        cslp = clpar[1];
        cslpErr = clparerr[1];
        ctick = clpar[0];
      } // end
    }
    else {
      // tcl cluster
      if (end == 0) {
        cwire = tcl[icl].BeginWir;
        cslp = tcl[icl].BeginSlp;
        cslpErr = tcl[icl].BeginSlpErr;
        ctick = tcl[icl].BeginTim;
      }
      else {
        cwire = tcl[icl].EndWir;
        cslp = tcl[icl].EndSlp;
        cslpErr = tcl[icl].EndSlpErr;
        ctick = tcl[icl].EndTim;
      } // end
    }

    // Closest approach wire
    float docaW =
      (vtx[ivx].Wire + cslp * (vtx[ivx].Time - ctick) + cwire * cslp * cslp) / (1 + cslp * cslp);
    float dW = docaW - vtx[ivx].Wire;
    float chi = dW / vtx[ivx].WireErr;
    float totChi = chi * chi;
    dW = vtx[ivx].Wire - cwire;
    float dT = ctick + dW * cslp - vtx[ivx].Time;
    if (cslpErr < 1E-3) cslpErr = 1E-3;
    // increase slope error for large angle clusters
    cslpErr *= AngleFactor(cslp);
    // cluster slope projection error
    float dTErr = dW * cslpErr;
    // squared
    dTErr *= dTErr;
    // add the vertex time error^2 to the cluster projection error^2
    dTErr += vtx[ivx].TimeErr * vtx[ivx].TimeErr;
    if (dTErr < 1E-3) dTErr = 1E-3;
    totChi += dT * dT / dTErr;
    totChi /= 2;

    return totChi;

  } // ClusterVertexChi

  /////////////////////////////////////////
  float ClusterCrawlerAlg::PointVertexChi(float wire, float tick, unsigned short ivx)
  {
    // Returns the Chisq/DOF between a (Wire, Tick) point and a vertex

    if (ivx > vtx.size()) return 9999;

    float dW = wire - vtx[ivx].Wire;
    float chi = dW / vtx[ivx].WireErr;
    float totChi = chi * chi;
    float dT = tick - vtx[ivx].Time;
    chi = dT / vtx[ivx].TimeErr;
    totChi += chi * chi;

    return totChi;

  } // PointVertexChi

  /////////////////////////////////////////
  std::string ClusterCrawlerAlg::PrintHit(unsigned int iht)
  {

    if (iht > fHits.size() - 1) return "Bad Hit";
    return std::to_string(fHits[iht].WireID().Plane) + ":" +
           std::to_string(fHits[iht].WireID().Wire) + ":" +
           std::to_string((int)fHits[iht].PeakTime());

  } // PrintHit

} // namespace cluster
