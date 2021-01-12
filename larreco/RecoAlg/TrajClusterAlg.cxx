//////////////////////////////////////////////////////////////////////
///
/// TrajClusterAlg
///
/// Bruce Baller, baller@fnal.gov
/// Citation: Liquid argon TPC signal formation, signal processing and reconstruction techniques
/// B. Baller 2017 JINST 12 P07010
///
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrajClusterAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
#include "larreco/RecoAlg/TCAlg/StepUtils.h"
#include "larreco/RecoAlg/TCAlg/PostStepUtils.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace tca {

  //------------------------------------------------------------------------------

  TrajClusterAlg::TrajClusterAlg(fhicl::ParameterSet const& pset)
    : fCaloAlg(pset.get<fhicl::ParameterSet>("CaloAlg")), fMVAReader("Silent")
  {
    tcc.showerParentReader = &fMVAReader;

    bool badinput = false;
    // set all configurable modes false
    tcc.modes.reset();

    if (pset.has_key("UseChannelStatus")) tcc.useChannelStatus = pset.get<bool>("UseChannelStatus");
    std::vector<std::string> skipAlgsVec;
    if (pset.has_key("SkipAlgs")) skipAlgsVec = pset.get<std::vector<std::string>>("SkipAlgs");
    std::vector<std::string> debugConfigVec;
    if (pset.has_key("DebugConfig"))
      debugConfigVec = pset.get<std::vector<std::string>>("DebugConfig");

    tcc.hitErrFac = pset.get<float>("HitErrFac", 0.4);
    // Allow the user to specify the typical hit rms for small-angle tracks
    std::vector<float> aveHitRMS;
    if (pset.has_key("AveHitRMS")) aveHitRMS = pset.get<std::vector<float>>("AveHitRMS");
    // Turn off the call to AnalyzeHits
    if (!aveHitRMS.empty()) {
      evt.aveHitRMSValid = true;
      evt.aveHitRMS = aveHitRMS;
    }
    tcc.angleRanges = pset.get<std::vector<float>>("AngleRanges");
    tcc.nPtsAve = pset.get<short>("NPtsAve", 20);
    tcc.minPtsFit = pset.get<std::vector<unsigned short>>("MinPtsFit");
    tcc.minPts = pset.get<std::vector<unsigned short>>("MinPts");
    tcc.maxAngleCode = pset.get<std::vector<unsigned short>>("MaxAngleCode");
    tcc.maxChi = pset.get<float>("MaxChi", 10);
    tcc.chargeCuts = pset.get<std::vector<float>>("ChargeCuts", {3, 0.15, 0.25});
    tcc.multHitSep = pset.get<float>("MultHitSep", 2.5);
    tcc.kinkCuts = pset.get<std::vector<float>>("KinkCuts");
    tcc.qualityCuts = pset.get<std::vector<float>>("QualityCuts", {0.8, 3});
    tcc.maxWireSkipNoSignal = pset.get<float>("MaxWireSkipNoSignal", 1);
    tcc.maxWireSkipWithSignal = pset.get<float>("MaxWireSkipWithSignal", 100);
    tcc.projectionErrFactor = pset.get<float>("ProjectionErrFactor", 2);
    tcc.VLAStepSize = pset.get<float>("VLAStepSize", 1.5);
    tcc.JTMaxHitSep2 = pset.get<float>("JTMaxHitSep", 2);
    tcc.deltaRayTag = pset.get<std::vector<short>>("DeltaRayTag", {-1, -1, -1});
    tcc.muonTag = pset.get<std::vector<short>>("MuonTag", {-1, -1, -1, -1});
    if (pset.has_key("ElectronTag")) tcc.electronTag = pset.get<std::vector<float>>("ElectronTag");
    tcc.showerTag = pset.get<std::vector<float>>("ShowerTag", {-1, -1, -1, -1, -1, -1});
    std::string fMVAShowerParentWeights = "NA";
    pset.get_if_present<std::string>("MVAShowerParentWeights", fMVAShowerParentWeights);
    tcc.chkStopCuts = pset.get<std::vector<float>>("ChkStopCuts", {-1, -1, -1});
    tcc.vtx2DCuts = pset.get<std::vector<float>>("Vertex2DCuts", {-1, -1, -1, -1, -1, -1, -1});
    tcc.vtx3DCuts = pset.get<std::vector<float>>("Vertex3DCuts", {-1, -1});
    tcc.vtxScoreWeights = pset.get<std::vector<float>>("VertexScoreWeights");
    tcc.match3DCuts = pset.get<std::vector<float>>("Match3DCuts", {-1, -1, -1, -1, -1});
    tcc.pfpStitchCuts = pset.get<std::vector<float>>("PFPStitchCuts", {-1});
    // Set the normal stepping mode Pos (lower wire number to higher wire number)
    tcc.modes[kModeStepPos] = true;
    // don't produce neutrino PFParticles, etc unless desired
    tcc.modes[kModeNeutrino] = pset.get<bool>("NeutrinoMode", false);
    pset.get_if_present<std::vector<float>>("NeutralVxCuts", tcc.neutralVxCuts);
    if (tcc.JTMaxHitSep2 > 0) tcc.JTMaxHitSep2 *= tcc.JTMaxHitSep2;

    // in the following section we ensure that the fcl vectors are appropriately sized so that later references are valid
    if (tcc.minPtsFit.size() != tcc.minPts.size()) badinput = true;
    if (tcc.maxAngleCode.size() != tcc.minPts.size()) badinput = true;
    if (badinput)
      throw art::Exception(art::errors::Configuration)
        << "Bad input from fcl file. Vector lengths for MinPtsFit and MaxAngleRange "
           "should be defined for each reconstruction pass";

    if (tcc.vtx2DCuts.size() < 10)
      throw art::Exception(art::errors::Configuration)
        << "Vertex2DCuts must be size 11\n 0 = Max length definition for short TJs\n 1 = Max "
           "vtx-TJ sep short TJs\n 2 = Max vtx-TJ sep long TJs\n 3 = Max position pull for >2 "
           "TJs\n 4 = Max vtx position error\n 5 = Min MCSMom for one of two TJs\n 6 = Min "
           "fraction of wires hit btw vtx and Tjs\n 7 = Min Score\n 8 = min ChgFrac at a vtx or "
           "merge point\n 9 = max MCSMom asymmetry, 10 = require chg on wires btw vtx and tjs in "
           "induction planes?";
    if (tcc.vtx2DCuts.size() == 10) {
      // User didn't specify a requirement for the presence of charge between a vertex and the start of the
      // vertex Tjs in induction planes. Assume that it is required
      tcc.vtx2DCuts.resize(11, 1.);
    }
    if (tcc.vtx3DCuts.size() < 3)
      throw art::Exception(art::errors::Configuration)
        << "Vertex3DCuts must be size > 2\n 0 = 2D Vtx max dX (cm)\n 1 = 2D Vtx max pull\n 2 = max "
           "3D separation (cm) btw PFP and vertex";
    if (tcc.vtx3DCuts.size() == 2) {
      std::cout << "WARNING: Vertex3DCuts is size 2 but should be size 3, where Vertex3DCuts[2] = "
                   "max 3D separation (cm) btw a PFP and a 3D vertex. Setting it to 3 cm\n";
      tcc.vtx3DCuts.resize(3, 3.);
    }
    if (tcc.kinkCuts.size() < 3) {
      throw art::Exception(art::errors::Configuration)
        << "KinkCuts must be size 3\n 0 = Number of points to fit at the end of the trajectory\n 1 "
           "= Minimum kink significance\n 2 = Use charge in significance calculation? (yes if > "
           "0). \nYou are using an out-of-date specification?\n";
    }

    if (tcc.chargeCuts.size() != 3)
      throw art::Exception(art::errors::Configuration)
        << "ChargeCuts must be size 3\n 0 = Charge pull cut\n 1 = Min allowed fractional chg RMS\n "
           "2 = Max allowed fractional chg RMS";
    // dressed muons - change next line
    if (tcc.muonTag.size() < 4)
      throw art::Exception(art::errors::Configuration)
        << "MuonTag must be size 4\n 0 = minPtsFit\n 1 = minMCSMom\n 2= maxWireSkipNoSignal\n 3 = "
           "min delta ray length for tagging\n 4 = dress muon window size (optional)";
    if (tcc.deltaRayTag.size() != 3)
      throw art::Exception(art::errors::Configuration)
        << "DeltaRayTag must be size 3\n 0 = Max endpoint sep\n 1 = min MCSMom\n 2 = max MCSMom";
    if (tcc.chkStopCuts.size() < 3)
      throw art::Exception(art::errors::Configuration)
        << "ChkStopCuts must be size 3\n 0 = Min Charge ratio\n 1 = Charge slope pull cut\n 2 = "
           "Charge fit chisq cut";
    if(tcc.showerTag.size() != 3) {
      std::cout
        <<  "Shower reconstruction code is not available. Reconfiguring ShowerTag to size 3:\n"
            " [0] = Min MCSMom, default = 100\n"
            " [1] = Max trajectory separation in WSE units, default = 2 units\n"
            " [2] = Min number of trajectories for a ShowerLike tag, default = 3\n";
      tcc.showerTag.resize(3);
      tcc.showerTag[0] = 100;
      tcc.showerTag[1] = 2;
      tcc.showerTag[2] = 3;
    } // tcc.showerTag.size() != 3
    if (tcc.match3DCuts.size() < 6)
      throw art::Exception(art::errors::Configuration)
        << "Match3DCuts must be size 7\n 0 = dx (cm) matching cut\n 1 = max number of 3D "
           "combinations\n 2 = min length for 2-view match\n 3 = number of TP3Ds in each plane to "
           "fit in each PFP section\n 4 = max pull for accepting TP3Ds in sections\n 5 = max "
           "ChiDOF for a SectionFit\n 6 = match limit for shower-like Tjs";
    if (tcc.match3DCuts.size() == 6) {
      // add another cut on limit the number of shower-like tjs to 3D match
      tcc.match3DCuts.resize(7);
      tcc.match3DCuts[6] = 1000;
    } 
    // check the angle ranges and convert from degrees to radians
    if (tcc.angleRanges.back() < 90) {
      mf::LogVerbatim("TC") << "Last element of AngleRange != 90 degrees. Fixing it\n";
      tcc.angleRanges.back() = 90;
    }

    // convert PFP stitch cuts
    if (tcc.pfpStitchCuts.size() > 1 && tcc.pfpStitchCuts[0] > 0) {
      // square the separation cut
      tcc.pfpStitchCuts[0] *= tcc.pfpStitchCuts[0];
      // convert angle to cos
      tcc.pfpStitchCuts[1] = cos(tcc.pfpStitchCuts[1]);
    }

    // configure algorithm debugging. Configuration for debugging standard stepping
    // is done in Utils/AnalyzeHits when the input hit collection is passed to SetInputHits
    tcc.modes[kModeDebug] = false;
    tcc.dbgAlg.reset();
    for (auto strng : debugConfigVec) {
      // try to interpret this as a C:T:P:W:Tick specification or something similar
      if (!DecodeDebugString(strng)) {
        // try to set a dbgAlg bit
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
          if (strng == AlgBitNames[ib]) {
            tcc.dbgAlg[ib] = true;
            tcc.modes[kModeDebug] = true;
            break;
          }
        } // ib
        // print some instructions and quit if there was a failure
        if (!tcc.modes[kModeDebug]) {
          std::cout << "DecodeDebugString failed: " << strng << "\n";
          DecodeDebugString("instruct");
          exit(1);
        }
      } // DecodeDebugString failed
    }   // strng

    for (auto& range : tcc.angleRanges) {
      if (range < 0 || range > 90)
        throw art::Exception(art::errors::Configuration)
          << "Invalid angle range " << range << " Must be 0 - 90 degrees";
      range *= M_PI / 180;
    } // range

    // Ensure that the size of the AlgBitNames vector is consistent with the AlgBit typedef enum
    if (kAlgBitSize != AlgBitNames.size())
      throw art::Exception(art::errors::Configuration)
        << "kAlgBitSize " << kAlgBitSize << " != AlgBitNames size " << AlgBitNames.size();
    if (kAlgBitSize > 128)
      throw art::Exception(art::errors::Configuration)
        << "Increase the size of UseAlgs to at least " << kAlgBitSize;
    fAlgModCount.resize(kAlgBitSize);

    if (kFlagBitSize != EndFlagNames.size())
      throw art::Exception(art::errors::Configuration)
        << "kFlagBitSize " << kFlagBitSize << " != EndFlagNames size " << EndFlagNames.size();

    if (kFlagBitSize > 8)
      throw art::Exception(art::errors::Configuration)
        << "Increase the size of EndFlag to at least " << kFlagBitSize;

    bool printHelp = false;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
      tcc.useAlg[ib] = true;

    for (auto strng : skipAlgsVec) {
      bool gotit = false;
      if (strng == "All") {
        // turn everything off
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          tcc.useAlg[ib] = false;
        gotit = true;
        break;
      } // All off
      for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
        if (strng == AlgBitNames[ib]) {
          tcc.useAlg[ib] = false;
          gotit = true;
          break;
        }
      } // ib
      if (!gotit) {
        std::cout << "******* Unknown SkipAlgs input string '" << strng << "'\n";
        printHelp = true;
      }
    } // strng
    if (printHelp) {
      std::cout << "Valid AlgNames:";
      for (auto strng : AlgBitNames)
        std::cout << " " << strng;
      std::cout << "\n";
      std::cout << "Or specify All to turn all algs off\n";
    }
    evt.eventsProcessed = 0;

    tcc.caloAlg = &fCaloAlg;
  }

  ////////////////////////////////////////////////
  bool
  TrajClusterAlg::SetInputHits(std::vector<recob::Hit> const& inputHits,
                               unsigned int run,
                               unsigned int event)
  {
    // defines the pointer to the input hit collection, analyzes them,
    // initializes global counters and refreshes service references
    ClearResults();
    evt.allHits = &inputHits;
    evt.run = run;
    evt.event = event;
    // refresh service references
    tcc.geom = lar::providerFrom<geo::Geometry>();
    evt.WorkID = 0;
    evt.globalT_UID = 0;
    evt.global2V_UID = 0;
    evt.global3V_UID = 0;
    evt.globalP_UID = 0;
    evt.global2S_UID = 0;
    evt.global3S_UID = 0;
    // find the average hit RMS using the full hit collection and define the
    // configuration for the current TPC

    if(tcc.modes[kModeDebug] && evt.eventsProcessed == 0) PrintDebugMode();
    ++evt.eventsProcessed;

    return AnalyzeHits();
  } // SetInputHits

  ////////////////////////////////////////////////
  void
  TrajClusterAlg::SetSourceHits(std::vector<recob::Hit> const& srcHits)
  {
    evt.srcHits = &srcHits;
    evt.tpcSrcHitRange.resize(tcc.geom->NTPC());
    for (auto& thr : evt.tpcSrcHitRange)
      thr = {UINT_MAX, UINT_MAX};
    for (unsigned int iht = 0; iht < (*evt.srcHits).size(); ++iht) {
      auto& hit = (*evt.srcHits)[iht];
      unsigned int tpc = hit.WireID().TPC;
      if (tpc >= evt.tpcSrcHitRange.size()) return;
      if (evt.tpcSrcHitRange[tpc].first == UINT_MAX) evt.tpcSrcHitRange[tpc].first = iht;
      evt.tpcSrcHitRange[tpc].second = iht;
    } // iht
  }   // SetSourceHits

  ////////////////////////////////////////////////
  void
  TrajClusterAlg::RunTrajClusterAlg(detinfo::DetectorClocksData const& clockData,
                                    detinfo::DetectorPropertiesData const& detProp,
                                    std::vector<unsigned int>& hitsInSlice,
                                    int sliceID)
  {
    // Reconstruct everything using the hits in a slice

    if(hitsInSlice.size() < 2) return;
    if(tcc.recoSlice > 0 && sliceID != tcc.recoSlice) return;

    if (!CreateSlice(clockData, detProp, hitsInSlice, sliceID)) return;

    seeds.resize(0);
    // get a reference to the stored slice
    auto& slc = slices[slices.size() - 1];
    // special debug mode reconstruction
    if (tcc.recoTPC > 0 && (short)slc.TPCID.TPC != tcc.recoTPC) {
      slices.pop_back();
      return;
    }

    if (evt.aveHitRMS.size() != slc.nPlanes)
      throw art::Exception(art::errors::Configuration)
        << " AveHitRMS vector size != the number of planes ";
    if (tcc.recoSlice)
      std::cout << "Reconstruct " << hitsInSlice.size() << " hits in Slice " << sliceID
                << " in TPC " << slc.TPCID.TPC << "\n";
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      ReconstructAllTraj(detProp, slc, inCTP);
      if (!slc.isValid) return;
    } // plane
    // Compare 2D vertices in each plane and try to reconcile T -> 2V attachments using
    // 2D and 3D(?) information
    Reconcile2Vs(slc);
    Find3Vs(detProp, slc);
    KillOrphan2Vs(detProp, slc);
    ScoreVertices(slc);
    // Define the ParentID of trajectories using the vertex score
//    DefineTjParents(slc, false);
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      if (!ChkVtxAssociations(slc, inCTP)) {
        std::cout << "RTC: ChkVtxAssociations found an error\n";
      }
    } // plane
    if (tcc.match3DCuts[0] > 0) {
      FindPFParticles(clockData, detProp, slc);
    } // 3D matching requested
    KillPoorVertices(slc);
    if (!slc.isValid) {
      mf::LogVerbatim("TC") << "RunTrajCluster failed in MakeAllTrajClusters";
      return;
    }

    // dump a trajectory?
    if (tcc.modes[kModeDebug] && tcc.dbgDump) DumpTj();

    // count algorithm usage
    for (auto& tj : slc.tjs) {
      for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
        if (tj.AlgMod[ib]) ++fAlgModCount[ib];
    } // tj

    // clear vectors that are not needed later
    slc.mallTraj.resize(0);

  } // RunTrajClusterAlg

  ////////////////////////////////////////////////
  void
  TrajClusterAlg::ReconstructAllTraj(detinfo::DetectorPropertiesData const& detProp,
                                     TCSlice& slc,
                                     CTP_t inCTP)
  {
    // Reconstruct trajectories in inCTP and put them in allTraj

    unsigned short plane = DecodeCTP(inCTP).Plane;
    if (slc.firstWire[plane] > slc.nWires[plane]) return;
    unsigned int nwires = slc.lastWire[plane] - slc.firstWire[plane] - 1;
    if (nwires > slc.nWires[plane]) return;

    // Make several passes through the hits with user-specified cuts for each
    // pass. In general these are to not reconstruct large angle trajectories on
    // the first pass
    float maxHitsRMS = 4 * evt.aveHitRMS[plane];
    for (unsigned short pass = 0; pass < tcc.minPtsFit.size(); ++pass) {
      // don't try to step with VLA hits
      if(tcc.useAlg[kNewCuts] && tcc.maxAngleCode[pass] == 2) break;
      for (unsigned int ii = 0; ii < nwires; ++ii) {
        // decide which way to step given the sign of StepDir
        unsigned int iwire = 0;
        unsigned int jwire = 0;
        if (tcc.modes[kModeStepPos]) {
          // step DS
          iwire = slc.firstWire[plane] + ii;
          jwire = iwire + 1;
        }
        else {
          // step US
          iwire = slc.lastWire[plane] - ii - 1;
          jwire = iwire - 1;
        }
        if (iwire > slc.wireHitRange[plane].size() - 1) continue;
        if (jwire > slc.wireHitRange[plane].size() - 1) continue;
        // skip bad wires or no hits on the wire
        if (slc.wireHitRange[plane][iwire].first == UINT_MAX) continue;
        if (slc.wireHitRange[plane][jwire].first == UINT_MAX) continue;
        unsigned int ifirsthit = slc.wireHitRange[plane][iwire].first;
        unsigned int ilasthit = slc.wireHitRange[plane][iwire].second;
        unsigned int jfirsthit = slc.wireHitRange[plane][jwire].first;
        unsigned int jlasthit = slc.wireHitRange[plane][jwire].second;
        if (ifirsthit > slc.slHits.size() || ilasthit > slc.slHits.size()) {
          std::cout << "RAT: bad hit range " << ifirsthit << " " << ilasthit << " size "
                    << slc.slHits.size() << " inCTP " << inCTP << "\n";
          return;
        }
        for (unsigned int iht = ifirsthit; iht <= ilasthit; ++iht) {
          tcc.dbgStp = (tcc.modes[kModeDebug] && (slc.slHits[iht].allHitsIndex == debug.Hit));
          if (tcc.dbgStp) {
            mf::LogVerbatim("TC") << "+++++++ Pass " << pass << " Found debug hit "
                                  << slices.size() - 1 << ":" << PrintHit(slc.slHits[iht])
                                  << " iht " << iht;
          }
          // Only consider hits that are available
          if (slc.slHits[iht].InTraj != 0) continue;
          // We hope to make a trajectory point at the hit position of iht in WSE units
          // with a direction pointing to jht
          if (slc.slHits[iht].allHitsIndex > (*evt.allHits).size() - 1) {
            std::cout << "RAT: Bad allHitsIndex\n";
            continue;
          }
          auto& iHit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
          if (LongPulseHit(iHit)) continue;
          unsigned int fromWire = iHit.WireID().Wire;
          float fromTick = iHit.PeakTime();
          float iqtot = iHit.Integral();
          float hitsRMSTick = iHit.RMS();
          std::vector<unsigned int> iHitsInMultiplet;
          if (pass == 0) {
            // only use the hit on the first pass
            iHitsInMultiplet.resize(1);
            iHitsInMultiplet[0] = iht;
          }
          else {
            GetHitMultiplet(slc, iht, iHitsInMultiplet, false);
            if (iHitsInMultiplet.empty()) continue;
            // ignore very high multiplicities
            if (iHitsInMultiplet.size() > 4) continue;
            if (tcc.useAlg[kNewCuts] && iHitsInMultiplet.size() > 3) continue;
          }
          if (iHitsInMultiplet.size() > 1) {
            fromTick = HitsPosTick(slc, iHitsInMultiplet, iqtot, kUnusedHits);
            hitsRMSTick = HitsRMSTick(slc, iHitsInMultiplet, kUnusedHits);
          }
          if (hitsRMSTick == 0) continue;
          bool fatIHit = (hitsRMSTick > maxHitsRMS);
          if (tcc.dbgStp)
            mf::LogVerbatim("TC") << " hit RMS " << iHit.RMS() << " BB Multiplicity "
                                  << iHitsInMultiplet.size() << " AveHitRMS[" << plane << "] "
                                  << evt.aveHitRMS[plane] << " HitsRMSTick " << hitsRMSTick
                                  << " fatIHit " << fatIHit;
          for (unsigned int jht = jfirsthit; jht <= jlasthit; ++jht) {
            // Only consider hits that are available
            if (slc.slHits[iht].InTraj != 0) break;
            if (slc.slHits[jht].InTraj != 0) continue;
            // clear out any leftover work inTraj's that weren't cleaned up properly
            for (auto& slHit : slc.slHits) if (slHit.InTraj < 0) slHit.InTraj = 0;
            unsigned int toWire = jwire;
            auto& jHit = (*evt.allHits)[slc.slHits[jht].allHitsIndex];
            if (LongPulseHit(jHit)) continue;
            float toTick = jHit.PeakTime();
            float jqtot = jHit.Integral();
            std::vector<unsigned int> jHitsInMultiplet;
            if (pass == 0) {
              // only use the hit on the first pass
              jHitsInMultiplet.resize(1);
              jHitsInMultiplet[0] = jht;
            }
            else {
              GetHitMultiplet(slc, jht, jHitsInMultiplet, false);
              if (jHitsInMultiplet.empty()) continue;
              // ignore very high multiplicities
              if (jHitsInMultiplet.size() > 4) continue;
              if (tcc.useAlg[kNewCuts] && jHitsInMultiplet.size() > 3) continue;
            }
            hitsRMSTick = HitsRMSTick(slc, jHitsInMultiplet, kUnusedHits);
            if (hitsRMSTick == 0) continue;
            bool fatJHit = (hitsRMSTick > maxHitsRMS);
            if (pass == 0) {
              // require both hits to be consistent
              if ((fatIHit && !fatJHit) || (!fatIHit && fatJHit)) { continue; }
            }
            else {
              // pass > 0
              if (jHitsInMultiplet.size() > 1)
                toTick = HitsPosTick(slc, jHitsInMultiplet, jqtot, kUnusedHits);
            }
            bool hitsOK = TrajHitsOK(slc, iHitsInMultiplet, jHitsInMultiplet);
            // Ensure that the hits StartTick and EndTick have the proper overlap
            if (!hitsOK) continue;
            // start a trajectory with direction from iht -> jht
            Trajectory work;
            if (!StartTraj(slc, work, fromWire, fromTick, toWire, toTick, inCTP, pass)) continue;
            // check for a major failure
            if (!slc.isValid) {
              std::cout << "RAT: StartTraj major failure\n";
              return;
            }
            if (work.Pts.empty()) {
              if (tcc.dbgStp) mf::LogVerbatim("TC") << "ReconstructAllTraj: StartTraj failed";
              continue;
            }
            work.Pts[0].DeltaRMS = tcc.hitErrFac * tcc.unitsPerTick * hitsRMSTick;
            // don't include the charge of iht since it will be low if this
            // is a starting/ending track
            work.AveChg = jqtot;
            // try to add close hits
            bool sigOK;
            AddHits(slc, work, 0, sigOK);
            // check for a major failure
            if (!slc.isValid) {
              std::cout << "RAT: AddHits major failure\n";
              return;
            }
            if (!sigOK || work.Pts[0].Chg == 0) {
              if (tcc.dbgStp) mf::LogVerbatim("TC") << " No hits at initial trajectory point ";
              ReleaseHits(slc, work);
              continue;
            }
            if(tcc.useAlg[kNewCuts] && work.Pts[0].Hits.size() > 2) {
              if (tcc.dbgStp) mf::LogVerbatim("TC") << " Too many hits " 
                    << work.Pts[0].Hits.size() << " at initial trajectory point ";
              ReleaseHits(slc, work);
              continue;
            }
            // move the TP position to the hit position but don't mess with the direction
            work.Pts[0].Pos = work.Pts[0].HitPos;
            // print the header and the first TP
            if (tcc.dbgStp) PrintTrajectory("RAT", slc, work, USHRT_MAX);
            // We can't update the trajectory yet because there is only one TP.
            work.EndPt[0] = 0;
            // now try stepping away
            StepAway(slc, work);
            // check for a major failure
            if (!slc.isValid) {
              std::cout << "RAT: StepAway major failure\n";
              return;
            }
            if (tcc.dbgStp)
              mf::LogVerbatim("TC") << " After first StepAway. IsGood " << work.IsGood;
            // Check the quality of the work trajectory
            CheckTraj(slc, work);
            // check for a major failure
            if (!slc.isValid) {
              std::cout << "RAT: CheckTraj major failure\n";
              return;
            }
            if (tcc.dbgStp)
              mf::LogVerbatim("TC") << "ReconstructAllTraj: After CheckTraj EndPt " << work.EndPt[0]
                                    << "-" << work.EndPt[1] << " IsGood " << work.IsGood 
                                    << " EndFlag[0] " << PackEndFlags(work, 0);
            if (tcc.dbgStp)
              mf::LogVerbatim("TC")
                << "StepAway done: IsGood " << work.IsGood << " NumPtsWithCharge "
                << NumPtsWithCharge(slc, work, true) << " cut " << tcc.minPts[work.Pass];
            // decide if the trajectory is long enough for this pass
            if (!work.IsGood || NumPtsWithCharge(slc, work, true) < tcc.minPts[work.Pass]) {
              if (tcc.dbgStp)
                mf::LogVerbatim("TC")
                  << " xxxxxxx Not enough points " << NumPtsWithCharge(slc, work, true)
                  << " minimum " << tcc.minPts[work.Pass] << " or !IsGood";
              ReleaseHits(slc, work);
              continue;
            }
            if (!StoreTraj(slc, work)) continue;
            if (tcc.dbgStp) PrintTrajectory("RAT", slc, slc.tjs.back(), USHRT_MAX);
            if (tcc.modes[kModeDebug] && !InTrajOK(slc, "RAT")) {
                std::cout << "RAT: InTrajOK major failure T" << slc.tjs.back().ID << "\n";
                return;
            }
            // This seems like the wrong place to do this
            if(!tcc.useAlg[kNewCuts]) ChkBeginChg(slc, slc.tjs.size() - 1);
            break;
          } // jht
        }   // iht
      }     // iwire

      // See if there are any seed trajectory points that were saved before reverse
      // propagation and try to make Tjs from them
      for (auto tp : seeds) {
        unsigned short nAvailable = 0;
        for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if (!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          if (slc.slHits[iht].InTraj == 0) ++nAvailable;
          tcc.dbgStp = (tcc.modes[kModeDebug] && (slc.slHits[iht].allHitsIndex == debug.Hit));
          if (tcc.dbgStp) {
            mf::LogVerbatim("TC") << "+++++++ Seed debug hit " << slices.size() - 1 << ":"
                                  << PrintHit(slc.slHits[iht]) << " iht " << iht;
          }
        } // ii
        if (nAvailable == 0) continue;
        Trajectory work;
        work.ID = evt.WorkID;
        for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if (!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          if (slc.slHits[iht].InTraj == 0) slc.slHits[iht].InTraj = work.ID;
        } // ii
        work.Pass = pass;
        work.StepDir = 1;
        if (tp.Dir[0] < 0) work.StepDir = -1;
        work.CTP = tp.CTP;
        work.ParentID = -1;
        work.Strategy.reset();
        work.Strategy[kNormal] = true;
        // don't allow yet another reverse propagation
        work.AlgMod[kRvPrp] = true;
        work.Pts.push_back(tp);
        StepAway(slc, work);
        CheckTraj(slc, work);
        // check for a major failure
        if (!slc.isValid) {
          std::cout << "RAT: CheckTraj major failure\n";
          return;
        }
        // decide if the trajectory is long enough for this pass
        if (!work.IsGood || NumPtsWithCharge(slc, work, true) < tcc.minPts[work.Pass]) {
          if (tcc.dbgStp)
            mf::LogVerbatim("TC") << " xxxxxxx Not enough points "
                                  << NumPtsWithCharge(slc, work, true) << " minimum "
                                  << tcc.minPts[work.Pass] << " or !IsGood";
          ReleaseHits(slc, work);
          continue;
        }
        if (!StoreTraj(slc, work)) continue;
        if (tcc.dbgStp) {
          auto& tj = slc.tjs[slc.tjs.size() - 1];
          mf::LogVerbatim("TC") << "TRP RAT Stored T" << tj.ID << " using seed TP "
                                << PrintPos(slc, tp);
        }
      } // seed

      seeds.resize(0);

      bool lastPass = (pass == tcc.minPtsFit.size() - 1);
      // don't use lastPass cuts if we will use LastEndMerge
      if (tcc.useAlg[kLastEndMerge]) lastPass = false;
      EndMerge(slc, inCTP, lastPass);
      if (!slc.isValid) return;

      Find2Vs(detProp, slc, inCTP, pass);

    } // pass

    // Last attempt to merge long straight Tjs that failed the EndMerge cuts
    LastEndMerge(slc, inCTP);
    // make junk trajectories using nearby un-assigned hits
    FindJunkTraj(slc, inCTP);
    // Merge junk Tjs with junk Tjs
    MergeJunk(slc, inCTP);
    // Merge short Tjs into junk Tjs
    MergeShortWithJunk(slc, inCTP);
    // dressed muons with halo trajectories
    if (tcc.muonTag.size() > 4 && tcc.muonTag[4] > 0) {
      for (auto& tj : slc.tjs) {
        if (tj.AlgMod[kKilled]) continue;
        if (tj.PDGCode != 13) continue;
        MakeHaloTj(slc, tj, tcc.dbgSlc);
      } // tj
    }   // dressed muons

    // temp check
    for (auto& tj : slc.tjs) {
      if(!tj.IsGood && !tj.AlgMod[kKilled]) {
        std::cout<<"T"<<tj.ID<<" is not good and is not killed";
        tj.AlgMod[kKilled] = true;
      }
    } // tj

    // Tag ShowerLike Tjs
    TagShowerLike(slc, inCTP);
    // Set TP Environment bits
    SetTPEnvironment(slc, inCTP);
    Find2Vs(detProp, slc, inCTP, USHRT_MAX);
    SplitTrajCrossingVertices(slc, inCTP);
    // Make vertices between long Tjs and junk Tjs
    MakeJunkVertices(slc, inCTP);
    // check for a major failure
    if (!slc.isValid) return;

    // last attempt to attach Tjs to vertices
    for (unsigned short ivx = 0; ivx < slc.vtxs.size(); ++ivx) {
      auto& vx2 = slc.vtxs[ivx];
      if (vx2.ID == 0) continue;
      if (vx2.CTP != inCTP) continue;
      AttachAnyTrajToVertex(slc, ivx, false);
    } // ivx

    // Set the kEnvOverlap bit true for all TPs that are close to other
    // trajectories that are close to vertices. The positions of TPs that
    // overlap are biased and shouldn't be used in a vertex fit. Also, these
    // TPs shouldn't be used to calculate dE/dx
    UpdateVxEnvironment(slc);

    // Check the Tj <-> vtx associations and define the vertex quality
    if (!ChkVtxAssociations(slc, inCTP)) {
      std::cout << "RAT: ChkVtxAssociations found an error. Events processed "
                << evt.eventsProcessed << " WorkID " << evt.WorkID << "\n";
    }

  } // ReconstructAllTraj

  //////////////////////////////////////////
  void
  TrajClusterAlg::FindJunkTraj(TCSlice& slc, CTP_t inCTP)
  {
    // Makes junk trajectories using unassigned hits

    if (tcc.JTMaxHitSep2 <= 0) return;
    if (!tcc.useAlg[kJunkTj]) return;
    unsigned short plane = DecodeCTP(inCTP).Plane;
    if ((int)slc.lastWire[plane] - 3 < (int)slc.firstWire[plane]) return;

    // shouldn't have to do this but...
    for (auto& slHit : slc.slHits)
      if (slHit.InTraj < 0) slHit.InTraj = 0;

    bool prt = false;

    std::vector<unsigned int> tHits;
    // Stay well away from the last wire in the plane
    for (unsigned int iwire = slc.firstWire[plane]; iwire < slc.lastWire[plane] - 3; ++iwire) {
      // skip bad wires or no hits on the wire
      if (slc.wireHitRange[plane][iwire].first > slc.slHits.size()) continue;
      if (slc.wireHitRange[plane][iwire].second > slc.slHits.size()) continue;
      unsigned int jwire = iwire + 1;
      if (slc.wireHitRange[plane][jwire].first > slc.slHits.size()) continue;
      if (slc.wireHitRange[plane][jwire].second > slc.slHits.size()) continue;
      unsigned int ifirsthit = slc.wireHitRange[plane][iwire].first;
      unsigned int ilasthit = slc.wireHitRange[plane][iwire].second;
      unsigned int jfirsthit = slc.wireHitRange[plane][jwire].first;
      unsigned int jlasthit = slc.wireHitRange[plane][jwire].second;
      for (unsigned int iht = ifirsthit; iht <= ilasthit; ++iht) {
        if(iht >= slc.slHits.size()) break;
        auto& islHit = slc.slHits[iht];
        if (islHit.InTraj != 0) continue;
        std::vector<unsigned int> iHits;
        if(tcc.useAlg[kNewCuts]) {
          iHits = FindJTHits(slc, iht);
        } else {
          GetHitMultiplet(slc, iht, iHits, true);
        }
        prt =
          (tcc.modes[kModeDebug] && std::find(iHits.begin(), iHits.end(), debug.Hit) != iHits.end());
        if (prt) mf::LogVerbatim("TC") << "FJT: debug iht multiplet size " << iHits.size();
        if (iHits.empty()) continue;
        for (unsigned int jht = jfirsthit; jht <= jlasthit; ++jht) {
          if(jht >= slc.slHits.size()) break;
          auto& jslHit = slc.slHits[jht];
          if (jslHit.InTraj != 0) continue;
          if (prt && HitSep2(slc, iht, jht) < 100)
            mf::LogVerbatim("TC") << " use " << PrintHit(jslHit) << " hitSep2 "
                                  << HitSep2(slc, iht, jht);
          if (HitSep2(slc, iht, jht) > tcc.JTMaxHitSep2) continue;
          std::vector<unsigned int> jHits;
          if(tcc.useAlg[kNewCuts]) {
            jHits = FindJTHits(slc, jht);
          } else {
            GetHitMultiplet(slc, jht, jHits, true);
          }
          if (jHits.empty()) continue;
          // check for hit overlap consistency
          if(tcc.useAlg[kNewCuts]) {
            if (!JTHitsOK(slc, iHits, jHits)) continue;
          } else {
            if (!TrajHitsOK(slc, iHits, jHits)) continue;
          }
          tHits.clear();
          // add the available hits and flag them
          for (auto iht : iHits)
            if (slc.slHits[iht].InTraj == 0) tHits.push_back(iht);
          for (auto jht : jHits)
            if (slc.slHits[jht].InTraj == 0) tHits.push_back(jht);
          for (auto tht : tHits)
            slc.slHits[tht].InTraj = -4;
          unsigned int loWire;
          if (iwire != 0) { loWire = iwire - 1; }
          else {
            loWire = 0;
          }
          unsigned int hiWire = jwire + 1;
          if (hiWire > slc.nWires[plane]) break;
          unsigned short nit = 0;
          while (nit < 100) {
            bool hitsAdded = false;
            for (unsigned int kwire = loWire; kwire <= hiWire; ++kwire) {
              if (slc.wireHitRange[plane][kwire].first == UINT_MAX) continue;
              if (slc.wireHitRange[plane][kwire].second == UINT_MAX) continue;
              unsigned int kfirsthit = slc.wireHitRange[plane][kwire].first;
              unsigned int klasthit = slc.wireHitRange[plane][kwire].second;
              for (unsigned int kht = kfirsthit; kht <= klasthit; ++kht) {
                if(kht >= slc.slHits.size()) continue;
                if (slc.slHits[kht].InTraj != 0) continue;
                // this shouldn't be needed but do it anyway
                if (std::find(tHits.begin(), tHits.end(), kht) != tHits.end()) continue;
                // re-purpose jHits and check for consistency
                if(tcc.useAlg[kNewCuts]) {
                  jHits = FindJTHits(slc, kht);
                } else {
                  GetHitMultiplet(slc, kht, jHits, true);
                }
                if(tcc.useAlg[kNewCuts]) {
                  if (!JTHitsOK(slc, tHits, jHits)) continue;
                } else {
                  if (!TrajHitsOK(slc, tHits, jHits)) continue;
                }
                // add them all and update the wire range
                for (auto jht : jHits) {
                  if (slc.slHits[jht].InTraj != 0) continue;
                  tHits.push_back(jht);
                  slc.slHits[jht].InTraj = -4;
                  if (kwire > hiWire) hiWire = kwire;
                  if (kwire < loWire) loWire = kwire;
                  hitsAdded = true;
                } // jht
                // allow continuing if a wire has hits that are already used
                if(tcc.useAlg[kNewCuts] && !jHits.empty()) hitsAdded = true;
              }   // kht
            }     // kwire
            if (!hitsAdded) break;
            ++nit;
            ++hiWire;
            if (hiWire >= slc.nWires[plane]) break;
          } // nit < 100
          // clear InTraj
          for (auto iht : tHits)
            slc.slHits[iht].InTraj = 0;
          if (tHits.size() < 3) continue;
          if (prt) {
            mf::LogVerbatim myprt("TC");
            myprt << "FJT: tHits";
            for (auto tht : tHits)
              myprt << " " << PrintHit(slc.slHits[tht]);
            myprt << "\n";
          } // prt
          // See if this is a ghost trajectory
          if (IsGhost(slc, tHits)) break;
          if (!MakeJunkTraj(slc, tHits)) {
            if (prt) mf::LogVerbatim("TC") << "FJT: MakeJunkTraj failed";
            break;
          }
          if (slc.slHits[jht].InTraj > 0) break;
        } // jht
      }   // iht
    }     // iwire
  }       // FindJunkTraj

  ////////////////////////////////////////////////
  std::vector<unsigned int>
  TrajClusterAlg::FindJTHits(const TCSlice& slc, unsigned int iht)
  {
    // a helper function for FindJunkTraj
    std::vector<unsigned int> hitList;
    if(iht >= slc.slHits.size()) return hitList;
    auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
    unsigned short plane = hit.WireID().Plane;
    int wire = hit.WireID().Wire;
    std::array<int, 2> wireWindow = {wire, wire};
    Point2_t timeWindow;
    timeWindow[0] = hit.PeakTime() * tcc.unitsPerTick;
    timeWindow[1] = timeWindow[0] + sqrt(tcc.JTMaxHitSep2);
    bool hitsNear = false;
    auto closeHits = FindCloseHits(slc, wireWindow, timeWindow, plane, kUnusedHits, true, hitsNear);
    if(closeHits.empty()) return closeHits;
    // include all hits that are in a multiplet even though that may
    // lie outside the timeWindow. We only need to check the first
    // and last hits since they are time-ordered
    for(unsigned short chk = 0; chk < 2; ++chk) {
      unsigned int iht = closeHits[0];
      if(chk > 0) iht = closeHits.back();
      if(iht >= slc.slHits.size()) continue;
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      if(hit.Multiplicity() == 1) continue;
      // index of the first hit, using the assumption that all hits in the
      // multiplet are correctly included in this slice
      if(hit.LocalIndex() > (int)iht) continue;
      unsigned int fht = iht - hit.LocalIndex();
      // ensure that this assumption is correct
      if((*evt.allHits)[slc.slHits[fht].allHitsIndex].LocalIndex() != 0) return closeHits;
      for(short li = 0; li < hit.Multiplicity(); ++li) {
        unsigned int jht = fht + li;
        if(std::find(closeHits.begin(), closeHits.end(), jht) == closeHits.end()) closeHits.push_back(jht);
      } // li
      if(closeHits.size() == 1) break;
    } // chk
    return closeHits;
  } // FindJTHits

  //////////////////////////////////////////
  void
  TrajClusterAlg::MergeTPHits(std::vector<unsigned int>& tpHits,
                              std::vector<recob::Hit>& newHitCol,
                              std::vector<unsigned int>& newHitAssns) const
  {
    // merge the hits indexed by tpHits into one or more hits with the requirement that the hits
    // are on different wires

    if (tpHits.empty()) return;

    // no merge required. Just put a close copy of the single hit in the output hit collection
    if (tpHits.size() == 1) {
      if (tpHits[0] > (*evt.allHits).size() - 1) {
        std::cout << "MergeTPHits Bad input hit index " << tpHits[0] << " allHits size "
                  << (*evt.allHits).size() << "\n";
        return;
      }
      newHitCol.push_back(MergeTPHitsOnWire(tpHits));
      newHitAssns[tpHits[0]] = newHitCol.size() - 1;
      return;
    } // tpHits.size() == 1

    // split the hit list into sub-lists of hits on a single wire
    std::vector<unsigned int> wires;
    std::vector<std::vector<unsigned int>> wireHits;
    auto& firstHit = (*evt.allHits)[tpHits[0]];
    wires.push_back(firstHit.WireID().Wire);
    std::vector<unsigned int> tmp(1, tpHits[0]);
    wireHits.push_back(tmp);
    for (unsigned short ii = 1; ii < tpHits.size(); ++ii) {
      auto& hit = (*evt.allHits)[tpHits[ii]];
      unsigned int wire = hit.WireID().Wire;
      unsigned short indx = 0;
      for (indx = 0; indx < wires.size(); ++indx)
        if (wires[indx] == wire) break;
      if (indx == wires.size()) {
        wires.push_back(wire);
        wireHits.resize(wireHits.size() + 1);
      }
      wireHits[indx].push_back(tpHits[ii]);
    } // ii

    // now merge hits in each sub-list.
    for (unsigned short indx = 0; indx < wireHits.size(); ++indx) {
      auto& hitsOnWire = wireHits[indx];
      newHitCol.push_back(MergeTPHitsOnWire(hitsOnWire));
      for (unsigned short ii = 0; ii < hitsOnWire.size(); ++ii) {
        newHitAssns[hitsOnWire[ii]] = newHitCol.size() - 1;
      }
    } // hitsOnWire

    return;

  } // MergeTPHits

  //////////////////////////////////////////
  recob::Hit
  TrajClusterAlg::MergeTPHitsOnWire(std::vector<unsigned int>& tpHits) const
  {
    // merge the hits indexed by tpHits into one hit

    if (tpHits.empty()) return recob::Hit();

    // no merge required. Just return a slightly modified copy of the single hit
    if (tpHits.size() == 1) {
      if (tpHits[0] > (*evt.allHits).size() - 1) {
        std::cout << "MergeTPHits Bad input hit index " << tpHits[0] << " allHits size "
                  << (*evt.allHits).size() << "\n";
        return recob::Hit();
      }
      auto& oldHit = (*evt.allHits)[tpHits[0]];
      raw::TDCtick_t startTick = oldHit.PeakTime() - 3 * oldHit.RMS();
      raw::TDCtick_t endTick = oldHit.PeakTime() + 3 * oldHit.RMS();

      return recob::Hit(oldHit.Channel(),
                        startTick,
                        endTick,
                        oldHit.PeakTime(),
                        oldHit.SigmaPeakTime(),
                        oldHit.RMS(),
                        oldHit.PeakAmplitude(),
                        oldHit.SigmaPeakAmplitude(),
                        oldHit.SummedADC(),
                        oldHit.Integral(),
                        oldHit.SigmaIntegral(),
                        1,
                        0, // Multiplicity, LocalIndex
                        1,
                        0, // GoodnessOfFit, DOF
                        oldHit.View(),
                        oldHit.SignalType(),
                        oldHit.WireID());
    } // tpHits.size() == 1

    double integral = 0;
    double sIntegral = 0;
    double peakTime = 0;
    double sPeakTime = 0;
    double peakAmp = 0;
    double sPeakAmp = 0;
    float sumADC = 0;
    raw::TDCtick_t startTick = INT_MAX;
    raw::TDCtick_t endTick = 0;
    for (auto allHitsIndex : tpHits) {
      if (allHitsIndex > (*evt.allHits).size() - 1) return recob::Hit();
      auto& hit = (*evt.allHits)[allHitsIndex];
      if (hit.StartTick() < startTick) startTick = hit.StartTick();
      if (hit.EndTick() > endTick) endTick = hit.EndTick();
      double intgrl = hit.Integral();
      sPeakTime += intgrl * hit.SigmaPeakTime();
      sPeakAmp += intgrl * hit.SigmaPeakAmplitude();
      sumADC += hit.SummedADC();
      integral += intgrl;
      sIntegral += intgrl * hit.SigmaIntegral();
      // Get the charge normalization from an input hit
    } // tpHit
    if (integral <= 0) {
      std::cout << "MergeTPHits found bad hit integral " << integral << "\n";
      return recob::Hit();
    }

    // Create a signal shape vector to find the rms and peak tick
    std::vector<double> shape(endTick - startTick + 1, 0.);
    for (auto allHitsIndex : tpHits) {
      auto& hit = (*evt.allHits)[allHitsIndex];
      double peakTick = hit.PeakTime();
      double rms = hit.RMS();
      double peakAmp = hit.PeakAmplitude();
      for (raw::TDCtick_t tick = startTick; tick <= endTick; ++tick) {
        double arg = ((double)tick - peakTick) / rms;
        unsigned short indx = tick - startTick;
        shape[indx] += peakAmp * exp(-0.5 * arg * arg);
      } // tick
    }   // allHitsIndex

    // find the peak tick
    double sigsum = 0;
    double sigsumt = 0;
    for (raw::TDCtick_t tick = startTick; tick <= endTick; ++tick) {
      unsigned short indx = tick - startTick;
      sigsum += shape[indx];
      sigsumt += shape[indx] * tick;
    } // tick

    peakTime = sigsumt / sigsum;
    // Use the sigma peak time calculated in the first loop
    sPeakTime /= integral;

    // find the RMS
    sigsumt = 0.;
    for (raw::TDCtick_t tick = startTick; tick <= endTick; ++tick) {
      double dTick = tick - peakTime;
      unsigned short indx = tick - startTick;
      sigsumt += shape[indx] * dTick * dTick;
    }
    double rms = std::sqrt(sigsumt / sigsum);
    // get a reference to the first hit to get the charge normalization, channel, view, etc
    auto& firstHit = (*evt.allHits)[tpHits[0]];
    double chgNorm = 2.507 * firstHit.PeakAmplitude() * firstHit.RMS() / firstHit.Integral();
    // find the amplitude from the integrated charge and the RMS
    peakAmp = (float)(integral * chgNorm / (2.507 * rms));
    // Use the sigma integral calculated in the first loop
    sPeakAmp /= integral;
    // reset the start and end tick
    startTick = peakTime - 3 * rms;
    endTick = peakTime + 3 * rms;

    // construct the hit
    return recob::Hit(firstHit.Channel(),
                      startTick,
                      endTick,
                      peakTime,
                      sPeakTime,
                      rms,
                      peakAmp,
                      sPeakAmp,
                      sumADC,
                      integral,
                      sIntegral,
                      1,
                      0, // Multiplicity, LocalIndex
                      1,
                      0, // GoodnessOfFit, DOF
                      firstHit.View(),
                      firstHit.SignalType(),
                      firstHit.WireID());

  } // MergeTPHits

  /////////////////////////////////////////
  bool
  TrajClusterAlg::CreateSlice(detinfo::DetectorClocksData const& clockData,
                              detinfo::DetectorPropertiesData const& detProp,
                              std::vector<unsigned int>& hitsInSlice,
                              int sliceID)
  {
    // Defines a TCSlice struct and pushes the slice onto slices.
    // Sets the isValid flag true if successful.
    if ((*evt.allHits).empty()) return false;
    if (hitsInSlice.size() < 2) return false;

    TCSlice slc;
    slc.ID = sliceID;
    slc.slHits.resize(hitsInSlice.size());
    bool first = true;
    unsigned int cstat = 0;
    unsigned int tpc = UINT_MAX;
    unsigned int cnt = 0;
    std::vector<unsigned int> nHitsInPln;
    for (auto iht : hitsInSlice) {
      if (iht >= (*evt.allHits).size()) return false;
      auto& hit = (*evt.allHits)[iht];
      if (first) {
        cstat = hit.WireID().Cryostat;
        tpc = hit.WireID().TPC;
        slc.TPCID = geo::TPCID(cstat, tpc);
        nHitsInPln.resize(tcc.geom->Nplanes(slc.TPCID));
        first = false;
      }
      if (hit.WireID().Cryostat != cstat || hit.WireID().TPC != tpc) return false;
      slc.slHits[cnt].allHitsIndex = iht;
      ++nHitsInPln[hit.WireID().Plane];
      ++cnt;
    } // iht
    // require at least two hits in each plane
    for (auto hip : nHitsInPln)
      if (hip < 2) return false;
    // Define the TCEvent wire hit range vector for this new TPC for ALL hits
    FillWireHitRange(slc.TPCID);
    // next define the Slice wire hit range vectors, UnitsPerTick, etc for this
    // slice
    if (!FillWireHitRange(clockData, detProp, slc)) return false;
    slc.isValid = true;
    slices.push_back(slc);
    if (tcc.modes[kModeDebug] && debug.Slice >= 0 && !tcc.dbgSlc) {
      tcc.dbgSlc = ((int)(slices.size() - 1) == debug.Slice);
      if (tcc.dbgSlc) std::cout << "Enabled debugging in sub-slice " << slices.size() - 1 << "\n";
      if (tcc.modes[kModeDebug] && (unsigned int)debug.Cryostat == cstat &&
          (unsigned int)debug.TPC == tpc && debug.Plane >= 0) {
        debug.CTP = EncodeCTP(
          (unsigned int)debug.Cryostat, (unsigned int)debug.TPC, (unsigned int)debug.Plane);
      }
    }
    // do a sanity check
    for (unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      unsigned int ahi = slc.slHits[iht].allHitsIndex;
      if (ahi >= (*evt.allHits).size()) {
        std::cout << "CreateSlice: Bad allHitsIndex " << ahi << " " << (*evt.allHits).size()
                  << "\n";
        return false;
      }
    } // iht
    return true;
  } // CreateSlice

  /////////////////////////////////////////
  void
  TrajClusterAlg::FinishEvent(detinfo::DetectorPropertiesData const& detProp)
  {
    // final steps that involve correlations between slices
    // Stitch PFParticles between TPCs

    // define the PFP TjUIDs vector before calling StitchPFPs
    for (auto& slc : slices) {
      if (!slc.isValid) continue;
      MakePFPTjs(slc);
      for (auto& pfp : slc.pfps) {
        if (pfp.ID <= 0) continue;
        pfp.TjUIDs.resize(pfp.TjIDs.size());
        for (unsigned short ii = 0; ii < pfp.TjIDs.size(); ++ii) {
          // do a sanity check while we are here
          if (pfp.TjIDs[ii] <= 0 || pfp.TjIDs[ii] > (int)slc.tjs.size()) {
            std::cout << "FinishEvent found an invalid T" << pfp.TjIDs[ii] << " in P" << pfp.UID
                      << "\n";
            slc.isValid = false;
            continue;
          } // sanity check
          auto& tj = slc.tjs[pfp.TjIDs[ii] - 1];
          pfp.TjUIDs[ii] = tj.UID;
        } // ii
      }   // pfp
    }     // slc

    StitchPFPs();
    // Ensure that all PFParticles have a start vertex
    for (auto& slc : slices) PFPVertexCheck(slc);
    if(tcc.dbgSummary) PrintAll(detProp, "Fin");
  } // FinishEvent

  /////////////////////////////////////////
  void TrajClusterAlg::MakeSpacePointsFromPFP(const tca::PFPStruct& pfp,
       const std::vector<unsigned int>& newHitIndex, std::vector<recob::SpacePoint>& spts, 
       std::vector<unsigned int>& sptsHit)
  {
    // Converts a PFPStruct into a set of recob::SpacePoints + hit assn

    spts.clear();
    sptsHit.clear();
    if(pfp.TP3Ds.empty()) return;
    if(pfp.ID <= 0) return;
    auto slcIndx = GetSliceIndex("P", pfp.UID);
    if(slcIndx.first == USHRT_MAX) return;
    auto& slc = slices[slcIndx.first];

    int id = 0;
    for(unsigned int pt = 0; pt < pfp.TP3Ds.size(); ++pt) {
      auto& tp3d = pfp.TP3Ds[pt];
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      // index of the hit in the new hit collection
      unsigned int nhi = UINT_MAX;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
        if(newHitIndex[ahi] == UINT_MAX) continue;
        nhi = newHitIndex[ahi];
        break;
      } // ii
      if(nhi == UINT_MAX) continue;
      Double32_t pos[3] = {tp3d.Pos[0], tp3d.Pos[1], tp3d.Pos[2]};
      Double32_t err[3] = {0., 0, 0.};
      ++id;
      spts.emplace_back(pos, err, 0., id);
      sptsHit.push_back(nhi);
    } // pt

  } // MakeSpacePointsFromPFP

  /////////////////////////////////////////
  void TrajClusterAlg::MakeTrackFromPFP(const tca::PFPStruct& pfp,
       const std::vector<unsigned int>& newHitIndex, recob::Track& trk, 
       std::vector<unsigned int>& trkHits)
  {
    trkHits.clear();
    // Converts a PFPStruct into a recob::Track and TrackHitMeta
    if(pfp.TP3Ds.size() < 4) return;
    if(pfp.ID <= 0) return;
    auto slcIndx = GetSliceIndex("P", pfp.UID);
    if(slcIndx.first == USHRT_MAX) return;
    auto& slc = slices[slcIndx.first];

    std::vector<recob::Track::Point_t> positions;
    std::vector<recob::Track::Vector_t> directions;
    std::vector<recob::TrajectoryPointFlags> tpFlags;
    // index of the hit in the new hit collection

    for(unsigned int pt = 0; pt < pfp.TP3Ds.size(); ++pt) {
      auto& tp3d = pfp.TP3Ds[pt];
      if(tp3d.TPIndex >= (*evt.allHits).size()) continue;
      // construct the flag
      unsigned int nhi = UINT_MAX;
      auto mask = recob::TrajectoryPointFlags::makeMask();

      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
        if(newHitIndex[ahi] == UINT_MAX) continue;
        nhi = newHitIndex[ahi];
        break;
      } // ii
      if(nhi == UINT_MAX) continue;
      // Map the Track flag bits from the TP Environment bitset
      for(unsigned short ib = 0; ib < 8; ++ib) {
        if(!tp.Environment[ib]) continue;
        if(ib == kEnvOverlap) mask.set(recob::TrajectoryPointFlagTraits::Shared);
        if(ib == kEnvNotGoodWire) mask.set(recob::TrajectoryPointFlagTraits::DetectorIssue);
        // Not all hits in the multiplet were used in the TP
        if(ib == kEnvUnusedHits) mask.set(recob::TrajectoryPointFlagTraits::Suspicious);
      } // ib

      recob::Track::Point_t pos = {tp3d.Pos[0], tp3d.Pos[1], tp3d.Pos[2]};
      recob::Track::Vector_t dir = {tp3d.Dir[0], tp3d.Dir[1], tp3d.Dir[2]};
      positions.push_back(pos);
      directions.push_back(dir);
      trkHits.push_back(nhi);
      tpFlags.push_back(recob::TrajectoryPointFlags(nhi, mask));
    } // pt
    if(trkHits.size() < 4) {
      trkHits.clear();
      return;
    }
    // All the SectionFit fits are independent so just do a simple sum
    int ndof = 4 * pfp.SectionFits.size();
    float chi = 0.;
    for(auto& sf : pfp.SectionFits) chi += sf.ChiDOF * 4;

    // construct the track, which has a trajectory with "momentum", a stab at a PDGCode,
    // a chisq and not-defined covariance matrix
    recob::tracking::SMatrixSym55 cov;
    trk = recob::Track(recob::TrackTrajectory(std::move(positions),
                      std::move(directions),std::move(tpFlags), true),
                      pfp.PDGCode, chi, ndof,
                      cov, cov, pfp.UID);
  } // MakeTrackFromPFP


} // namespace cluster
