////////////////////////////////////////////////////////////////////////
// Class:       ClusterTrackAna
// Plugin Type: analyzer (art v3_03_01)
// File:        ClusterTrackAna_module.cc
//
// Calculates the performance of 2D cluster and 3D track reconstruction modules.
// The metrics Efficiency, Purity, and Efficiency * Purity (EP) are calculated using
// MC truth matched hits in each TPC and Plane. Note that the metrics are calculated in
// each TPC and plane for clusters or tracks that span multiple TPCs. The metrics are summed and
// reported at the end of the job for each selected true particle type. The metrics are also
// summed over ALL selected true particle types. Note that the term "cluster" refers to
// a subset of hits in a recob::Track or recob::Cluster that reside in one TPC and one plane
// and may be matched to a simb::MCParticle
//
// Generated at Wed Jan  8 10:33:20 2020 by Bruce Baller using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cluster {
  class ClusterTrackAna;
}

class cluster::ClusterTrackAna : public art::EDAnalyzer {
public:
  explicit ClusterTrackAna(fhicl::ParameterSet const& p);

  ClusterTrackAna(ClusterTrackAna const&) = delete;
  ClusterTrackAna(ClusterTrackAna&&) = delete;
  ClusterTrackAna& operator=(ClusterTrackAna const&) = delete;
  ClusterTrackAna& operator=(ClusterTrackAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  void endJob() override;

  art::InputTag fHitModuleLabel;
  art::InputTag fClusterModuleLabel;
  art::InputTag fTrackModuleLabel;
  simb::Origin_t fTruthOrigin;
  std::vector<int> fSkipPDGCodes;
  short fPrintLevel;
  unsigned int fInTPC;
  float fBadEP;

  // Count of the number of events considered
  unsigned int fEventCnt{0};
  // Count of the number of poorly reconstructed clusters
  unsigned int fNBadEP{0};
  // count of EP entries for electrons(0), muons(1), pions(2), kaons(3), protons(4)
  std::array<float, 5> Cnts{{0}};
  std::array<float, 5> EPSums{{0}};
  // same for Efficiency
  std::array<float, 5> ESums{{0}};
  // and for Purity
  std::array<float, 5> PSums{{0}};

  art::Handle<std::vector<recob::Hit>> fHitHandle;
  std::vector<unsigned int> fHitMCPIndex;

  bool fCompareProductIDs{
    true}; ///< compare Hit and Cluster-> Hit art product IDs on the first event
  bool fFirstPrint{true};
  unsigned int fCurrentRun{0};

  void FirstLastHitInPlane(unsigned int tpc,
                           unsigned int plane,
                           unsigned int mcpi,
                           unsigned int& firstHitIndex,
                           unsigned int& lastHitIndex);
  std::string HitLocation(
    unsigned int iht); ///< packs hit WireID and PeakTime into a compact format
};

cluster::ClusterTrackAna::ClusterTrackAna(fhicl::ParameterSet const& pset) : EDAnalyzer{pset}
{
  fHitModuleLabel = pset.get<art::InputTag>("HitModuleLabel");
  bool validModuleLabel = false;
  fClusterModuleLabel = "NA";
  fTrackModuleLabel = "NA";
  fClusterModuleLabel = pset.get<art::InputTag>("ClusterModuleLabel", "NA");
  if (fClusterModuleLabel != "NA") validModuleLabel = true;
  fTrackModuleLabel = pset.get<art::InputTag>("TrackModuleLabel", "NA");
  if (validModuleLabel && fTrackModuleLabel != "NA")
    throw cet::exception("ClusterTrackAna")
      << "You must specify either a ClusterModuleLabel OR a TrackModuleLabel\n";
  if (!validModuleLabel && fTrackModuleLabel != "NA") validModuleLabel = true;
  // origin = 0 (anything), 1(nu), 2(cosmics), 3(SN nu), 4(SingleParticle)
  int tmp = pset.get<int>("TruthOrigin", 0);
  fTruthOrigin = (simb::Origin_t)tmp;
  fPrintLevel = pset.get<short>("PrintLevel", 0);
  if (pset.has_key("SkipPDGCodes")) fSkipPDGCodes = pset.get<std::vector<int>>("SkipPDGCodes");
  fBadEP = pset.get<float>("BadEP", 0.);
  fInTPC = UINT_MAX;
  int intpc = pset.get<int>("InTPC", -1);
  if (intpc >= 0) fInTPC = intpc;
  // do some initialization
  Cnts.fill(0.);
  EPSums.fill(0.);
  ESums.fill(0.);
  PSums.fill(0.);
} // ClusterTrackAna constructor

////////////////////////////////////////////////
void cluster::ClusterTrackAna::analyze(art::Event const& evt)
{
  // Match hits to MCParticles, then consider reconstructed hits in each TPC and plane
  // to calculate Efficiency, Purity and Efficiency * Purity (aka EP).

  ++fEventCnt;
  auto const* geom = lar::providerFrom<geo::Geometry>();
  //  auto hitsHandle = art::Handle<std::vector<recob::Hit>>();
  if (!evt.getByLabel(fHitModuleLabel, fHitHandle))
    throw cet::exception("ClusterTrackAna")
      << "Failed to get a handle to hit collection '" << fHitModuleLabel.label() << "'\n";
  // get a reference to the MCParticles
  auto mcpHandle = art::Handle<std::vector<simb::MCParticle>>();
  if (!evt.getByLabel("largeant", mcpHandle))
    throw cet::exception("ClusterTrackAna")
      << "Failed to get a handle to MCParticles using largeant\n";

  // decide whether to consider cluster -> hit -> MCParticle for any MCParticle origin or for
  // a specific user-specified origin
  bool anySource = (fTruthOrigin == simb::kUnknown);

  if (fPrintLevel > 0 && evt.run() != fCurrentRun) {
    mf::LogVerbatim("ClusterTrackAna") << "Run: " << evt.run();
    fCurrentRun = evt.run();
  }

  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
  sim::ParticleList const& plist = pi_serv->ParticleList();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  // make a list of Hit -> MCParticle assns in all TPCs. The first step is
  // to make a list of Geant TrackIDs whose origin was requested by the user.
  std::vector<int> trkIDs;
  // and a vector of MCParticle indices into the mcpHandle vector indexed by trkIDs
  std::vector<unsigned int> MCPIs;
  for (sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
    const simb::MCParticle* part = (*ipart).second;
    int trackID = part->TrackId();
    art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
    if (!anySource && theTruth->Origin() != fTruthOrigin) continue;
    int pdg = abs(part->PdgCode());
    bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
    if (!isCharged) continue;
    if (std::find(fSkipPDGCodes.begin(), fSkipPDGCodes.end(), pdg) != fSkipPDGCodes.end()) continue;
    float TMeV = 1000 * (part->E() - part->Mass());
    // ignore low energy particles
    if (TMeV < 1) continue;
    // find the MCParticle index and stash it in trkMCP
    unsigned int mcpi = UINT_MAX;
    for (mcpi = 0; mcpi < (*mcpHandle).size(); ++mcpi)
      if ((*mcpHandle)[mcpi].TrackId() == trackID) break;
    if (mcpi == UINT_MAX) {
      std::cout << "Failed to find a MCParticle from ParticleList\n";
      return;
    }
    trkIDs.push_back(trackID);
    MCPIs.push_back(mcpi);
  } // ipart
  if (trkIDs.empty()) return;

  // next construct a companion vector of MCParticle indices indexed to the full hit collection
  fHitMCPIndex.resize((*fHitHandle).size(), UINT_MAX);
  // Make a list of the <first, last> MC-matched hit in each TPC. This will be used
  // to only iterate through the range of hits that are interesting
  std::vector<std::pair<unsigned int, unsigned int>> hitRange(geom->NTPC() + 1);
  size_t noMatch = 0;
  size_t nMatch = 0;
  size_t nInTPC = 0;
  for (auto& hr : hitRange)
    hr = std::make_pair(UINT_MAX, UINT_MAX);
  for (size_t iht = 0; iht < (*fHitHandle).size(); ++iht) {
    auto& hit = (*fHitHandle)[iht];
    // only consider hits in a single TPC?
    unsigned int tpc = hit.WireID().TPC;
    if (fInTPC != UINT_MAX && tpc != fInTPC) continue;
    ++nInTPC;
    auto tides = bt_serv->HitToTrackIDEs(clockData, hit);
    if (tides.empty()) {
      ++noMatch;
      continue;
    }
    for (auto& tide : tides) {
      // declare a match to a MCParticle if > 50% of energy is from it
      if (tide.energyFrac < 0.5) continue;
      int trackID = tide.trackID;
      // find the MCParticle index and define the Hit -> MCParticle assn
      for (size_t indx = 0; indx < trkIDs.size(); ++indx) {
        if (trkIDs[indx] == trackID) {
          fHitMCPIndex[iht] = MCPIs[indx];
          break;
        }
      } // indx
      break;
    } // tide
    if (fHitMCPIndex[iht] == UINT_MAX) continue;
    ++nMatch;
    // populate hitRange
    if (hitRange[tpc].first == UINT_MAX) hitRange[tpc].first = iht;
    hitRange[tpc].second = iht;
  } // iht

  if (nMatch == 0) {
    fHitMCPIndex.resize(0);
    return;
  }

  if (fPrintLevel > 3) {
    // Print the gory details
    mf::LogVerbatim myprt("ClusterTrackAna");
    myprt << "Checking " << trkIDs.size() << " selected MCParticles with the requested TruthOrigin";
    if (fInTPC == UINT_MAX) { myprt << " in all TPCs\n"; }
    else {
      myprt << " in TPC " << fInTPC;
    }
    myprt << " in Run " << evt.run() << " Event " << evt.event();
    myprt << "\n";
    myprt << "Found " << nMatch << " MC-matched hits with the requested TruthOrigin";
    myprt << " out of " << nInTPC << " total hits.\n";
    myprt << "Found " << noMatch << " hits not matched to ANY MCParticle.\n";
    // count the number of hits matched to each MCParticle in the list
    //                     mcpi          count
    std::vector<std::pair<unsigned int, unsigned int>> mcpCnt;
    for (auto mcpi : fHitMCPIndex) {
      if (mcpi == UINT_MAX) continue;
      unsigned short mIndx = 0;
      for (mIndx = 0; mIndx < mcpCnt.size(); ++mIndx)
        if (mcpCnt[mIndx].first == mcpi) break;
      if (mIndx == mcpCnt.size()) mcpCnt.push_back(std::make_pair(mcpi, 0));
      // increment the count of MC-matched hits
      ++mcpCnt[mIndx].second;
    } // mcpi
    myprt << " MCPI TrackID  PDGCode  T(MeV)   nHits    Process";
    for (auto mcpcnt : mcpCnt) {
      myprt << "\n";
      unsigned int mcpi = mcpcnt.first;
      auto& mcp = (*mcpHandle)[mcpi];
      myprt << std::setw(5) << mcpi;
      myprt << std::setw(8) << mcp.TrackId();
      myprt << std::setw(8) << abs(mcp.PdgCode());
      int TMeV = 1000 * (mcp.E() - mcp.Mass());
      myprt << std::setw(9) << TMeV;
      myprt << std::setw(8) << mcpcnt.second;
      myprt << "    " << mcp.Process();
    } // mcpcnt
  }   // fPrintLevel > 3

  // fill a vector of hits indices from all clusters or all tracks
  std::vector<std::vector<unsigned int>> recoHits;
  // and a companion vector of indices into the cluster or track collection
  std::vector<unsigned int> recoIndex;
  // handles to these collections to enable printing the cluster or track ID. The
  // handle will be valid for either clusters or tracks
  art::Handle<std::vector<recob::Cluster>> inputClusters;
  art::Handle<std::vector<recob::Track>> inputTracks;
  if (fClusterModuleLabel.label() != "NA") {
    // get MC-matched hits from clusters
    if (!evt.getByLabel(fClusterModuleLabel, inputClusters))
      throw cet::exception("ClusterTrackAna")
        << "Failed to get a handle to the cluster collection '" << fClusterModuleLabel.label()
        << "'\n";
    art::FindManyP<recob::Hit> hitsFromCls(inputClusters, evt, fClusterModuleLabel);
    if (!hitsFromCls.isValid())
      throw cet::exception("ClusterTrackAna") << "Failed to get a handle to Cluster -> Hit assns\n";
    for (unsigned int icl = 0; icl < (*inputClusters).size(); ++icl) {
      std::vector<art::Ptr<recob::Hit>> cluhits = hitsFromCls.at(icl);
      if (cluhits.empty()) continue;
      if (fCompareProductIDs) {
        if (cluhits[0].id() != fHitHandle.id())
          throw cet::exception("ClusterTrackAna")
            << "Hits associated with ClusterModuleLabel are in a different collection than "
               "HitModuleLabel.\n";
        fCompareProductIDs = false;
      } // fCompareProductIDs
      // only load MC-matched hits. Hits that are not MC-matched were either matched to a
      // MCParticle that was not selected or were mis-reconstructed by the hit finder. Neither
      // of these conditions warrant penalizing the reconstruction module performance
      std::vector<unsigned int> hitsIndex;
      for (auto& cluhit : cluhits) {
        if (fHitMCPIndex[cluhit.key()] != UINT_MAX) hitsIndex.push_back(cluhit.key());
      }
      if (hitsIndex.empty()) continue;
      recoIndex.push_back(icl);
      recoHits.push_back(hitsIndex);
    } // icl
  }
  else {
    // get MC-matched hits from tracks
    if (!evt.getByLabel(fTrackModuleLabel, inputTracks))
      throw cet::exception("ClusterTrackAna")
        << "Failed to get a handle to the track collection '" << fTrackModuleLabel.label() << "'\n";
    art::FindManyP<recob::Hit> hitsFromTrk(inputTracks, evt, fTrackModuleLabel);
    if (!hitsFromTrk.isValid())
      throw cet::exception("ClusterTrackAna") << "Failed to get a handle to Track -> Hit assns\n";
    for (unsigned int itk = 0; itk < (*inputTracks).size(); ++itk) {
      std::vector<art::Ptr<recob::Hit>> trkhits = hitsFromTrk.at(itk);
      if (trkhits.empty()) continue;
      if (fCompareProductIDs) {
        if (trkhits[0].id() != fHitHandle.id())
          throw cet::exception("ClusterTrackAna")
            << "Hits associated with TrackModuleLabel are in a different collection than "
               "HitModuleLabel.\n";
        fCompareProductIDs = false;
      } // fCompareProductIDs
      std::vector<unsigned int> hitsIndex;
      for (auto& trkhit : trkhits) {
        if (fHitMCPIndex[trkhit.key()] != UINT_MAX) hitsIndex.push_back(trkhit.key());
      }
      if (hitsIndex.empty()) continue;
      recoIndex.push_back(itk);
      recoHits.push_back(hitsIndex);
    } // itk
  }   // get hits from tracks

  for (const auto& tpcid : geom->Iterate<geo::TPCID>()) {
    unsigned int tpc = tpcid.TPC;
    if (hitRange[tpc].first == UINT_MAX) continue;
    // iterate over planes
    for (unsigned short plane = 0; plane < geom->Nplanes(); ++plane) {
      unsigned int tpcMatHitCnt = 0;
      unsigned int tpcTotHitCnt = 0;
      // create a list of (MCParticle index, matched hit count> pairs
      //  mcParticle  plane
      std::vector<std::pair<unsigned int, float>> mcpCnt;
      // count MCParticle, Cluster/Track hit counts - size matched to mcpCnt
      std::vector<std::vector<std::pair<unsigned int, float>>> mcpClsCnt;
      for (unsigned int iht = hitRange[tpc].first; iht <= hitRange[tpc].second; ++iht) {
        auto& hit = (*fHitHandle)[iht];
        if (hit.WireID().TPC != tpc) continue;
        if (hit.WireID().Plane != plane) continue;
        ++tpcTotHitCnt;
        // ignore hits that were either unmatched to the selected set of PDGCodes
        // or have a not-requested origin
        if (fHitMCPIndex[iht] == UINT_MAX) continue;
        ++tpcMatHitCnt;
        unsigned int mcpi = fHitMCPIndex[iht];
        unsigned short mIndx = 0;
        for (mIndx = 0; mIndx < mcpCnt.size(); ++mIndx)
          if (mcpCnt[mIndx].first == mcpi) break;
        if (mIndx == mcpCnt.size()) {
          // found a MCParticle not in the list yet, so add it
          mcpCnt.push_back(std::make_pair(mcpi, 0));
          mcpClsCnt.resize(mcpCnt.size());
        }
        // increment the number of MC-matched hits
        ++mcpCnt[mIndx].second;
      } // iht
      // ignore TPCs/planes with few MC-matched hits
      if (fPrintLevel > 2) {
        mf::LogVerbatim myprt("ClusterTrackAna");
        myprt << "TPC:" << tpc << " Plane:" << plane << " has " << tpcMatHitCnt << "/"
              << tpcTotHitCnt;
        myprt << " (MC-matched hits) / (all hits)";
      }
      if (tpcMatHitCnt < 3) continue;
      // next iterate over all clusters/tracks and count mc-matched hits that are in this TPC/plane
      for (unsigned int ii = 0; ii < recoHits.size(); ++ii) {
        float nRecoHitsInPlane = 0;
        float nRecoMatHitsInPlane = 0;
        for (auto iht : recoHits[ii]) {
          auto& hit = (*fHitHandle)[iht];
          if (hit.WireID().TPC != tpc) continue;
          if (hit.WireID().Plane != plane) continue;
          ++nRecoHitsInPlane;
          if (fHitMCPIndex[iht] == UINT_MAX) continue;
          ++nRecoMatHitsInPlane;
          unsigned int mcpi = fHitMCPIndex[iht];
          // find the MCParticle index in mcpCnt and use it to count the match entry
          // in mcpClsCnt
          unsigned short mIndx = 0;
          for (mIndx = 0; mIndx < mcpCnt.size(); ++mIndx)
            if (mcpCnt[mIndx].first == mcpi) break;
          if (mIndx == mcpCnt.size()) {
            std::cout << "Logic error: fHitMCPIndex = " << fHitMCPIndex[iht]
                      << " is valid but isn't in the list of MCParticles to use. Please send email "
                         "to baller@fnal.gov.\n";
            continue;
          }
          unsigned short cIndx = 0;
          for (cIndx = 0; cIndx < mcpClsCnt[mIndx].size(); ++cIndx)
            if (mcpClsCnt[mIndx][cIndx].first == ii) break;
          if (cIndx == mcpClsCnt[mIndx].size()) mcpClsCnt[mIndx].push_back(std::make_pair(ii, 0));
          ++mcpClsCnt[mIndx][cIndx].second;
        } // cluhit
        if (nRecoMatHitsInPlane == 0) continue;
      } // ii
      // find the cluster that has the most hits matched to each MC Particle
      for (unsigned short mIndx = 0; mIndx < mcpCnt.size(); ++mIndx) {
        // require at least 3 MC-matched hits
        if (mcpCnt[mIndx].second < 3) continue;
        unsigned int mcpi = mcpCnt[mIndx].first;
        auto& mcp = (*mcpHandle)[mcpi];
        unsigned int pdgCode = abs(mcp.PdgCode());
        unsigned short pIndx = USHRT_MAX;
        if (pdgCode == 11) pIndx = 0;
        if (pdgCode == 13) pIndx = 1;
        if (pdgCode == 211) pIndx = 2;
        if (pdgCode == 321) pIndx = 3;
        if (pdgCode == 2212) pIndx = 4;
        if (mcpClsCnt[mIndx].empty()) {
          // Un-reconstructed MCParticle hits
          ++Cnts[pIndx];
          ++fNBadEP;
          if (fPrintLevel > 0) {
            mf::LogVerbatim myprt("ClusterTrackAna");
            myprt << " MCPI " << mcpi;
            int TMeV = 1000 * (mcp.E() - mcp.Mass());
            myprt << " " << TMeV << " MeV";
            std::string partName = std::to_string(pdgCode);
            if (pdgCode == 11) partName = "El";
            if (pdgCode == 13) partName = "Mu";
            if (pdgCode == 211) partName = "Pi";
            if (pdgCode == 311) partName = "Ka";
            if (pdgCode == 2212) partName = "Pr";
            myprt << std::setw(5) << partName;
            // print out the range of truth-matched hits
            unsigned int firstHitIndex = UINT_MAX;
            unsigned int lastHitIndex = UINT_MAX;
            FirstLastHitInPlane(tpc, plane, mcpi, firstHitIndex, lastHitIndex);
            myprt << " Failed to reconstruct. Truth-matched hit range from ";
            myprt << HitLocation(firstHitIndex);
            myprt << " to ";
            myprt << HitLocation(lastHitIndex);
            myprt << " <- EP = 0!";
          } // fPrintLevel > 0
          continue;
        } // (mcpClsCnt[mIndx].empty()
        std::pair<unsigned int, float> big = std::make_pair(UINT_MAX, 0);
        for (unsigned short cIndx = 0; cIndx < mcpClsCnt[mIndx].size(); ++cIndx) {
          auto& mcc = mcpClsCnt[mIndx][cIndx];
          if (mcc.second > big.second) big = mcc;
        } // cIndx
        if (big.first == UINT_MAX) {
          if (fPrintLevel > 2) {
            mf::LogVerbatim myprt("ClusterTrackAna");
            unsigned int mcpi = mcpCnt[mIndx].first;
            auto& mcp = (*mcpHandle)[mcpi];
            myprt << "Match failed: mcpi " << mcpi << " pdg " << mcp.PdgCode();
          }
          std::cout << "match failed mIndx " << mIndx << "\n";
          continue;
        } // big.first == UINT_MAX
        unsigned int ii = big.first;
        float eff = big.second / mcpCnt[mIndx].second;
        float nRecoHitsInPlane = 0;
        // define some variables to print the range of the cluster/track
        unsigned int firstRecoHitIndex = UINT_MAX;
        unsigned int lastRecoHitIndex = UINT_MAX;
        for (auto iht : recoHits[ii]) {
          auto& hit = (*fHitHandle)[iht];
          if (hit.WireID().TPC != tpc) continue;
          if (hit.WireID().Plane != plane) continue;
          ++nRecoHitsInPlane;
          if (firstRecoHitIndex == UINT_MAX) {
            firstRecoHitIndex = iht;
            lastRecoHitIndex = iht;
          }
          unsigned int wire = (*fHitHandle)[iht].WireID().Wire;
          if (wire < (*fHitHandle)[firstRecoHitIndex].WireID().Wire) firstRecoHitIndex = iht;
          if (wire > (*fHitHandle)[lastRecoHitIndex].WireID().Wire) lastRecoHitIndex = iht;
        } // iht
        float pur = big.second / nRecoHitsInPlane;
        ++Cnts[pIndx];
        float ep = eff * pur;
        EPSums[pIndx] += ep;
        ESums[pIndx] += eff;
        PSums[pIndx] += pur;
        bool hasBadEP = (ep < fBadEP);
        if (hasBadEP) ++fNBadEP;
        bool prt = fPrintLevel > 1 || (fPrintLevel == 1 && hasBadEP);
        if (prt) {
          mf::LogVerbatim myprt("ClusterTrackAna");
          if (fFirstPrint) {
            myprt << "Hit location format is TPC:Plane:Wire:Tick\n";
            myprt << " MCPI     Ptcl T(MeV)  nMCPHits  ____From____ _____To_____";
            if (inputClusters.isValid()) { myprt << "  clsID"; }
            else {
              myprt << "  trkID";
            }
            myprt
              << "  ____From____  _____To_____   nRecoHits nMCPRecoHits    Eff      Pur      EP\n";
            fFirstPrint = false;
          } // fFirstPrint
          myprt << std::setw(5) << mcpi;
          // convert the PDG code into nicer format
          std::string partName = std::to_string(pdgCode);
          if (pdgCode == 11) partName = " Electron";
          if (pdgCode == 13) partName = "     Muon";
          if (pdgCode == 211) partName = "     Pion";
          if (pdgCode == 311) partName = "     Kaon";
          if (pdgCode == 2212) partName = "   Proton";
          myprt << partName;
          int TMeV = 1000 * (mcp.E() - mcp.Mass());
          myprt << std::setw(7) << TMeV;
          // nMCPHits
          myprt << std::setw(10) << mcpCnt[mIndx].second;
          unsigned int firstTruHitIndex = UINT_MAX;
          unsigned int lastTruHitIndex = UINT_MAX;
          FirstLastHitInPlane(tpc, plane, mcpi, firstTruHitIndex, lastTruHitIndex);
          myprt << std::setw(14) << HitLocation(firstTruHitIndex);
          myprt << std::setw(14) << HitLocation(lastTruHitIndex);
          int id = -1;
          if (inputClusters.isValid()) {
            // print cluster info
            auto& cls = (*inputClusters)[recoIndex[ii]];
            id = cls.ID();
          }
          else if (inputTracks.isValid()) {
            // print track info
            auto& trk = (*inputTracks)[recoIndex[ii]];
            id = trk.ID();
          }
          myprt << std::setw(6) << id;
          myprt << std::setw(14) << HitLocation(firstRecoHitIndex);
          myprt << std::setw(14) << HitLocation(lastRecoHitIndex);
          myprt << std::setw(12) << (int)nRecoHitsInPlane;
          myprt << std::setw(13) << (int)big.second;
          myprt << std::fixed << std::setprecision(2);
          myprt << std::setw(8) << eff;
          myprt << std::setw(8) << pur;
          myprt << std::setw(8) << ep;
          if (hasBadEP) myprt << " <- BadEP";
          myprt << " Evt: " << evt.event();
          myprt << " Evt Cnt " << fEventCnt;
        } // prt
      }   // mIndx
    }     // plane
  }       // tpcid

  fHitMCPIndex.resize(0);

} // analyze

////////////////////////////////////////////////
std::string cluster::ClusterTrackAna::HitLocation(unsigned int iht)
{
  // Put the hit location into a compact human-readable format
  if (iht >= (*fHitHandle).size()) return "NA";
  auto& hit = (*fHitHandle)[iht];
  return std::to_string(hit.WireID().TPC) + ":" + std::to_string(hit.WireID().Plane) + ":" +
         std::to_string(hit.WireID().Wire) + ":" + std::to_string((int)hit.PeakTime());
} // HitLocation

////////////////////////////////////////////////
void cluster::ClusterTrackAna::FirstLastHitInPlane(unsigned int tpc,
                                                   unsigned int plane,
                                                   unsigned int mcpi,
                                                   unsigned int& firstHitIndex,
                                                   unsigned int& lastHitIndex)
{
  // Returns the index of the first hit (lowest wire number) and last hit (highest wire number)
  // matched to the MCParticle indexed by mcpi in the requested tpc, plane
  firstHitIndex = UINT_MAX;
  lastHitIndex = UINT_MAX;
  for (unsigned int iht = 0; iht < (*fHitHandle).size(); ++iht) {
    if (fHitMCPIndex[iht] != mcpi) continue;
    auto& hit = (*fHitHandle)[iht];
    if (hit.WireID().TPC != tpc) continue;
    if (hit.WireID().Plane != plane) continue;
    if (firstHitIndex == UINT_MAX) {
      firstHitIndex = iht;
      lastHitIndex = iht;
    }
    unsigned int wire = (*fHitHandle)[iht].WireID().Wire;
    if (wire < (*fHitHandle)[firstHitIndex].WireID().Wire) firstHitIndex = iht;
    if (wire > (*fHitHandle)[lastHitIndex].WireID().Wire) lastHitIndex = iht;
  } // iht
} // FirstLastHitInPlane

////////////////////////////////////////////////
void cluster::ClusterTrackAna::endJob()
{
  // output results
  mf::LogVerbatim myprt("ClusterTrackAna");
  myprt << "ClusterTrackAna summary results for " << fEventCnt;
  if (fClusterModuleLabel.label() != "NA") {
    myprt << " events using ClusterModuleLabel: " << fClusterModuleLabel.label();
  }
  else {
    myprt << " events using TrackModuleLabel: " << fTrackModuleLabel.label();
  }
  myprt << " Origin: " << fTruthOrigin;
  if (fInTPC != UINT_MAX) { myprt << " in TPC " << fInTPC; }
  else {
    myprt << " in all TPCs";
  }
  myprt << "\n";
  float cnts = 0;
  for (unsigned short pIndx = 0; pIndx < 5; ++pIndx)
    cnts += Cnts[pIndx];
  if (cnts == 0) {
    myprt << "No ClusterTrackAna results";
    return;
  }
  float sumEP = 0;
  float sumE = 0;
  float sumP = 0;
  myprt << "Efficiency (Eff), Purity (Pur) and Eff * Pur (EP) by selected truth particle types\n";
  std::array<std::string, 5> pName = {{"El", "Mu", "Pi", "K ", "P "}};
  myprt << "particle     Eff     Pur     EP\n";
  myprt << "-------------------------------";
  for (unsigned short pIndx = 0; pIndx < 5; ++pIndx) {
    if (Cnts[pIndx] == 0) continue;
    float ave;
    myprt << "\n   " << pName[pIndx] << "   ";
    myprt << std::fixed << std::setprecision(3);
    ave = ESums[pIndx] / Cnts[pIndx];
    myprt << std::setw(8) << ave;
    ave = PSums[pIndx] / Cnts[pIndx];
    myprt << std::setw(8) << ave;
    ave = EPSums[pIndx] / Cnts[pIndx];
    myprt << std::setw(8) << ave;
    if (pIndx == 0) continue;
    sumEP += EPSums[pIndx];
    sumE += ESums[pIndx];
    sumP += PSums[pIndx];
  } // pIndx
  if (cnts == 0) return;
  myprt << "\n";
  myprt << "Averages for all selected truth particles\n";
  myprt << "Ave Eff " << sumE / cnts;
  myprt << " Ave Pur " << sumP / cnts;
  myprt << " Ave EP " << sumEP / cnts;
  myprt << " nBadEP " << fNBadEP;
  myprt << " (EP < " << std::fixed << std::setprecision(2) << fBadEP << ")";
  myprt << "\n";
  myprt << "MCParticle counts in all TPCs and Planes:";
  for (unsigned short pIndx = 0; pIndx < 5; ++pIndx) {
    if (Cnts[pIndx] == 0) continue;
    myprt << " " << pName[pIndx] << " " << (int)Cnts[pIndx];
  } // pIndx
} // endJob
DEFINE_ART_MODULE(cluster::ClusterTrackAna)
