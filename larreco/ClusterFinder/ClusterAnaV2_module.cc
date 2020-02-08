////////////////////////////////////////////////////////////////////////
// Class:       ClusterAnaV2
// Plugin Type: analyzer (art v3_03_01)
// File:        ClusterAnaV2_module.cc
//
// Version 2 of ClusterAna that supports multi-TPC detectors
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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cluster {
  class ClusterAnaV2;
}


class cluster::ClusterAnaV2 : public art::EDAnalyzer {
public:
  explicit ClusterAnaV2(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ClusterAnaV2(ClusterAnaV2 const&) = delete;
  ClusterAnaV2(ClusterAnaV2&&) = delete;
  ClusterAnaV2& operator=(ClusterAnaV2 const&) = delete;
  ClusterAnaV2& operator=(ClusterAnaV2&&) = delete;

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
  short fInTPC;
  float fBadEP;

  // count of EP entries for electrons(0), muons(1), pions(2), kaons(3), protons(4)
  unsigned int fEventCnt {0};
  unsigned int fNBadEP {0};
  std::array<float, 5> Cnts {{0}};
  std::array<float, 5> EPSums {{0}};
  // same for Efficiency
  std::array<float, 5> ESums {{0}};
  // and for Purity
  std::array<float, 5> PSums {{0}};

  bool fCompareProductIDs {true};     ///< compare Hit and Cluster-> Hit art product IDs on the first event
  bool fFirstPrint {true};
};


cluster::ClusterAnaV2::ClusterAnaV2(fhicl::ParameterSet const& pset) : EDAnalyzer{pset}
{
  fHitModuleLabel = pset.get<art::InputTag>("HitModuleLabel");
  bool validModuleLabel = false;
  fClusterModuleLabel = "NA";
  fTrackModuleLabel = "NA";
  fClusterModuleLabel = pset.get<art::InputTag>("ClusterModuleLabel", "NA");
  if(fClusterModuleLabel != "NA") validModuleLabel = true;
  fTrackModuleLabel = pset.get<art::InputTag>("TrackModuleLabel", "NA");
  if(validModuleLabel && fTrackModuleLabel != "NA") throw cet::exception("ClusterAnaV2")<<"You must specify either a ClusterModuleLabel OR a TrackModuleLabel\n";
  if(!validModuleLabel && fTrackModuleLabel != "NA") validModuleLabel = true;
  // origin = 0 (anything), 1(nu), 2(cosmics), 3(SN nu), 4(SingleParticle)
  int tmp = pset.get<int>("TruthOrigin", 0);
  fTruthOrigin = (simb::Origin_t)tmp;
  fPrintLevel = pset.get<short>("PrintLevel", 0);
  if(pset.has_key("SkipPDGCodes")) fSkipPDGCodes = pset.get<std::vector<int>>("SkipPDGCodes");
    fBadEP = pset.get<float>("BadEP", 0.);
    fInTPC = pset.get<short>("InTPC", -1);
    // do some initialization
    Cnts.fill(0.);
    EPSums.fill(0.);
    ESums.fill(0.);
    PSums.fill(0.);
    } // ClusterAnaV2 constructor

////////////////////////////////////////////////
void cluster::ClusterAnaV2::analyze(art::Event const& evt)
{
  // Match hits to MCParticles, then consider reconstructed hits in each TPC and plane
  // to calculate Efficiency, Purity and Efficiency * Purity (aka EP).

  ++fEventCnt;
  auto const* geom = lar::providerFrom<geo::Geometry>();
  auto inputHits = art::Handle<std::vector<recob::Hit>>();
  if(!evt.getByLabel(fHitModuleLabel, inputHits)) throw cet::exception("ClusterAnaV2")<<"Failed to get a handle to hit collection '"<<fHitModuleLabel.label()<<"'\n";
//  unsigned int nInputHits = (*inputHits).size();
  // get a reference to the MCParticles
  auto mcpHandle = art::Handle<std::vector<simb::MCParticle>>();
  if(!evt.getByLabel("largeant", mcpHandle)) throw cet::exception("ClusterAnaV2")<<"Failed to get a handle to MCParticles using largeant\n";

  if(fFirstPrint) {
    mf::LogVerbatim("ClusterAna")<<"Reconstructed cluster hit range format is TPC:Plane:Wire:Tick";
    fFirstPrint = false;
  }

  // decide whether to consider cluster -> hit -> MCParticle for any MCParticle origin or for
  // a specific user-specified origin
  bool anySource = (fTruthOrigin == simb::kUnknown);

  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
  sim::ParticleList const& plist = pi_serv->ParticleList();
  
  bool firstPrt = true;

  // make a list of Hit -> MCParticle assns in all TPCs. The first step is
  // to make a list of Geant TrackIDs whose origin was requested by the user
  std::vector<int> trackIDs;
  for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
    const simb::MCParticle* part = (*ipart).second;
    art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(part->TrackId());
    if(!anySource && theTruth->Origin() != fTruthOrigin) continue;
    int pdg = abs(part->PdgCode());
    bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
    if(!isCharged) continue;
    if(std::find(fSkipPDGCodes.begin(), fSkipPDGCodes.end(), pdg) != fSkipPDGCodes.end()) continue;
    float TMeV = 1000 * (part->E() - part->Mass());
    // ignore low energy particles
    if(TMeV < 1) continue;
    trackIDs.push_back(part->TrackId());
    if(fPrintLevel > 3) {
      mf::LogVerbatim myprt("ClusterAna");
      myprt<<"TrackId "<<std::setw(8)<<part->TrackId();
      unsigned int mcpi = UINT_MAX;
      for(mcpi = 0; mcpi < (*mcpHandle).size(); ++mcpi) {
        auto& mcp = (*mcpHandle)[mcpi];
        if(mcp.TrackId() == part->TrackId()) break;
      } // indx
      myprt<<" mcpi "<<std::setw(8)<<mcpi;
      myprt<<" pdg "<<std::setw(8)<<pdg;
      myprt<<" origin "<<theTruth->Origin();
      myprt<<" T = "<<(int)TMeV<<" MeV "<<part->Process();
    }
  } // ipart
  if(trackIDs.empty()) {
    if(fPrintLevel > 0) std::cout<<"ClusterAnaV2: No MC particles found\n";
    return;
  }
  // next construct a companion vector of MCParticle indices indexed to the full hit collection
  std::vector<unsigned int> hitMCPIndex((*inputHits).size(), UINT_MAX);
  // Make a list of the <first, last> MC-matched hit in each TPC. This will be used
  // to only iterate through the range of hits that are interesting
  std::vector<std::pair<unsigned int, unsigned int>> hitRange(geom->NTPC() + 1);
  for(auto& hr : hitRange) hr = std::make_pair(UINT_MAX, UINT_MAX);
  for(size_t iht = 0; iht < (*inputHits).size(); ++iht) {
    auto& hit = (*inputHits)[iht];
    // only consider hits in a single TPC?
    if(fInTPC >= 0 && hit.WireID().TPC != fInTPC) continue;
    auto tides = bt_serv->HitToTrackIDEs(hit);
    if(tides.empty()) continue;
    bool gottaMatch = false;
    for(auto& tide : tides) {
      // declare a match to a MCParticle if > 50% of energy is from it
      if(tide.energyFrac < 0.5) continue;
      int trackID = tide.trackID;
      // see if this is a MCParticle we want to use
      if(std::find(trackIDs.begin(), trackIDs.end(), trackID) == trackIDs.end()) continue;
      // find the MCParticle index and define the Hit -> MCParticle assn
      for(unsigned int mcpi = 0; mcpi < (*mcpHandle).size(); ++mcpi) {
        auto& mcp = (*mcpHandle)[mcpi];
        if(mcp.TrackId() != trackID) continue;
        hitMCPIndex[iht] = mcpi;
        break;
      } // indx
      gottaMatch = true;
      break;
    } // tide
    if(!gottaMatch) continue;
    // populate hitRange
    unsigned short tpc = hit.WireID().TPC;
    if(tpc >= hitRange.size()) throw cet::exception("ClusterAnaV2")<<"Invalid hit TPC\n";
    if(hitRange[tpc].first == UINT_MAX) hitRange[tpc].first = iht;
    hitRange[tpc].second = iht;
  } // iht

  // fill a vector of hits indices from all clusters or all tracks
  std::vector<std::vector<unsigned int>> recoHits;
  // and a companion vector of indices into the cluster or track collection
  std::vector<unsigned int> recoIndex;
  // handles to these collections to enable printing the cluster or track ID. The
  // handle will be valid for either clusters or tracks
  art::Handle<std::vector<recob::Cluster>> inputClusters;
  art::Handle<std::vector<recob::Track>> inputTracks;
  if(fClusterModuleLabel.label() != "NA") {
    // get hits from clusters
    if(!evt.getByLabel(fClusterModuleLabel,inputClusters)) throw cet::exception("ClusterAnaV2")<<"Failed to get a handle to the cluster collection '"<<fClusterModuleLabel.label()<<"'\n";
    art::FindManyP<recob::Hit> hitsFromCls(inputClusters, evt, fClusterModuleLabel);
    if(!hitsFromCls.isValid()) throw cet::exception("ClusterAnaV2")<<"Failed to get a handle to Cluster -> Hit assns\n";
    for(unsigned int icl = 0; icl < (*inputClusters).size(); ++icl) {
      std::vector<art::Ptr<recob::Hit>> cluhits = hitsFromCls.at(icl);
      if(cluhits.empty()) continue;
      if(fCompareProductIDs) {
        if(cluhits[0].id() != inputHits.id()) throw cet::exception("ClusterAnaV2")<<"Hits associated with ClusterModuleLabel are in a different collection than HitModuleLabel.\n";
        fCompareProductIDs = false;
      } // fCompareProductIDs
      recoIndex.push_back(icl);
      std::vector<unsigned int> hitsIndex;
      for(auto& cluhit : cluhits) hitsIndex.push_back(cluhit.key());
      recoHits.push_back(hitsIndex);
    } // icl
  } else {
    // get hits from tracks
    if(!evt.getByLabel(fTrackModuleLabel, inputTracks)) throw cet::exception("ClusterAnaV2")<<"Failed to get a handle to the track collection '"<<fTrackModuleLabel.label()<<"'\n";
    art::FindManyP<recob::Hit> hitsFromTrk(inputTracks, evt, fTrackModuleLabel);
    if(!hitsFromTrk.isValid()) throw cet::exception("ClusterAnaV2")<<"Failed to get a handle to Track -> Hit assns\n";
    for(unsigned int itk = 0; itk < (*inputTracks).size(); ++itk) {
      std::vector<art::Ptr<recob::Hit>> trkhits = hitsFromTrk.at(itk);
      if(trkhits.empty()) continue;
      if(fCompareProductIDs) {
        if(trkhits[0].id() != inputHits.id()) throw cet::exception("ClusterAnaV2")<<"Hits associated with TrackModuleLabel are in a different collection than HitModuleLabel.\n";
        fCompareProductIDs = false;
      } // fCompareProductIDs
      recoIndex.push_back(itk);
      std::vector<unsigned int> hitsIndex;
      for(auto& trkhit : trkhits) hitsIndex.push_back(trkhit.key());
      recoHits.push_back(hitsIndex);
    } // icl
  }// get hits from tracks
  if(recoHits.empty()) {
    std::cout<<"recoHits is empty. Does this make sense with "<<(*inputHits).size()<<" hits?";
    return;
  }

  if(fPrintLevel > 1) mf::LogVerbatim("ClusterAnaV2")<<"loaded "<<recoHits.size()<<" recoHits collections";

  for(const auto& tpcid : geom->IterateTPCIDs()) {
    unsigned int tpc = tpcid.TPC;
    if(hitRange[tpc].first == UINT_MAX) continue;
    // iterate over planes
    for(unsigned short plane = 0; plane < geom->Nplanes(); ++plane) {
      unsigned int tpcHitCnt = 0;
      // create a list of (MCParticle index, matched hit count> pairs
      //  mcParticle  plane
      std::vector<std::pair<unsigned int, float>> mcpCnt;
      // count MCParticle, Cluster/Track hit counts - size matched to mcpCnt
      std::vector<std::vector<std::pair<unsigned int, float>>> mcpClsCnt;
      for(unsigned int iht = hitRange[tpc].first; iht <= hitRange[tpc].second; ++iht) {
        // ignore unmatched or not-requested-origin hits
        if(hitMCPIndex[iht] == UINT_MAX) continue;
        auto& hit = (*inputHits)[iht];
        if(hit.WireID().TPC != tpc) continue;
        if(hit.WireID().Plane != plane) continue;
        ++tpcHitCnt;
        unsigned int mcpi = hitMCPIndex[iht];
        unsigned short mIndx = 0;
        for(mIndx = 0; mIndx < mcpCnt.size(); ++mIndx) if(mcpCnt[mIndx].first == mcpi) break;
        if(mIndx == mcpCnt.size()) {
          // found a MCParticle not in the list yet, so add it
          mcpCnt.push_back(std::make_pair(mcpi, 0));
          mcpClsCnt.resize(mcpCnt.size());
        }
        // increment the number of MC-matched hits
        ++mcpCnt[mIndx].second;
      } // iht
      // ignore TPCs/planes with few MC-matched hits
      if(tpcHitCnt < 3) continue;
      if(fPrintLevel > 1) mf::LogVerbatim("ClusterAna")<<"TPC:"<<tpc<<" Plane:"<<plane<<" has "<<tpcHitCnt<<" MC-matched hits";
      // next iterate over all clusters/tracks and count mc-matched hits that are in this TPC/plane
      for(unsigned int ii = 0; ii < recoHits.size(); ++ii) {
        float nRecoHitsInPlane = 0;
        float nRecoMatHitsInPlane = 0;
        for(auto iht : recoHits[ii]) {
          auto& hit = (*inputHits)[iht];
          if(hit.WireID().TPC != tpc) continue;
          if(hit.WireID().Plane != plane) continue;
          ++nRecoHitsInPlane;
          if(hitMCPIndex[iht] == UINT_MAX) continue;
          ++nRecoMatHitsInPlane;
          unsigned int mcpi = hitMCPIndex[iht];
          // find the MCParticle index in mcpCnt and use it to count the match entry
          // in mcpClsCnt
          unsigned short mIndx = 0;
          for(mIndx = 0; mIndx < mcpCnt.size(); ++mIndx) if(mcpCnt[mIndx].first == mcpi) break;
          if(mIndx == mcpCnt.size()) {
            std::cout<<"Logic error: hitMCPIndex = "<<hitMCPIndex[iht]<<" is valid but isn't in the list of MCParticles to use. Please send email to baller@fnal.gov.\n";
            continue;
          }
          unsigned short cIndx = 0;
          for(cIndx = 0; cIndx < mcpClsCnt[mIndx].size(); ++cIndx)
            if(mcpClsCnt[mIndx][cIndx].first == ii) break;
          if(cIndx == mcpClsCnt[mIndx].size()) mcpClsCnt[mIndx].push_back(std::make_pair(ii, 0));
          ++mcpClsCnt[mIndx][cIndx].second;
        } // cluhit
        if(nRecoMatHitsInPlane == 0) continue;
      } // ii
      // find the cluster that has the most hits matched to each MC Particle
      for(unsigned short mIndx = 0; mIndx < mcpCnt.size(); ++mIndx) {
        // require at least 3 MC-matched hits
        if(mcpCnt[mIndx].second < 3) continue;
        unsigned int mcpi = mcpCnt[mIndx].first;
        auto& mcp = (*mcpHandle)[mcpi];
        unsigned int pdgCode = abs(mcp.PdgCode());
        unsigned short pIndx = USHRT_MAX;
        if(pdgCode == 11) pIndx = 0;
        if(pdgCode == 13) pIndx = 1;
        if(pdgCode == 211) pIndx = 2;
        if(pdgCode == 321) pIndx = 3;
        if(pdgCode == 2212) pIndx = 4;
        if(mcpClsCnt[mIndx].empty()) {
          // Un-reconstructed MCParticle hits
          ++Cnts[pIndx];
          if(fPrintLevel > 0) {
            mf::LogVerbatim myprt("ClusterAna");
            if(firstPrt) {
              myprt<<"Run: "<<evt.run();
              myprt<<" Event: "<<evt.event()<<"\n";
              firstPrt = false;
            }
            myprt<<" MCPI "<<mcpi<<" PDG Code "<<pdgCode;
            myprt<<" Failed to reconstruct in plane "<<plane<<". Truth-matched hit range from ";
            // print out the range of truth-matched hits
            unsigned int firstHitIndex = UINT_MAX;
            unsigned int lastHitIndex = UINT_MAX;
            for(unsigned int iht = 0; iht < (*inputHits).size(); ++iht) {
              if(hitMCPIndex[iht] != mcpi) continue;
              auto& hit = (*inputHits)[iht];
              if(hit.WireID().TPC != tpc) continue;
              if(hit.WireID().Plane != plane) continue;
              if(firstHitIndex == UINT_MAX) firstHitIndex = iht;
              lastHitIndex = iht;
            } // iht
            if((*inputHits)[firstHitIndex].WireID().Wire > (*inputHits)[lastHitIndex].WireID().Wire) std::swap(firstHitIndex, lastHitIndex);
            auto& fHit = (*inputHits)[firstHitIndex];
            myprt<<fHit.WireID().TPC<<":"<<fHit.WireID().Plane<<":"<<fHit.WireID().Wire;
            myprt<<(int)fHit.PeakTime();
            auto& lHit = (*inputHits)[lastHitIndex];
            myprt<<" to ";
           myprt<<lHit.WireID().TPC<<":"<<lHit.WireID().Plane<<":"<<lHit.WireID().Wire;
           myprt<<(int)lHit.PeakTime();
           myprt<<" <- EP = 0!";
         } // fPrintLevel > 0
          continue;
        } // (mcpClsCnt[mIndx].empty()
        std::pair<unsigned int, float> big = std::make_pair(UINT_MAX, 0);
        for(unsigned short cIndx = 0; cIndx < mcpClsCnt[mIndx].size(); ++cIndx) {
          auto& mcc = mcpClsCnt[mIndx][cIndx];
          if(mcc.second > big.second) big = mcc;
        } // cIndx
        if(big.first == UINT_MAX) {
          if(fPrintLevel > 1) {
            mf::LogVerbatim myprt("ClusterAna");
            unsigned int mcpi = mcpCnt[mIndx].first;
            auto& mcp = (*mcpHandle)[mcpi];
            myprt<<"Match failed: mcpi "<<mcpi<<" pdg "<<mcp.PdgCode();
          }
          std::cout<<"match failed mIndx "<<mIndx<<"\n";
          continue;
        } // big.first == UINT_MAX
        unsigned int ii = big.first;
        float eff = big.second / mcpCnt[mIndx].second;
        float nRecoHitsInPlane = 0;
        // define some variables to print the range of the cluster/track
        unsigned int firstHitIndex = UINT_MAX;
        unsigned int lastHitIndex = UINT_MAX;
        for(auto iht : recoHits[ii]) {
          auto& hit = (*inputHits)[iht];
          if(hit.WireID().TPC != tpc) continue;
          if(hit.WireID().Plane != plane) continue;
          if(firstHitIndex == UINT_MAX) firstHitIndex = iht;
          lastHitIndex = iht;
          ++nRecoHitsInPlane;
        } // iht
        if((*inputHits)[firstHitIndex].WireID().Wire >
           (*inputHits)[lastHitIndex].WireID().Wire) std::swap(firstHitIndex, lastHitIndex);
        float pur = big.second / nRecoHitsInPlane;
        ++Cnts[pIndx];
        float ep = eff * pur;
        EPSums[pIndx] += ep;
        ESums[pIndx] += eff;
        PSums[pIndx] += pur;
        bool hasBadEP = (ep < fBadEP);
        if(hasBadEP) ++fNBadEP;
        if(fPrintLevel > 0 || hasBadEP) {
          mf::LogVerbatim myprt("ClusterAna");
          if(firstPrt) {
            myprt<<"Run: "<<evt.run();
            myprt<<" Event: "<<evt.event()<<"\n";
            firstPrt = false;
          }
          myprt<<" MCPI "<<mcpi;
          int TMeV = 1000 * (mcp.E() - mcp.Mass());
          myprt<<" T "<<TMeV<<" MeV";
          if(inputClusters.isValid()) {
            // print cluster info
            auto& cls = (*inputClusters)[recoIndex[ii]];
            myprt<<" cls ID "<<cls.ID();
          } else if(inputTracks.isValid()) {
            // print track info
            auto& trk = (*inputTracks)[recoIndex[ii]];
            myprt<<" trk ID "<<trk.ID();
          }
          myprt<<" nMCPHits "<<mcpCnt[mIndx].second;
          myprt<<" nRecoHitsInPlane "<<nRecoHitsInPlane<<" nMCPRecoHits "<<big.second;
          myprt<<" eff "<<std::fixed<<std::setprecision(2)<<eff;
          myprt<<" pur "<<pur;
          auto& fHit = (*inputHits)[firstHitIndex];
          auto& lHit = (*inputHits)[lastHitIndex];
          myprt<<" from ";
          myprt<<fHit.WireID().TPC<<":"<<fHit.WireID().Plane<<":"<<fHit.WireID().Wire<<":"<<(int)fHit.PeakTime();
          myprt<<" to ";
          myprt<<lHit.WireID().TPC<<":"<<lHit.WireID().Plane<<":"<<lHit.WireID().Wire<<":"<<(int)lHit.PeakTime();
          if(hasBadEP) myprt<<" <- BadEP";
        } // fPrintLevel > 1
      } // mIndx
    } // plane
  } // tpcid

} // analyze

////////////////////////////////////////////////
void cluster::ClusterAnaV2::endJob()
{
  // output results
  mf::LogVerbatim myprt("ClusterAna");
  myprt<<"ClusterAnaV2 summary results for "<<fEventCnt;
  if(fClusterModuleLabel.label() != "NA"){
    myprt<<" events using ClusterModuleLabel: "<<fClusterModuleLabel.label();
  } else {
    myprt<<" events using TrackModuleLabel: "<<fTrackModuleLabel.label();
  }
  myprt<<" Origin: "<<fTruthOrigin;
  if(fInTPC >= 0) {
    myprt<<" in TPC "<<fInTPC;
  } else {
    myprt<<" in all TPCs";
  }
  myprt<<"\n";
  float cnts = 0;
  float sumEP = 0;
  float sumE = 0;
  float sumP = 0;
  myprt<<"Efficiency (Eff), Purity (Pur) and Eff * Pur (EP) by selected truth particle types\n";
  std::array<std::string, 5> pName = {{"El", "Mu", "Pi", "K ", "P "}};
  myprt<<"particle     Eff     Pur     EP\n";
  myprt<<"-------------------------------";
  for(unsigned short pIndx = 0; pIndx < 5; ++pIndx) {
    if(Cnts[pIndx] == 0) continue;
    float ave;
    myprt<<"\n   "<<pName[pIndx]<<"   ";
    myprt<<std::fixed<<std::setprecision(3);
    ave = ESums[pIndx] / Cnts[pIndx];
    myprt<<std::setw(8)<<ave;
    ave = PSums[pIndx] / Cnts[pIndx];
    myprt<<std::setw(8)<<ave;
    ave = EPSums[pIndx] / Cnts[pIndx];
    myprt<<std::setw(8)<<ave;
    if(pIndx == 0) continue;
    sumEP += EPSums[pIndx];
    sumE += ESums[pIndx];
    sumP += PSums[pIndx];
    cnts += Cnts[pIndx];
  } // pIndx
  if(cnts == 0) return;
  myprt<<"\n";
  myprt<<"Averages for all selected truth particles\n";
  myprt<<" Ave Eff "<<sumE/cnts;
  myprt<<" Ave Pur "<<sumP/cnts;
  myprt<<" Ave EP "<<sumEP/cnts;
  myprt<<" nBadEP "<<fNBadEP;
  myprt<<"\n";
  myprt<<" Cnts";
  for(unsigned short pIndx = 0; pIndx < 5; ++pIndx) {
    if(Cnts[pIndx] == 0) continue;
    myprt<<" "<<pName[pIndx]<<" "<<(int)Cnts[pIndx];
  } // pIndx
} // endJob
DEFINE_ART_MODULE(cluster::ClusterAnaV2)
