////////////////////////////////////////////////////////////////////////
//
// ClusterAna class
//
// Bruce Baller
//
////////////////////////////////////////////////////////////////////////

#include <TProfile.h>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>

//Framework includes
#include "fhiclcpp/fwd.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Core/EDAnalyzer.h"

/// Cluster finding and building
namespace cluster {

    // MCParticle -> cluster matching struct and indices
    struct MatchStruct {
      unsigned int mcpi {UINT_MAX};       ///< MCParticle index
      unsigned int tpc {UINT_MAX};        ///< TPC in which this MatchStruct is used
      std::vector<unsigned int> clsIndex; ///< index of the MC-matched cluster in each plane
      std::vector<float> mcpTruHitCount;  ///< number of MC-Matched hits to mcpi in each plane
      std::vector<float> clsTruHitCount;  ///< number of MC-matched hits to mcpi in clsIndex in each plane
      std::vector<unsigned int> firstHit; ///< first MC-matched hit in each plane
      std::vector<unsigned int> lastHit;  ///< last MC-matched hit in each plane
    };

  class ClusterAna : public art::EDAnalyzer {
  public:
    explicit ClusterAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void beginJob() override;
    void endJob() override;
    std::string PrintHit(const recob::Hit& hit);
    int PDGCodeIndex(int pdgCode);
    void FindFirstLastWire(art::Handle<std::vector<recob::Hit>> const& allHits, 
                           std::vector<unsigned int> const& hitMCPIndex,
                           MatchStruct const& match,
                           unsigned short inPlane, unsigned int& first, unsigned int& last);

    TH1F* fNClusters;
    TH1F* fNHitInCluster;
    // Cosmic Rays
    TH1F* fCREP2;
    TH1F* fCRE;
    TH1F* fCRP;
    // Neutrino interactions
    TH1F* fT_elec;
    TH1F* fT_muon;
    TH1F* fT_pion;
    TH1F* fT_kaon;
    TH1F* fT_prot;
    //  TH1F* fEP;
    TH1F* fEP_elec;
    TH1F* fEP_muon;
    TH1F* fEP_pion;
    TH1F* fEP_kaon;
    TH1F* fEP_prot;
    TH1F* fE_elec;
    TH1F* fE_muon;
    TH1F* fE_pion;
    TH1F* fE_kaon;
    TH1F* fE_prot;
    TH1F* fP_elec;
    TH1F* fP_muon;
    TH1F* fP_pion;
    TH1F* fP_kaon;
    TH1F* fP_prot;

    TH1F* fNuVtx_dx;
    TH1F* fNuVtx_dy;
    TH1F* fNuVtx_dz;

    TProfile* fEP_T_elec;
    TProfile* fEP_T_muon;
    TProfile* fEP_T_pion;
    TProfile* fEP_T_kaon;
    TProfile* fEP_T_prot;

    art::InputTag fClusterModuleLabel;
    art::ProductID fHitCollectionProductID;
    std::array<float, 2> fElecKERange;
    std::array<float, 2> fMuonKERange;
    std::array<float, 2> fPionKERange;
    std::array<float, 2> fKaonKERange;
    std::array<float, 2> fProtKERange;
    short fTrackWeightOption;
    bool fMergeDaughters;
    bool fSkipCosmics;
    short fPrintLevel;
    short moduleID;

    std::array<std::string, 5> fNames {{"electron", "muon", "pion", "kaon", "proton"}};
    std::array<float, 5> fEffSum {{0}};
    std::array<float, 5> fPurSum {{0}};
    std::array<float, 5> fEffPurSum {{0}};
    std::array<float, 5> fSum {{0}};
    unsigned int fNBadEP {0};
    unsigned int fEventsProcessed {0};

  }; // class ClusterAna

}

namespace cluster {

  //--------------------------------------------------------------------
  ClusterAna::ClusterAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fClusterModuleLabel(pset.get<art::InputTag>("ClusterModuleLabel"))
    , fElecKERange(pset.get<std::array<float,2>>("ElecKERange"))
    , fMuonKERange(pset.get<std::array<float,2>>("MuonKERange"))
    , fPionKERange(pset.get<std::array<float,2>>("PionKERange"))
    , fKaonKERange(pset.get<std::array<float,2>>("KaonKERange"))
    , fProtKERange(pset.get<std::array<float,2>>("ProtKERange"))
    , fSkipCosmics(pset.get<bool>("SkipCosmics", true))
    , fPrintLevel(pset.get<short>("PrintLevel"))
  {

  }

  //------------------------------------------------------------------
  void
  ClusterAna::beginJob()
  {

    if (fPrintLevel > 0) {
      mf::LogVerbatim myprt("ClusterAna");
      myprt << "ClusterAna: MCParticle selection";
      if (fSkipCosmics) myprt << ": ignoring Cosmics.";
      if(fPrintLevel > 1) myprt << " Note: Hit format is <TPC>:<Plane>:<Wire>:<PeakTime>";
    }

    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    fT_elec = tfs->make<TH1F>("T_elec", "T(MeV) electron", 500, 0, 2000);
    fT_muon = tfs->make<TH1F>("T_muon", "T(MeV) muon", 500, 0, 2000);
    fT_pion = tfs->make<TH1F>("T_pion", "T(MeV) pion", 500, 0, 2000);
    fT_kaon = tfs->make<TH1F>("T_kaon", "T(MeV) kaon", 500, 0, 2000);
    fT_prot = tfs->make<TH1F>("T_prot", "T(MeV) proton", 500, 0, 2000);

    fEP_elec = tfs->make<TH1F>("EP_elec", "EP electron", 50, 0, 1);
    fEP_muon = tfs->make<TH1F>("EP_muon", "EP muon", 50, 0, 1);
    fEP_pion = tfs->make<TH1F>("EP_pion", "EP pion", 50, 0, 1);
    fEP_kaon = tfs->make<TH1F>("EP_kaon", "EP kaon", 50, 0, 1);
    fEP_prot = tfs->make<TH1F>("EP_prot", "EP proton", 50, 0, 1);

    fE_elec = tfs->make<TH1F>("E_elec", "Efficiency electron", 50, 0, 1);
    fE_muon = tfs->make<TH1F>("E_muon", "Efficiency muon", 50, 0, 1);
    fE_pion = tfs->make<TH1F>("E_pion", "Efficiency pion", 50, 0, 1);
    fE_kaon = tfs->make<TH1F>("E_kaon", "Efficiency kaon", 50, 0, 1);
    fE_prot = tfs->make<TH1F>("E_prot", "Efficiency proton", 50, 0, 1);

    fP_elec = tfs->make<TH1F>("P_elec", "Purity electron", 50, 0, 1);
    fP_muon = tfs->make<TH1F>("P_muon", "Purity muon", 50, 0, 1);
    fP_pion = tfs->make<TH1F>("P_pion", "Purity pion", 50, 0, 1);
    fP_kaon = tfs->make<TH1F>("P_kaon", "Purity kaon", 50, 0, 1);
    fP_prot = tfs->make<TH1F>("P_prot", "Purity proton", 50, 0, 1);

    fEP_T_elec = tfs->make<TProfile>("EP_T_elec", "EP electron vs T", 200, 0, 2000);
    fEP_T_muon = tfs->make<TProfile>("EP_T_muon", "EP muon vs T", 200, 0, 2000);
    fEP_T_pion = tfs->make<TProfile>("EP_T_pion", "EP pion vs T", 200, 0, 2000);
    fEP_T_kaon = tfs->make<TProfile>("EP_T_kaon", "EP kaon vs T", 200, 0, 2000);
    fEP_T_prot = tfs->make<TProfile>("EP_T_prot", "EP proton vs T", 200, 0, 2000);
  }

  void
  ClusterAna::endJob()
  {
    mf::LogVerbatim myprt("ClusterAna");
    myprt<<"ClusterAna results for " <<fEventsProcessed << " events. MCParticle counts:";
    for(unsigned short indx = 0; indx < 5; ++indx) {
      if (fSum[indx] == 0) continue;
      myprt <<" " << fNames[indx] <<": " << (int)fSum[indx];
    } // indx
    myprt<<"\n";
    float totSum = 0;
    float totEPSum = 0;
    myprt << fClusterModuleLabel.label() << " EP:";
    for(unsigned short indx = 0; indx < 5; ++indx) {
      if (fSum[indx] == 0) continue;
      float aveEP = fEffPurSum[indx] / fSum[indx];
      myprt << " " << fNames[indx] <<" " << std::fixed << std::setprecision(2) << aveEP;
      totSum += fSum[indx];
      totEPSum += fEffPurSum[indx];
    } // indx
    if (totSum > 0) {
      totEPSum /= totSum;
      myprt << " AveEP: " << totEPSum;
      myprt << " nBad/nMCPInPln " << fNBadEP << "/" << (int)totSum;
    }
  } // endJob

  void
  ClusterAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()) return;

    // get a handle to the MCParticle collection
    art::InputTag mcpLabel = "largeant";
    auto mcps = art::Handle<std::vector<simb::MCParticle>>();
    if (!evt.getByLabel(mcpLabel, mcps)) return;
    if ((*mcps).empty()) return;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    // limit the number of print output lines
    unsigned short nPrtMCP = 0;
    unsigned short nPrtSel = 0;
    // temporary list of MCParticle to be used
    std::vector<unsigned int> mcpiList;
    for(unsigned int mcpi = 0; mcpi < (*mcps).size(); ++mcpi) {
      auto& mcp = (*mcps)[mcpi];
      int tid = mcp.TrackId();
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(tid);
      if (fSkipCosmics && theTruth->Origin() == simb::kCosmicRay) continue;
      float TMeV = 1000 * (mcp.E() - mcp.Mass());
      // consider electrons, muons, pions, kaons and protons
      int pdgIndx = PDGCodeIndex(mcp.PdgCode());
      if (pdgIndx < 0) continue;
      bool useIt = true;
      if (pdgIndx == 0) {
        if (fElecKERange[0] < 0 || TMeV < fElecKERange[0] || TMeV > fElecKERange[1]) useIt = false;
      } else if (pdgIndx == 1) {
        if (fMuonKERange[0] < 0 || TMeV < fMuonKERange[0] || TMeV > fMuonKERange[1]) useIt = false;
      } else if (pdgIndx == 2) {
        if (fPionKERange[0] < 0 || TMeV < fPionKERange[0] || TMeV > fPionKERange[1]) useIt = false;
      } else if (pdgIndx == 3) {
        if (fKaonKERange[0] < 0 || TMeV < fKaonKERange[0] || TMeV > fKaonKERange[1]) useIt = false;
      } else if (pdgIndx == 4) {
        if (fProtKERange[0] < 0 || TMeV < fProtKERange[0] || TMeV > fProtKERange[1]) useIt = false;
      }
      if (fPrintLevel > 2) {
        mf::LogVerbatim myprt("ClusterAna");
        if(nPrtMCP < 500) {
          ++nPrtMCP;
          myprt << "TrackId " << mcp.TrackId();
          myprt << ", PDG code " << mcp.PdgCode();
          myprt << ", T " << (int)TMeV << " MeV";
          myprt << ", Mother " << mcp.Mother();
          myprt << ", Process " << mcp.Process();
          if (useIt) myprt<<" <<< useIt";
        } else {
          myprt << "--> truncated print of " << (*mcps).size() << " MCParticles";
        }
      } // fPrintLevel > 2
      if (useIt) mcpiList.push_back(mcpi);
    } // mcpi
    if (mcpiList.empty()) return;

    // get clusters and cluster-hit associations
    art::Handle<std::vector<recob::Cluster>> allCls;
    if (!evt.getByLabel(fClusterModuleLabel, allCls))
      throw cet::exception("ClusterAna")<<"Failed to get a handle to cluster collection '"
      <<fClusterModuleLabel.label()<<"'\n";
    art::FindManyP<recob::Hit> fmch(allCls, evt, fClusterModuleLabel);
    // find the hit collection product ID on the first event 
    // from the first cluster hit (if clusters exist)
    if(!fHitCollectionProductID.isValid() && !(*allCls).empty())
        fHitCollectionProductID = fmch.at(0)[0].id();
    if(!fHitCollectionProductID.isValid()) return;
    // get a handle to the correct hit collection (allHits)
    auto allHits = art::Handle<std::vector<recob::Hit>>();
    // by checking the art product ID of all hit collections
    std::vector<art::Handle<std::vector<recob::Hit>>> hitHandles;
    evt.getManyByType(hitHandles);
    for (auto& hitHandle : hitHandles) {
      if(hitHandle.id() == fHitCollectionProductID) {
        allHits = hitHandle;
        break;
      } // product id's match
    } // hitHandle
    if(!allHits.isValid()) {
      throw cet::exception("ClusterAna")
          << "Failed to get a handle to the hits collection";
    }
    // size the vector of Hit -> MCParticle index
    if((*allHits).size() < 3) return;
    std::vector<unsigned int> hitMCPIndex((*allHits).size(), UINT_MAX);

    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<geo::Geometry const> geom;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Create/populate a MatchStruct in each TPC if there is a MC-matched hit in it
    std::vector<MatchStruct> matches;
    // next match the hits and populate the hitMCPIndex vector
    unsigned int nMatched = 0;
    for(unsigned int iht = 0; iht < (*allHits).size(); ++iht) {
      auto& hit = (*allHits)[iht];
      auto tides = bt_serv->HitToTrackIDEs(clockData, hit);
      for(auto& tide : tides) {
        int tid = tide.trackID;
        if (tide.energyFrac > 0.5) {
          for(auto mcpi : mcpiList) {
            if((*mcps)[mcpi].TrackId() != tid) continue;
            hitMCPIndex[iht] = mcpi;
            ++nMatched;
            break;
          }
        } // tide.energyFrac > 0.5
      } // tide
      if(hitMCPIndex[iht] == UINT_MAX) continue;
      // Update or create the appropriate MatchStruct
      unsigned int tpc = hit.WireID().TPC;
      unsigned int indx = 0;
      for(indx = 0; indx < matches.size(); ++indx) {
        auto& match = matches[indx];
        if(match.mcpi == hitMCPIndex[iht] && match.tpc == tpc) break;
      } // indx
      // Create a new matches entry if one wasn't found
      if(indx == matches.size()) {
        MatchStruct newMatch;
        newMatch.mcpi = hitMCPIndex[iht];
        newMatch.tpc = tpc;
        newMatch.clsIndex.resize(geom->Nplanes(), UINT_MAX);
        newMatch.clsTruHitCount.resize(geom->Nplanes(), 0);
        newMatch.mcpTruHitCount.resize(geom->Nplanes(), 0);
        newMatch.firstHit.resize(geom->Nplanes(), UINT_MAX);
        newMatch.lastHit.resize(geom->Nplanes(), UINT_MAX);
        matches.push_back(newMatch);
      } // not found
      // update the entry
      auto& match = matches[indx];
      unsigned short plane = hit.WireID().Plane;
      // note that firstHit (lastHit) is the first (last) hit in the
      // hit collection - not the first (last) in the plane
      if (match.firstHit[plane] == UINT_MAX) match.firstHit[plane] = iht;
      match.lastHit[plane] = iht;
      ++match.mcpTruHitCount[plane];
    } // iht
    if (fPrintLevel > 1) mf::LogVerbatim("ClusterAna")<<"Matched "<<nMatched<<" hits to "
            << matches.size() << " selected MCParticles ";
    if(nMatched == 0) return;

    // Update the sum before returning if there were at least 3 MC-matched hits
    // in a plane and no cluster was reconstructed
    if ((*allCls).empty()) {
      for(auto& match : matches) {
        if (match.mcpi == UINT_MAX) continue;
        auto& mcp = (*mcps)[match.mcpi];
        int pdgIndx = PDGCodeIndex(mcp.PdgCode());
        for(unsigned short plane = 0; plane < geom->Nplanes(); ++plane) {
          if (match.mcpTruHitCount[plane] > 2) ++fSum[pdgIndx];
        } // plane
      } // match
      return;
    } // no clusters

    // Match MCParticles to Clusters in each TPC
    for(const auto& tpcid : geom->IterateTPCIDs()) {
      // ignore bogus DUNE TPCs
      if (geom->TPC(tpcid).DriftDistance() < 25.0) continue;
      unsigned short tpc = tpcid.TPC;
      // Find the best match of this cluster to one MCParticle in this TPC.
      // This code assumes that a cluster only resides in one plane
      unsigned short plane = 0;
      for (unsigned int icl = 0; icl < (*allCls).size(); ++icl) {
        // Temporary vector of (mcp index, MC-matched hit count) pairs for this cluster
        std::vector<std::pair<unsigned int, float>> mcpiCnts;
        auto clsHits = fmch.at(icl);
        for (auto& clsHit : clsHits) {
          unsigned int iht = clsHit.key();
          // ignore hits that aren't MC-matched
          if (hitMCPIndex[iht] == UINT_MAX) continue;
          auto& hit = (*allHits)[iht];
          // ignore hits that aren't in this tpc
          if(hit.WireID().TPC != tpc) continue;
          plane = hit.WireID().Plane;
          unsigned int mcpi = hitMCPIndex[iht];
          // look for an mcpiCnts entry
          unsigned short indx = 0;
          for (indx = 0; indx < mcpiCnts.size(); ++indx) if(mcpiCnts[indx].first == mcpi) break;
          if(indx == mcpiCnts.size()) mcpiCnts.push_back(std::make_pair(mcpi, 0));
          ++mcpiCnts[indx].second;
        } // clsHit
        if (mcpiCnts.empty()) continue;
        // Find the MCParticle that has the most matched hits to this cluster
        float mostHitsMatched = 0;
        unsigned int bestMCPIndex = 0;
        for(auto& mcpiCnt : mcpiCnts) {
          if (mcpiCnt.second < mostHitsMatched) continue;
          mostHitsMatched = mcpiCnt.second;
          bestMCPIndex = mcpiCnt.first;
        } // mcpCnt
        // compare the number of matched hits for this cluster with an existing
        // MCParticle -> cluster match
        for(auto& match : matches) {
          if (match.tpc != tpc) continue;
          if (match.mcpi != bestMCPIndex) continue;
          // Found an existing match. Are there more matched hits?
          if (match.clsTruHitCount[plane] > mostHitsMatched) continue;
          // a better match
          match.clsTruHitCount[plane] = mostHitsMatched;
          match.clsIndex[plane] = icl;
          break;
        } // match
      } // icl
    } // tpc
    if (fPrintLevel > 1) {
      mf::LogVerbatim myprt("ClusterAna");
      if(nPrtSel < 500) {
        ++nPrtSel;
        for(auto& match : matches) {
          if (match.mcpi == UINT_MAX) continue;
          auto& mcp = (*mcps)[match.mcpi];
          myprt << "TrackId " << mcp.TrackId();
          myprt << ", PDG code " << mcp.PdgCode();
          int TMeV = 1000 * (mcp.E() - mcp.Mass());
          myprt << ", T " << TMeV << " MeV";
          myprt << ", Process " << mcp.Process();
          myprt << ", in TPC " << match.tpc;
          myprt << "\n";
          for(unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
            if (match.firstHit[plane] >= (*allHits).size()) continue;
            if (match.mcpTruHitCount[plane] < 3) continue;
            unsigned int first, last;
            FindFirstLastWire(allHits, hitMCPIndex, match, plane, first, last);
            if (first >= (*allHits).size()) {
              myprt<<" oops "<<match.mcpi<<" "<<plane<<" first "<<match.firstHit[plane]<<" "<<match.lastHit[plane]<<"\n";
              continue;
            }
            myprt << "    Plane " << plane;
            auto& fhit = (*allHits)[first];
            myprt << " true hits range " << PrintHit(fhit);
            auto& lhit = (*allHits)[last];
            myprt << " - " << PrintHit(lhit);
            myprt << " mcpTruHitCount " << match.mcpTruHitCount[plane];
            if (match.clsTruHitCount[plane] > 0) {
              auto& cls = (*allCls)[match.clsIndex[plane]];
              myprt << " -> Cls " <<cls.ID();
              auto clsHits = fmch.at(match.clsIndex[plane]);
              myprt << " nRecoHits " << clsHits.size();
              myprt << " nTruRecoHits " << match.clsTruHitCount[plane] << "\n";
            } // match.clsTruHitCount[plane] > 0
            else {
              myprt<<" *** No Cluster match\n";
            }
          } // plane
        } // match
      } else {
        myprt << "--> truncated print of " << matches.size() << " selected MCParticles";
      }
    } // fPrintLevel > 1

    // Calculate Efficiency and Purity
    ++fEventsProcessed;
    for(auto& match : matches) {
      if (match.mcpi == UINT_MAX) continue;
      auto& mcp = (*mcps)[match.mcpi];
      float TMeV = 1000 * (mcp.E() - mcp.Mass());
      int indx = PDGCodeIndex(mcp.PdgCode());
      if (indx < 0 || indx > 4) continue;
      for(unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
        if (match.mcpTruHitCount[plane] < 3) continue;
        float eff = 0;
        float pur = 0;
        if (match.clsIndex[plane] != UINT_MAX) {
          auto clsHits = fmch.at(match.clsIndex[plane]);
          if (clsHits.size() > 0) {
            eff = match.clsTruHitCount[plane] / match.mcpTruHitCount[plane];
            pur = match.clsTruHitCount[plane] / (float)clsHits.size();
          }
        }
        float effpur = eff * pur;
        // accumulate
        fEffSum[indx] += eff;
        fPurSum[indx] += pur;
        fEffPurSum[indx] += effpur;
        ++fSum[indx];
        // fill the histograms
        if (indx == 0) {
          fT_elec->Fill(TMeV);
          fE_elec->Fill(eff);
          fP_elec->Fill(pur);
          fEP_elec->Fill(effpur);
          fEP_T_elec->Fill(TMeV, effpur);
        } else if (indx == 1) {
          fT_muon->Fill(TMeV);
          fE_muon->Fill(eff);
          fP_muon->Fill(pur);
          fEP_muon->Fill(effpur);
          fEP_T_muon->Fill(TMeV, effpur);
        } else if (indx == 2) {
          fT_pion->Fill(TMeV);
          fE_pion->Fill(eff);
          fP_pion->Fill(pur);
          fEP_pion->Fill(effpur);
          fEP_T_pion->Fill(TMeV, effpur);
        } else if (indx == 3) {
          fT_kaon->Fill(TMeV);
          fE_kaon->Fill(eff);
          fP_kaon->Fill(pur);
          fEP_kaon->Fill(effpur);
          fEP_T_kaon->Fill(TMeV, effpur);
        } else if (indx == 4) {
          fT_prot->Fill(TMeV);
          fE_prot->Fill(eff);
          fP_prot->Fill(pur);
          fEP_prot->Fill(effpur);
          fEP_T_prot->Fill(TMeV, effpur);
        }
        if (effpur < 0.7) ++fNBadEP;
        if (fPrintLevel > 0 && effpur < 0.7) {
          mf::LogVerbatim myprt("ClusterAna");
          myprt << "BadEP evt " << evt.event();
          myprt << " " << fNames[indx];
          if (match.clsIndex[plane] != UINT_MAX) {
            auto& cls = (*allCls)[match.clsIndex[plane]];
            myprt << " clsID " << cls.ID();
          }
          myprt << " EP " <<std::fixed << std::setprecision(2) << effpur;
          unsigned int first, last;
          FindFirstLastWire(allHits, hitMCPIndex, match, plane, first, last);
          if (first != UINT_MAX) {
            auto& fhit = (*allHits)[first];
            auto& lhit = (*allHits)[last];
            myprt << " true hit range " << PrintHit(fhit) << " - " << PrintHit(lhit);
            myprt << " nTrueHits " << (int)match.mcpTruHitCount[plane];
          }
        } // fPrintLevel > 0 && effpur < 0.7
      } // plane
    } // match

  } // analyze

  void
  ClusterAna::FindFirstLastWire(art::Handle<std::vector<recob::Hit>> const& allHits,
                                std::vector<unsigned int> const& hitMCPIndex,
                                MatchStruct const& match,
                                unsigned short inPlane, unsigned int& first, unsigned int& last)
  {
    // Returns the index of the first hit and last hit in the specified plane
    // where first (last) = the lowest (highest) wire number
    first = UINT_MAX;
    last = 0;
    if (match.mcpi == UINT_MAX) return;
    unsigned int fWire = UINT_MAX;
    unsigned int lWire = 0;
    for(unsigned int iht = match.firstHit[inPlane]; iht <= match.lastHit[inPlane]; ++iht) {
      if (hitMCPIndex[iht] != match.mcpi) continue;
      auto& hit = (*allHits)[iht];
      if (hit.WireID().Plane != inPlane) continue;
      if (hit.WireID().TPC != match.tpc) continue;
      unsigned int wire = hit.WireID().Wire;
      if (wire < fWire) {
        fWire = wire;
        first = iht;
      } // wire > fWire
      if (wire > lWire) {
        lWire = wire;
        last = iht;
      } // wire > fWire
    } // iht
  } // FindFirstLastWire

  int
  ClusterAna::PDGCodeIndex(int pdgCode)
  {
    pdgCode = abs(pdgCode);
    if (pdgCode == 11) return 0;
    if (pdgCode == 13) return 1;
    if (pdgCode == 211) return 2;
    if (pdgCode == 321) return 3;
    if (pdgCode == 2212) return 4;
    return -1;
  }

  std::string
  ClusterAna::PrintHit(const recob::Hit& hit)
  {
    return std::to_string(hit.WireID().TPC) + ":" + std::to_string(hit.WireID().Plane) + ":" 
         + std::to_string(hit.WireID().Wire) + ":" + std::to_string((int)hit.PeakTime());
  } // PrintHit

} // end namespace

DEFINE_ART_MODULE(cluster::ClusterAna)
