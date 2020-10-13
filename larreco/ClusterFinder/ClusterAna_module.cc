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

  class ClusterAna : public art::EDAnalyzer {
  public:
    explicit ClusterAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void beginJob() override;
    void endJob() override;
    std::string PrintHit(const recob::Hit& hit);
    int EncodePDGCode(int pdgCode);

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

    art::InputTag fHitsModuleLabel;
    art::InputTag fClusterModuleLabel;
    art::InputTag fVertexModuleLabel;
    std::vector<float> fElecKERange;
    std::vector<float> fMuonKERange;
    std::vector<float> fPionKERange;
    std::vector<float> fKaonKERange;
    std::vector<float> fProtKERange;
    short fTrackWeightOption;
    bool fMergeDaughters;
    bool fSkipCosmics;
    bool fSkipMultiTPC;
    short fPrintLevel;
    short moduleID;

    std::array<std::string, 5> fNames {{"elec", "muon", "pion", "kaon", "prot"}};
    std::array<float, 5> fEffSum {{0}};
    std::array<float, 5> fPurSum {{0}};
    std::array<float, 5> fEffPurSum {{0}};
    std::array<float, 5> fSum {{0}};
    unsigned int fNBadEP {0};

  }; // class ClusterAna

}

namespace cluster {

  //--------------------------------------------------------------------
  ClusterAna::ClusterAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fHitsModuleLabel(pset.get<art::InputTag>("HitsModuleLabel"))
    , fClusterModuleLabel(pset.get<art::InputTag>("ClusterModuleLabel"))
    , fVertexModuleLabel(pset.get<art::InputTag>("VertexModuleLabel"))
    , fElecKERange(pset.get<std::vector<float>>("ElecKERange"))
    , fMuonKERange(pset.get<std::vector<float>>("MuonKERange"))
    , fPionKERange(pset.get<std::vector<float>>("PionKERange"))
    , fKaonKERange(pset.get<std::vector<float>>("KaonKERange"))
    , fProtKERange(pset.get<std::vector<float>>("ProtKERange"))
    , fSkipCosmics(pset.get<bool>("SkipCosmics", true))
    , fSkipMultiTPC(pset.get<bool>("SkipMultiTPC", true))
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
      if(fSkipCosmics) myprt << ": ignoring Cosmics";
      if(fSkipMultiTPC) myprt << ": ignoring MCParticles spanning multiple TPCs";
      myprt << "\n";
      myprt << "Hit format is <TPC>:<Plane>:<Wire>:<PeakTime>";
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
    myprt<<"ClusterAna results\n";
    float totSum = 0;
    float totEPSum = 0;
    for(unsigned short indx = 0; indx < 4; ++indx) {
      if(fSum[indx] == 0) continue;
        myprt << fNames[indx];
        float aveEP = fEffPurSum[indx] / fSum[indx];
        myprt << "EP " <<std::fixed << std::setprecision(2) << aveEP;
        myprt << " cnt " << (int)fSum[indx] << " ";
        totSum += fSum[indx];
        totEPSum += fEffPurSum[indx];
    } // indx
    if(totSum > 0) {
      totEPSum /= totSum;
      myprt << "AllEP " << totEPSum;
      myprt << " nBadEP " << fNBadEP;
    }
  } // endJob

  void
  ClusterAna::analyze(const art::Event& evt)
  {

    if(evt.isRealData()) return;

    // get a reference to the hit collection
    auto allHits = art::Handle<std::vector<recob::Hit>>();
    if(!evt.getByLabel(fHitsModuleLabel, allHits)) 
      throw cet::exception("ClusterAna")<<"Failed to get a handle to hit collection '"
      <<fHitsModuleLabel.label()<<"'\n";
    if((*allHits).empty()) return;
    // define a Hit -> MCParticle index assn
    std::vector<unsigned int> mcpIndex((*allHits).size(), UINT_MAX);

    // get true particles
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    art::ServiceHandle<geo::Geometry const> geom;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // MCParticle -> cluster matching struct and indices
    // of the first and last hit in ANY plane in the hit collection
    struct matchStruct {
      unsigned int mcpi {UINT_MAX};
      int mcpTrackId {INT_MAX};
      int tpc {-1};  // -1 = undefined, -2 = don't care
      std::vector<unsigned int> clsIndex {{UINT_MAX}};
      std::vector<float> clsTruHitCount {{0}};
      std::vector<float> mcpTruHitCount {{0}};
      unsigned int firstHit {UINT_MAX};
      unsigned int lastHit {UINT_MAX};
      unsigned short plane;
    };
    std::vector<matchStruct> matches;

    // get a handle to the MCParticle collection
    art::InputTag mcpLabel = "largeant";
    auto mcps = art::Handle<std::vector<simb::MCParticle>>();
    if(!evt.getByLabel(mcpLabel, mcps)) {
      std::cout<<"Failed to get a handle to largeant MCParticle collection\n";
      return;
    }

    for(unsigned int mcpi = 0; mcpi < (*mcps).size(); ++mcpi) {
      auto& mcp = (*mcps)[mcpi];
      int tid = mcp.TrackId();
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(tid);
      if(fSkipCosmics && theTruth->Origin() == simb::kCosmicRay) continue;
      int pdg = abs(mcp.PdgCode());
      float TMeV = 1000 * (mcp.E() - mcp.Mass());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      bool useIt = true;
      if (pdg == 11) {
        if (fElecKERange[0] < 0) useIt = false;
        // only consider primary electrons
        if (mcp.Process() != "primary") useIt = false;
        if (TMeV < fElecKERange[0] || TMeV > fElecKERange[1]) useIt = false;
      } else if (pdg == 13) {
        if (fMuonKERange[0] < 0) useIt = false;
        if (TMeV < fMuonKERange[0] || TMeV > fMuonKERange[1]) useIt = false;
      } else if (pdg == 211) {
        if (fPionKERange[0] < 0) useIt = false;
        if (TMeV < fPionKERange[0] || TMeV > fPionKERange[1]) useIt = false;
      } else if (pdg == 321) {
        if (fKaonKERange[0] < 0) useIt = false;
        if (TMeV < fKaonKERange[0] || TMeV > fKaonKERange[1]) useIt = false;
      } else if (pdg == 2212) {
        if (fProtKERange[0] < 0) useIt = false;
        if (TMeV < fProtKERange[0] || TMeV > fProtKERange[1]) useIt = false;
      }
      if(useIt) {
        matchStruct aMatch;
        aMatch.mcpi = mcpi;
        aMatch.mcpTrackId = tid;
        aMatch.clsIndex.resize(geom->Nplanes());
        aMatch.clsTruHitCount.resize(geom->Nplanes());
        aMatch.mcpTruHitCount.resize(geom->Nplanes());
        matches.push_back(aMatch);        
      } // useIT
    } // mcpi
    if(matches.empty()) return;

    // next match the hits. Populate mcpIndex and matches hit ranges
    for(unsigned int iht = 0; iht < (*allHits).size(); ++iht) {
      auto& hit = (*allHits)[iht];
      auto tides = bt_serv->HitToTrackIDEs(clockData, hit);
      for(auto& tide : tides) {
        if(tide.energyFrac > 0.5) {
          bool gotit = false;
          for(auto& match : matches) {
            if(match.mcpTrackId != tide.trackID) continue;
            mcpIndex[iht] = match.mcpi;
            if(match.firstHit == UINT_MAX) match.firstHit = iht;
            match.lastHit = iht;
            unsigned int plane = hit.WireID().Plane;
            ++match.mcpTruHitCount[plane];
            gotit = true;
          } // match
          if(gotit) break;          
        } // tide.energyFrac > 0.5
      } // tide
    } // iht

    // ignore MCParticles that don't have the first/last hit defined or 
    // that span several TPCs if the user has selected that option. This is done
    // by setting the match tpc variable to -1 -> ignore this MCParticle
    for(auto& match : matches) {
      if(match.mcpi == UINT_MAX) continue;
      if(match.firstHit == UINT_MAX) continue;
      auto& fhit = (*allHits)[match.firstHit];
      match.tpc = fhit.WireID().TPC;
      if(fSkipMultiTPC) {
        auto& lhit = (*allHits)[match.lastHit];
        if((int)lhit.WireID().TPC != match.tpc) {
          if(fPrintLevel > 2) mf::LogVerbatim("ClusterAna")<<"mcpi "<<match.mcpi<<" isn't Good ";
          // update the mcpIndex vector
          for(auto& mcpi : mcpIndex) if(mcpi == match.mcpi) mcpi = UINT_MAX;
          match.mcpi = UINT_MAX;
        } // lhit.WireID().TPC != match.tpc
      } // fSkipMultiTPC
    } // match

    // get clusters and cluster-hit associations
    art::Handle<std::vector<recob::Cluster>> allCls;
    evt.getByLabel(fClusterModuleLabel, allCls);
    art::FindManyP<recob::Hit> fmh(allCls, evt, fClusterModuleLabel);
    if (allCls->size() == 0) {
      std::cout << "No clusters in event. Hits size " << (*allHits).size() << "\n";
      // don't return or the statistics will be wrong...
    }

    // Match MCParticles -> Clusters
    for (unsigned int icl = 0; icl < allCls->size(); ++icl) {
      // Find the best match of this cluster to a MCParticle.
      // Temporary vector of (mcp index, match count) pairs for this cluster
      std::vector<std::pair<unsigned int, float>> mcpCnts;
      auto clsHits = fmh.at(icl);
      unsigned short plane = USHRT_MAX;
      for(auto& clsHit : clsHits) {
        unsigned int iht = clsHit.key();
        if(plane == USHRT_MAX) plane = (*allHits)[iht].WireID().Plane;
        // ignore matches to MCParticles that aren't relevant
        if(mcpIndex[iht] <= 0) continue;
        // find the index of an existing Cluster -> MCParticle in the temporary vector
        unsigned short indx = 0;
        for(indx = 0; indx < mcpCnts.size(); ++ indx) if(mcpCnts[indx].first == mcpIndex[iht]) break;
        // no match found so add one
        if(indx == mcpCnts.size()) mcpCnts.push_back(std::make_pair(mcpIndex[iht], 0));
        ++mcpCnts[indx].second;
      } // clsHit
      if(mcpCnts.empty()) continue;
      // Find the best Cluster -> MCParticle match
      float bestHitsMatched = 0;
      unsigned int bestMCPIndex = 0;
      for(auto& mcpCnt : mcpCnts) {
        if(mcpCnt.second < bestHitsMatched) continue;
        bestHitsMatched = mcpCnt.second;
        bestMCPIndex = mcpCnt.first;
      } // mcpCnt
      // compare the number of matched hits for this cluster with an existing
      // MCParticle -> cluster match
      unsigned int mcpi = bestMCPIndex;
      for(auto& match : matches) {
        if(match.mcpi != mcpi) continue;
        // Found an existing match. Are there more matched hits?
        if(match.clsTruHitCount[plane] > bestHitsMatched) continue;
        // a better match
        match.clsTruHitCount[plane] = bestHitsMatched;
        match.clsIndex[plane] = icl;
        break;
      } // match
    } // icl
    if (fPrintLevel > 2) {
      mf::LogVerbatim myprt("ClusterAna");
      for(auto& match : matches) {
        if(match.mcpi == UINT_MAX) continue;
        auto& mcp = (*mcps)[match.mcpi];
        myprt << "TrackId " << mcp.TrackId();
        myprt << " PDG code " << mcp.PdgCode();
        int TMeV = 1000 * (mcp.E() - mcp.Mass());
        myprt << " T " << TMeV << " MeV,";
        myprt << " Process " << mcp.Process();
        myprt << " in TPC " << match.tpc;
        myprt << "\n";
        if(match.firstHit < (*allHits).size()) {
          for(unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
            if(match.mcpTruHitCount[plane] < 3) continue;
            unsigned int first = UINT_MAX;
            unsigned int last = UINT_MAX;
            for(unsigned int iht = match.firstHit; iht <= match.lastHit; ++iht) {
              if(mcpIndex[iht] != match.mcpi) continue;
              auto& hit = (*allHits)[iht];
              if(hit.WireID().Plane != plane) continue;
              if(first == UINT_MAX) first = iht;
              last = iht;
            } // iht
            myprt << "    Plane " << plane;
            auto& fhit = (*allHits)[first];
            myprt << " true hits range " << PrintHit(fhit);
            auto& lhit = (*allHits)[last];
            myprt << " - " << PrintHit(lhit);
            myprt << " mcpTruHitCount " << match.mcpTruHitCount[plane];
            if(match.clsTruHitCount[plane] > 0) {
              auto& cls = (*allCls)[match.clsIndex[plane]];
              myprt << " Cls " <<cls.ID();
              auto clsHits = fmh.at(match.clsIndex[plane]);
              myprt << " nRecoHits " << clsHits.size();
              myprt << " nTruRecoHits " << match.clsTruHitCount[plane] << "\n";
            } // match.clsTruHitCount[plane] > 0
            else {
              myprt<<" *** No MCParticle -> Cluster match\n";
            }
          } // plane
        } // valid firstHit match
      } // match
    } // fPrintLevel > 2

    // Calculate Efficiency and Purity
    for(auto& match : matches) {
      if(match.mcpi == UINT_MAX) continue;
      auto& mcp = (*mcps)[match.mcpi];
      float TMeV = 1000 * (mcp.E() - mcp.Mass());
      int indx = EncodePDGCode(mcp.PdgCode());
      if(indx < 0 || indx > 4) continue;
      for(unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
        if(match.mcpTruHitCount[plane] < 3) continue;
        float eff = 0;
        float pur = 0;
        if(match.clsIndex[plane] != UINT_MAX) {
          auto clsHits = fmh.at(match.clsIndex[plane]);
          if(clsHits.size() > 0) {
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
        if(indx == 0) {
          fT_elec->Fill(TMeV);
          fE_elec->Fill(eff);
          fP_elec->Fill(pur);
          fEP_elec->Fill(effpur);
          fEP_T_elec->Fill(TMeV, effpur);
        } else if(indx == 1) {
          fT_muon->Fill(TMeV);
          fE_muon->Fill(eff);
          fP_muon->Fill(pur);
          fEP_muon->Fill(effpur);
          fEP_T_muon->Fill(TMeV, effpur);
        } else if(indx == 2) {
          fT_pion->Fill(TMeV);
          fE_pion->Fill(eff);
          fP_pion->Fill(pur);
          fEP_pion->Fill(effpur);
          fEP_T_pion->Fill(TMeV, effpur);
        } else if(indx == 3) {
          fT_kaon->Fill(TMeV);
          fE_kaon->Fill(eff);
          fP_kaon->Fill(pur);
          fEP_kaon->Fill(effpur);
          fEP_T_kaon->Fill(TMeV, effpur);
        } else if(indx == 4) {
          fT_prot->Fill(TMeV);
          fE_prot->Fill(eff);
          fP_prot->Fill(pur);
          fEP_prot->Fill(effpur);
          fEP_T_prot->Fill(TMeV, effpur);
        }
        if(effpur < 0.7) ++fNBadEP;
        if(fPrintLevel > 0 && effpur < 0.7) {
          mf::LogVerbatim myprt("ClusterAna");
          myprt << "BadEP evt " << evt.event();
          myprt << " " << fNames[indx];
          myprt << " EP " <<std::fixed << std::setprecision(2) << effpur;
          unsigned int first = UINT_MAX;
          unsigned int last = UINT_MAX;
          for(unsigned int iht = match.firstHit; iht <= match.lastHit; ++iht) {
            auto& hit = (*allHits)[iht];
            if(hit.WireID().Plane != plane) continue;
            if(first == UINT_MAX) first = iht;
            last = iht;
          } // iht
          if(first != UINT_MAX) {
            auto& fhit = (*allHits)[first];
            auto& lhit = (*allHits)[last];
            myprt << " true hit range " << PrintHit(fhit) << " - " << PrintHit(lhit);
            myprt << " nTrueHits " << (int)match.mcpTruHitCount[plane];
          }
        } // fPrintLevel > 0 && effpur < 0.7
      } // plane
    } // match

  } // analyze

  int
  ClusterAna::EncodePDGCode(int pdgCode)
  {
    pdgCode = abs(pdgCode);
    if(pdgCode == 11) return 0;
    if(pdgCode == 13) return 1;
    if(pdgCode == 211) return 2;
    if(pdgCode == 321) return 3;
    if(pdgCode == 2212) return 4;
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
