#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"


#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

namespace tca {

  //////////////////////////////////////////
  void TruthMatcher::Initialize()
  {
     // Initialize the variables used to calculate Efficiency * Purity (aka EP) for matching to truth
    EPCnts.fill(0); 
    TSums.fill(0.0);
    EPTSums.fill(0.0);
    TruVxCounts.fill(0);
    nBadEP = 0;
  } // Initialize
/*
  //////////////////////////////////////////
  void TruthMatcher::MatchTrueHits(TCSlice& slc)
  {
    // Matches reco hits to MC true tracks and puts the association into
    // TCHit TruTrkID. This code is almost identical to the first part of MatchTruth.
    

    if(tcc.matchTruth[0] < 0) return;
    if(slc.slHits.empty()) return;
    // ensure that the hit vector and mcpList vector sizes are the same
    slc.mcpListIndex.resize(slc.slHits.size());
    
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // list of all true particles
    sim::ParticleList const& plist = pi_serv->ParticleList();
    if(plist.empty()) return;
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;

    // first find the source particle 
    simb::Origin_t sourceOrigin = simb::kUnknown;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* mcp = (*ipart).second;
      int trackID = mcp->TrackId();
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
      if(tcc.matchTruth[0] == 1) {
        // Look for beam neutrino or single particle
        if(theTruth->Origin() == simb::kBeamNeutrino) {
          sourcePtclTrackID = trackID;
          sourceOrigin = simb::kBeamNeutrino;
          if(tcc.matchTruth[1] > 0) {
            Vector3_t dir {{mcp->Px(), mcp->Py(), mcp->Pz()}};
            SetMag(dir, 1);
            std::cout<<"Found beam neutrino sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode();
            std::cout<<" Vx "<<std::fixed<<std::setprecision(1)<<mcp->Vx()<<" Vy "<<mcp->Vy()<<" Vz "<<mcp->Vz();
            std::cout<<" dir "<<std::setprecision(3)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
          }
          break;
        } // beam neutrino
        if(theTruth->Origin() == simb::kSingleParticle) {
          sourcePtclTrackID = trackID;
          sourceOrigin = simb::kSingleParticle;
          if(tcc.matchTruth[1] > 0) {
            Vector3_t dir {{mcp->Px(), mcp->Py(), mcp->Pz()}};
            SetMag(dir, 1);
            std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode();
            std::cout<<" Vx "<<std::fixed<<std::setprecision(1)<<mcp->Vx()<<" Vy "<<mcp->Vy()<<" Vz "<<mcp->Vz();
            std::cout<<" dir "<<std::setprecision(3)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            break;
          }
        } // single particle
      } else if(tcc.matchTruth[0] == 2 && theTruth->Origin() == simb::kCosmicRay) {
        sourcePtclTrackID = trackID;
        sourceOrigin = simb::kCosmicRay;
      }
    } // ipart
    
    if(sourcePtclTrackID == -1) {
      if(tcc.matchTruth[1] > 0) std::cout<<"MatchTrueHits: SourcePtcl not found\n";
      return;
    }
    
    if(tcc.matchTruth[1] > 2) {
      // print out a whole bunch of information
      mf::LogVerbatim myprt("TC");
      myprt<<"Displaying all neutrino origin MCParticles with T > 50 MeV\n";
      myprt<<" trackID PDGCode Mother T(MeV) ________dir_______            Process\n";
      for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
        simb::MCParticle* mcp = (*ipart).second;
        int trackID = mcp->TrackId();
        art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
        if(theTruth->Origin() != simb::kBeamNeutrino) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (mcp->E() - mcp->Mass());
        if(TMeV < 50) continue;
        myprt<<std::setw(8)<<trackID;
        myprt<<std::setw(8)<<mcp->PdgCode();
        myprt<<std::setw(8)<<mcp->Mother();
        myprt<<std::setw(6)<<TMeV;
        Point3_t pos;
        pos[0] = mcp->Vx();
        pos[1] = mcp->Vy();
        pos[2] = mcp->Vz();
        Vector3_t dir {{mcp->Px(), mcp->Py(), mcp->Pz()}};
        SetMag(dir, 1);
        myprt<<std::setprecision(2);
        for(unsigned short xyz = 0; xyz < 3; ++xyz) myprt<<std::setw(6)<<dir[xyz];
        myprt<<std::setw(22)<<mcp->Process();
        myprt<<"\n";
      } // ipart
    } // big print
      
    // flag MCParticles to select for measuring performance. 
    std::vector<bool> select(plist.size(), false);
    // ensure that we get all of the primary particles
    unsigned int indx = 0;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* mcp = (*ipart).second;
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(mcp->TrackId());
      if(theTruth->Origin() != sourceOrigin) continue;
      select[indx] = true;
      ++indx;
    }

    // make a temp vector of hit -> geant trackID
    std::vector<int> gtid(slc.slHits.size(), 0);
    // find hits that match to the source particle
    for(unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      std::vector<sim::TrackIDE> ides;
      // HitToTrackIDEs accepts either a recob::Hit or art::Ptr<recob::Hit> so we need to
      // reference the hit. 
      recob::Hit hit = evt.allHits[slc.slHits[iht].allHitsIndex];
      ides = bt_serv->HitToTrackIDEs(hit);
      if(ides.empty()) continue;
      float energy = 0;
      for(auto ide : ides) energy += ide.energy;
      if(energy == 0) continue;
      // require 1/2 of the energy be due to one MC particle
      energy /= 2;
      int hitTruTrkID = 0;
      for(auto ide : ides) {
        if(ide.energy > energy) {
          hitTruTrkID = ide.trackID;
          break;
        }
      } // ide
      if(hitTruTrkID == 0) continue;
      // ensure it has the correct source
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(hitTruTrkID);
      if(theTruth->Origin() != sourceOrigin) continue;
      unsigned int indx = 0;
      for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
        ++indx;
        simb::MCParticle* mcp = (*ipart).second;
        if(mcp->TrackId() != hitTruTrkID) continue;
        select[indx - 1] = true;
        gtid[iht] = mcp->TrackId();
        // find the eve particle and select it as well
        const simb::MCParticle* momMCP = pi_serv->TrackIdToMotherParticle_P(mcp->TrackId());
        if(!momMCP) continue;
        if(momMCP->TrackId() == mcp->TrackId()) break;
        unsigned int mindx = 0;
        for(sim::ParticleList::const_iterator mpart = plist.begin(); mpart != plist.end(); ++mpart) {
          simb::MCParticle* mmcp = (*mpart).second;
          if(mmcp->TrackId() == momMCP->TrackId()) {
            select[mindx] = true;
            break;
          }
          ++mindx;
        } // mpart
        break;
      } // ipart
    } // iht

    // save the selected MCParticles
    indx = 0;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      if(select[indx]) {
        simb::MCParticle* mcp = (*ipart).second;
        slc.mcpList.push_back(mcp);
      }
      ++indx;
    } // ipart
    
    if(slc.mcpList.size() > USHRT_MAX) {
      std::cout<<"MatchTrueHits: mcpList size too large\n";
      slc.mcpList.clear();
      return;
    } // mcpList too large
    
    // define mcpListIndex for the hits
    unsigned short nMatch = 0;
    for(unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      if(gtid[iht] == 0) continue;
      for(unsigned short indx = 0; indx < slc.mcpList.size(); ++indx) {
        auto& mcp = slc.mcpList[indx];
        if(mcp->TrackId() != gtid[iht]) continue;
        slc.mcpListIndex[iht] = indx;
        ++nMatch;
        break;
      } // indx
    } // iht
    
    if(tcc.matchTruth[1] > 0) std::cout<<"MatchTrueHits: mcpList size "<<slc.mcpList.size()<<" slHits size "<<slc.slHits.size()<<" num MC-matched hits "<<nMatch<<"\n";

  } // MatchTrueHits
  
  //////////////////////////////////////////
  void TruthMatcher::StudyShowerParents(TCSlice& slc, HistStuff& hist)
  {
    // study characteristics of shower parent pfps. This code is adapted from TCShower FindParent
    if(slc.pfps.empty()) return;
    if(slc.mcpList.empty()) return;
    
    // Look for truth pfp primary electron
    Point3_t primVx {{-666.0, -666.0, -666.0}};
    // the primary should be the first one in the list as selected in GetHitCollection
    auto& primMCP = slc.mcpList[0];
    primVx[0] = primMCP->Vx();
    primVx[1] = primMCP->Vy();
    primVx[2] = primMCP->Vz();
    geo::Vector_t posOffsets;
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    posOffsets = SCE->GetPosOffsets({primVx[0], primVx[1], primVx[2]});
    posOffsets.SetX(-posOffsets.X());
    primVx[0] += posOffsets.X();
    primVx[1] += posOffsets.Y();
    primVx[2] += posOffsets.Z();
    geo::TPCID inTPCID;
    // ignore if the primary isn't inside a TPC
    if(!InsideTPC(primVx, inTPCID)) return;
    // or if it is inside the wrong tpc
    if(inTPCID != slc.TPCID) return;

    std::string fcnLabel = "SSP";
    // Create a truth shower for each primary electron
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    MCParticleListUtils mcpu{slc};
    for(unsigned int part = 0; part < slc.mcpList.size(); ++part) {
      auto& mcp = slc.mcpList[part];
      // require electron or photon
      if(abs(mcp->PdgCode()) != 11 && abs(mcp->PdgCode()) != 111) continue;
      int eveID = pi_serv->ParticleList().EveId(mcp->TrackId());
      // require that it is primary
      if(mcp->TrackId() != eveID) continue;
      int truPFP = 0;
      auto ss3 = mcpu.MakeCheatShower(slc, part, primVx, truPFP);
      if(ss3.ID == 0) continue;
      if(truPFP == 0) continue;
      if(!StoreShower(fcnLabel, slc, ss3)) {
        std::cout<<"Failed to store 3S"<<ss3.ID<<"\n";
        break;
      } // store failed
      // now fill the TTree
      float ss3Energy = ShowerEnergy(ss3);
      for(auto& pfp : slc.pfps) {
        if(pfp.TPCID != ss3.TPCID) continue;
        // ignore neutrinos
        if(pfp.PDGCode == 12 || pfp.PDGCode == 14) continue;
        // ignore shower pfps
        if(pfp.PDGCode == 1111) continue;
        float pfpEnergy = 0;
        float minEnergy = 1E6;
        for(auto tid : pfp.TjIDs) {
          auto& tj = slc.tjs[tid - 1];
          float energy = ChgToMeV(tj.TotChg);
          pfpEnergy += energy;
          if(energy < minEnergy) minEnergy = energy;
        }
        pfpEnergy -= minEnergy;
        pfpEnergy /= (float)(pfp.TjIDs.size() - 1);
        // find the end that is farthest away
        unsigned short pEnd = FarEnd(slc, pfp, ss3.ChgPos);
        auto pToS = PointDirection(pfp.XYZ[pEnd], ss3.ChgPos);
        // take the absolute value in case the shower direction isn't well known
        float costh1 = std::abs(DotProd(pToS, ss3.Dir));
        float costh2 = DotProd(pToS, pfp.Dir[pEnd]);
        // distance^2 between the pfp end and the shower start, charge center, and shower end
        float distToStart2 = PosSep2(pfp.XYZ[pEnd], ss3.Start);
        float distToChgPos2 = PosSep2(pfp.XYZ[pEnd], ss3.ChgPos);
        float distToEnd2 = PosSep2(pfp.XYZ[pEnd], ss3.End);
//        mf::LogVerbatim("TC")<<" 3S"<<ss3.ID<<" P"<<pfp.ID<<"_"<<pEnd<<" distToStart "<<sqrt(distToStart2)<<" distToChgPos "<<sqrt(distToChgPos2)<<" distToEnd "<<sqrt(distToEnd2);
        // find the end of the shower closest to the pfp
        unsigned short shEnd = 0;
        if(distToEnd2 < distToStart2) shEnd = 1;
        if(shEnd == 0 && distToChgPos2 < distToStart2) continue;
        if(shEnd == 1 && distToChgPos2 < distToEnd2) continue;
//      mf::LogVerbatim("TC")<<" 3S"<<ss3.ID<<"_"<<shEnd<<" P"<<pfp.ID<<"_"<<pEnd<<" costh1 "<<costh1;
        Point2_t alongTrans;
        // find the longitudinal and transverse components of the pfp start point relative to the
        // shower center
        FindAlongTrans(ss3.ChgPos, ss3.Dir, pfp.XYZ[pEnd], alongTrans);
//      mf::LogVerbatim("TC")<<"   alongTrans "<<alongTrans[0]<<" "<<alongTrans[1];
        hist.fSep = sqrt(distToChgPos2);
        hist.fShEnergy = ss3Energy;
        hist.fPfpEnergy = pfpEnergy;
        hist.fPfpLen = PosSep(pfp.XYZ[0], pfp.XYZ[1]);
        hist.fMCSMom = MCSMom(slc, pfp.TjIDs);
        hist.fDang1 = acos(costh1);
        hist.fDang2 = acos(costh2);
        hist.fChgFrac = 0;
        float chgFrac = 0;
        float totSep = 0;
        // find the charge fraction btw the pfp start and the point that is 
        // half the distance to the charge center in each plane
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(ss3.TPCID.Cryostat, ss3.TPCID.TPC, plane);
          int ssid = 0;
          for(auto cid : ss3.CotIDs) {
            auto& ss = slc.cots[cid - 1];
            if(ss.CTP != inCTP) continue;
            ssid = ss.ID;
            break;
          } // cid
          if(ssid == 0) continue;
          auto tpFrom = MakeBareTP(slc, pfp.XYZ[pEnd], pToS, inCTP);
          auto& ss = slc.cots[ssid - 1];
          auto& stp1 = slc.tjs[ss.ShowerTjID - 1].Pts[1];
          float sep = PosSep(tpFrom.Pos, stp1.Pos);
          float toPos = tpFrom.Pos[0] + 0.5 * tpFrom.Dir[0] * sep;
          float cf = ChgFracBetween(slc, tpFrom, toPos, false);
          // weight by the separation in the plane
          totSep += sep;
          chgFrac += sep * cf;
        } // plane
        if(totSep > 0) hist.fChgFrac = chgFrac / totSep;
        hist.fAlong = alongTrans[0];
        hist.fTrans = alongTrans[1];
        hist.fInShwrProb = InShowerProbLong(ss3Energy, -hist.fSep);
        bool isBad = (hist.fDang1 > 2 || hist.fChgFrac < 0.5 || hist.fInShwrProb < 0.05);
        if(pfp.ID == truPFP && isBad) {
          mf::LogVerbatim myprt("TC");
          myprt<<"SSP: 3S"<<ss3.ID<<" shEnergy "<<(int)ss3Energy<<" P"<<pfp.ID<<" pfpEnergy "<<(int)pfpEnergy;
          myprt<<" MCSMom "<<hist.fMCSMom<<" len "<<hist.fPfpLen;
          myprt<<" Dang1 "<<hist.fDang1<<" Dang2 "<<hist.fDang2<<" chgFrac "<<hist.fChgFrac;
          myprt<<" fInShwrProb "<<hist.fInShwrProb;
          myprt<<" EventsProcessed "<<evt.eventsProcessed;
        }
        if(pfp.ID == truPFP) {
          hist.fShowerParentSig->Fill();
        } else {
          hist.fShowerParentBkg->Fill();
        }
      } // pfp
    } // part
//    PrintShowers(fcnLabel, tjs);
//    Print2DShowers(fcnLabel, slc, USHRT_MAX, false);
    // kill the cheat showers
    for(auto& ss3 : slc.showers) {
      if(ss3.ID == 0) continue;
      if(!ss3.Cheat) continue;
      for(auto cid : ss3.CotIDs) {
        auto& ss = slc.cots[cid - 1];
        ss.ID = 0;
        auto& stj = slc.tjs[ss.ShowerTjID - 1];
        stj.AlgMod[kKilled] = true;
      } // cid
      ss3.ID = 0;
    } // ss3
  } // StudyShowerParents
  
  //////////////////////////////////////////
  void TruthMatcher::StudyElectrons(TCSlice& slc, const HistStuff& hist)
  {
    // study tjs matched to electrons to develop an electron tag
    
//    float likely;
//    bool flipDirection;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.mcpListIndex == UINT_MAX) continue;
      unsigned short npts = tj.EndPt[1] - tj.EndPt[0] + 1;
      if(npts < 10) continue;
//      PrimaryElectronLikelihood(slc, tj, likely, flipDirection, true);
      auto& mcp = slc.mcpList[tj.mcpListIndex];
      unsigned short pdgIndex = PDGCodeIndex(slc, mcp->PdgCode());
      if(pdgIndex > 4) continue;
      unsigned short midpt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
      float mom1 = MCSMom(slc, tj, tj.EndPt[0], midpt);
      float mom2 = MCSMom(slc, tj, midpt, tj.EndPt[1]);
      float asym = std::abs(mom1 - mom2) / (mom1 + mom2);
      hist.fChgRMS[pdgIndex]->Fill(tj.ChgRMS);
      hist.fMomAsym[pdgIndex]->Fill(asym);
      float elike = asym * tj.ChgRMS;
      hist.fElectronLike[pdgIndex]->Fill(elike);
      float len = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
      hist.fElectronLike_Len[pdgIndex]->Fill(len, elike);
    } // tj
  } // StudyElectrons
  
  //////////////////////////////////////////
  void TruthMatcher::StudyPiZeros(TCSlice& slc, const HistStuff& hist)
  {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    sim::ParticleList const& plist = pi_serv->ParticleList();

    unsigned short cnt = 0;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      ++cnt;
      const simb::MCParticle* p = (*ipart).second;
      int pdg = abs(p->PdgCode());
      if(cnt == 1 && pdg != 111) {
        std::cout<<"First MC particle isn't a pizero\n";
        return;
      }
      if(cnt == 1) continue;
      if(pdg != 22) break;
      double photE = 1000 * p->E();
      if(photE < 50) continue;
//      Point3_t photStart {p->Vx(), p->Vy(), p->Vz()};
      Vector3_t photDir {{p->Px(), p->Py(), p->Pz()}};
      SetMag(photDir, 1);
      // sum up the charge for all hits that are daughters of this photon
      std::vector<float> chgSum(3);
      for(unsigned int iht = 0; iht < slc.mcpListIndex.size(); ++iht) {
        auto mcpIndex = slc.mcpListIndex[iht];
        if(mcpIndex ==  USHRT_MAX) continue;
        auto& mcp = slc.mcpList[mcpIndex];
        int eveID = pi_serv->ParticleList().EveId(mcp->TrackId());
        if(eveID != p->TrackId()) continue;
        auto& hit = evt.allHits[slc.slHits[iht].allHitsIndex];
        unsigned short plane = hit->WireID().Plane;
        chgSum[plane] += hit->Integral();
      } // iht
      for(unsigned short plane = 0; plane < 3; ++plane) {
        if(chgSum[plane] == 0) continue;
        float chgCal = photE / chgSum[plane];
        hist.fChgToMeV[plane]->Fill(chgCal);
        hist.fChgToMeV_Etru->Fill(photE, chgCal);
      } // plane
    } // ipart

  } // StudyPiZeros

  //////////////////////////////////////////
  void TruthMatcher::MatchTruth(TCSlice& slc, const HistStuff& hist, bool fStudyMode)
  {
    // The hits have already been matched to the truth in MatchTrueHits. Here we match reconstructed objects
    // to the truth-matched hits to measure performance    
    
    if(tcc.matchTruth[0] < 0) return;
    if(slc.mcpList.empty()) return;
    for(auto& pfp : slc.pfps) pfp.EffPur = 0;

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    geo::Vector_t posOffsets;
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
//    if(!SCE->EnableSimSpatialSCE()) std::cout<<">>>> EnableSimSpatialSCE isn't active...\n";

    // these are only applicable to neutrinos
    bool neutrinoVxReconstructable = false;
    bool vxReconstructedNearNuVx = false;
    bool neutrinoPFPCorrect = false;
    // Locate the primary vertex
    Point3_t primVx {{-666.0, -666.0, -666.0}};
    // the primary should be the first one in the list as selected in GetHitCollection
    auto& primMCP = slc.mcpList[0];
    primVx[0] = primMCP->Vx();
    primVx[1] = primMCP->Vy();
    primVx[2] = primMCP->Vz();
    posOffsets = SCE->GetPosOffsets({primVx[0], primVx[1], primVx[2]});
    posOffsets.SetX(-posOffsets.X());
    primVx[0] += posOffsets.X();
    primVx[1] += posOffsets.Y();
    primVx[2] += posOffsets.Z();
    if(tcc.matchTruth[1] > 1) std::cout<<"Prim PDG code  "<<primMCP->PdgCode()<<" primVx "<<std::fixed<<std::setprecision(1)<<primVx[0]<<" "<<primVx[1]<<" "<<primVx[2]<<"\n";
    geo::TPCID inTPCID = slc.TPCID;
    if(!InsideTPC(primVx, inTPCID)) {
      if(tcc.matchTruth[1] > 0) std::cout<<"Found a primary particle but it is not inside any TPC\n";
      return;
    }
    neutrinoVxReconstructable = true;
    
    // Look for the MC truth process that should be considered (beam neutrino,
    // single particle, cosmic rays), then form a list of selected MCParticles 
    // that will be used to measure performance. Feb 16: Changed GetHitCollection to only
    // store the MCParticles for the desired MCTruth collection
    std::vector<unsigned int> mcpSelect;
    // vector of reconstructable primary particles
    std::vector<unsigned int> primMCPs;
    for(unsigned int part = 0; part < slc.mcpList.size(); ++part) {
      auto& mcp = slc.mcpList[part];
      // require it is charged
      int pdg = abs(mcp->PdgCode());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      // require that it can be reconstructed in 3D
      if(!CanReconstruct(part, 3, inTPCID)) continue;
      mcpSelect.push_back(part);
      // Now require MCParticle primaries
      if(mcp->Mother() != 0) continue;
      // add to the list of primaries
      primMCPs.push_back(part);
    } // part
    if(mcpSelect.empty()) return;

    if(neutrinoVxReconstructable) ++TruVxCounts[0];
    
    // Form a list of mother-daughter pairs that should be considered as a single particle
    std::vector<std::pair<unsigned int, unsigned int>> moda;
    for(unsigned int part = 0; part < mcpSelect.size(); ++part) {
      auto& mcp = slc.mcpList[mcpSelect[part]];
      if(mcp->NumberDaughters() == 0) continue;
      unsigned int ndtr = 0;
      unsigned int dtrSelIndex = 0;
      for(int idtr = 0; idtr < mcp->NumberDaughters(); ++idtr) {
        int dtrTrackId = mcp->Daughter(idtr);
        // ignore it if it's not in the list
        bool ignore = true;
        for(unsigned short dpart = part + 1; dpart < mcpSelect.size(); ++dpart) {
          auto& dmcp = slc.mcpList[mcpSelect[dpart]];
          if(dmcp->TrackId() == dtrTrackId) {
            dtrSelIndex = dpart;
            ignore = false;
          }
        } // dpart
        if(ignore) continue;
        ++ndtr;
      } // idtr
      // require only one daughter
      if(ndtr != 1) continue;
      // require that the daughter have the same PDG code
      auto& dtr = slc.mcpList[mcpSelect[dtrSelIndex]];
      if(dtr->PdgCode() != mcp->PdgCode()) continue;
      if(tcc.matchTruth[1] > 1) mf::LogVerbatim("TC")<<"Daughter MCP "<<dtrSelIndex<<" -> Mother MCP "<<part;
      moda.push_back(std::make_pair(mcpSelect[part], mcpSelect[dtrSelIndex]));
      // delete the daughter entry
      mcpSelect.erase(mcpSelect.begin() + dtrSelIndex);
    } // part
    
    // NOTE: The PutMCPHitsInVector function will return an incomplete set of hits
    // matched to a MCParticle if it is called before mother-daughter assocations are corrected below.
    if(!moda.empty()) {
      // over-write the daughter hit -> MCParticle association with the mother.
      // Note that mother-daughter pairs are ordered by increasing generation. Reverse
      // moda so that the grand-daughters come first. Grand-daughters will then be
      // over-written by daughters in the previous generation
      if(moda.size() > 1) std::reverse(moda.begin(), moda.end());
      for(auto& mcpi : slc.mcpListIndex) {
        // see if this is a daughter
        for(auto& md : moda) if(md.second == mcpi) mcpi = md.first;
      } // hit
    } // !moda.empty
    
    // True histograms
    for(unsigned int part = 0; part < mcpSelect.size(); ++part) {
      auto& mcp = slc.mcpList[mcpSelect[part]];
      unsigned short pdgIndex = PDGCodeIndex(slc, mcp->PdgCode());
      if(pdgIndex > 4) continue;
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      hist.fTruT[pdgIndex]->Fill(TMeV);
    } // part

    // Match tjs to MC particles. Declare it a match to the MC particle that has the most
    // hits in the tj
    std::vector<std::vector<int>> tjid(slc.nPlanes);
    std::vector<std::vector<unsigned short>> ntht(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      tjid[plane].resize(slc.mcpList.size());
      ntht[plane].resize(slc.mcpList.size());
    }
    
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.Cryostat != inTPCID.Cryostat) continue;
      if(planeID.TPC != inTPCID.TPC) continue;
      tj.mcpListIndex = UINT_MAX;
      std::vector<unsigned int> mcpIndex, cnt;
      auto tjhits = PutTrajHitsInVector(tj, kUsedHits);
      for(auto iht : tjhits) {
        unsigned short mcpli = slc.mcpListIndex[iht];
        if(mcpli > slc.mcpList.size() - 1) continue;
        // ignore Tj hits that aren't in mcpSelect
        if(std::find(mcpSelect.begin(), mcpSelect.end(), mcpli) == mcpSelect.end()) continue;
        unsigned int indx = 0;
        for(indx = 0; indx < slc.mcpListIndex.size(); ++indx) if(mcpli == mcpListIndex[indx]) break;
        if(indx == slc.mcpListIndex.size()) {
          mcpIndex.push_back(mcpli);
          cnt.push_back(1);
        } else {
          ++cnt[indx];
        }
      } // iht
      unsigned short maxcnt = 0;
      unsigned int tmpIndex = UINT_MAX;
      for(unsigned short ii = 0; ii < cnt.size(); ++ii) {
        if(cnt[ii] > maxcnt) {
          maxcnt = cnt[ii];
          tmpIndex = mcpIndex[ii];
        }
      } // ii
      if(tmpIndex == UINT_MAX) continue;
      // calculate efficiency and purity.
      // Put the truth-matched hits into a vector for counting
      auto mcpHits = PutMCPHitsInVector(tmpIndex, tj.CTP);
      if(mcpHits.size() < 3) continue;
      float eff = (float)maxcnt / (float)mcpHits.size();
      float pur = (float)maxcnt / (float)tjhits.size();
      float ep = eff * pur;
      // see if there is a previous match with more truth-matched hits
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      if(tjid[plane][tmpIndex] > 0) {
        if(maxcnt > ntht[plane][tmpIndex]) {
          // clear the old one
          auto& mtj = slc.tjs[tjid[plane][tmpIndex] - 1];
          mtj.EffPur = 0;
          mtj.mcpListIndex = UINT_MAX;
          // use the new one
          tj.EffPur = ep;
          tj.mcpListIndex = tmpIndex;
          tjid[plane][tmpIndex] = tj.ID;
          ntht[plane][tmpIndex] = maxcnt;
        } else {
          // use the old one
          continue;
        }
      } else {
        // no previous match
        tj.EffPur = eff * pur;
        tj.mcpListIndex = tmpIndex;
        tjid[plane][tmpIndex] = tj.ID;
        ntht[plane][tmpIndex] = maxcnt;
      }
    } // tj

    if(neutrinoVxReconstructable) {
      // find the vertex closest to the true primary vertex
      float best = 1;
      unsigned short vx3ID = 0;
      if(!slc.pfps.empty()) {
        auto& pfp = slc.pfps[0];
        if((pfp.PDGCode == 14 || pfp.PDGCode == 12) && pfp.Vx3ID[0] > 0) {
          // Found a neutrino pfp with a start vertex
          auto& vx3 = slc.vtx3s[pfp.Vx3ID[0] - 1];
          // check the proximity to the true vertex
          Point3_t vpos = {{vx3.X, vx3.Y, vx3.Z}};
          if(PosSep(vpos, primVx) < 1) neutrinoPFPCorrect = true;
        } // neutrino pfp 
      } // PFParticles exist
      for(auto& vx3 : slc.vtx3s) {
        if(vx3.ID == 0) continue;
        if(vx3.TPCID != inTPCID) continue;
        Point3_t vpos = {{vx3.X, vx3.Y, vx3.Z}};
        float sep = PosSep(vpos, primVx);
        if(sep < best) {
          best = sep;
          vx3ID = vx3.ID;
        }
      } // vx3
      if(vx3ID > 0) vxReconstructedNearNuVx = true;
    } // neutrinoVxReconstructable

    if(tcc.matchTruth[1] > 0) {
      // print out
      mf::LogVerbatim myprt("TC");
      myprt<<"Number of primary particles "<<primMCPs.size()<<" Vtx";
      for(unsigned short xyz = 0; xyz < 3; ++xyz) myprt<<" "<<std::fixed<<std::setprecision(1)<<primVx[xyz];
      myprt<<" nuVx Reconstructable? "<<neutrinoVxReconstructable<<" vx near nuVx? "<<vxReconstructedNearNuVx;
      myprt<<" neutrinoPFPCorrect? "<<neutrinoPFPCorrect<<"\n";
      myprt<<"MCPIndx TrackId   PDG   eveID    KE    Len _______Dir_______  ____ProjInPln___               Process      StartHit-EndHit_nTruHits";
      for(unsigned int ipart = 0; ipart < slc.mcpList.size(); ++ipart) {
        bool doPrt = (std::find(mcpSelect.begin(), mcpSelect.end(), ipart) != mcpSelect.end());
        auto& mcp = slc.mcpList[ipart];
        // also print if this is a pizero or decay photon > 30 MeV
        if(mcp->PdgCode() == 111) doPrt = true;
        // Kinetic energy in MeV
        int TMeV = 1000 * (mcp->E() - mcp->Mass());
        if(mcp->PdgCode() == 22 && mcp->Process() == "Decay" && TMeV > 30) doPrt = true;
        if(!doPrt) continue;
        myprt<<"\n";
        myprt<<std::setw(7)<<ipart;
        myprt<<std::setw(8)<<mcp->TrackId();
        myprt<<std::setw(6)<<mcp->PdgCode();
        myprt<<std::setw(8)<<pi_serv->ParticleList().EveId(mcp->TrackId());
        myprt<<std::setw(6)<<TMeV;
        Point3_t start {{mcp->Vx(), mcp->Vy(), mcp->Vz()}};
        posOffsets = SCE->GetPosOffsets({start[0], start[1], start[2]});
        posOffsets.SetX(-posOffsets.X());
        start[0] += posOffsets.X();
        start[1] += posOffsets.Y();
        start[2] += posOffsets.Z();
        Point3_t end {{mcp->EndX(), mcp->EndY(), mcp->EndZ()}};
        posOffsets = SCE->GetPosOffsets({end[0], end[1], end[2]});
        posOffsets.SetX(-posOffsets.X());
        end[0] += posOffsets.X();
        end[1] += posOffsets.Y();
        end[2] += posOffsets.Z();
        myprt<<std::setw(7)<<std::setprecision(1)<<PosSep(start, end);
        Vector3_t dir {{mcp->Px(), mcp->Py(), mcp->Pz()}};
        SetMag(dir, 1);
        for(unsigned short xyz = 0; xyz < 3; ++xyz) myprt<<std::setw(6)<<std::setprecision(2)<<dir[xyz];
        std::vector<float> startWire(slc.nPlanes);
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = plane;
          auto tp = MakeBareTP(slc, start, dir, inCTP);
          myprt<<std::setw(6)<<tp.Delta;
          startWire[plane] = tp.Pos[0];
        } // plane
        myprt<<std::setw(22)<<mcp->Process();
        // print the extent of the particle in each TPC plane
        for(const geo::TPCID& tpcid : tcc.geom->IterateTPCIDs()) {
          geo::TPCGeo const& TPC = tcc.geom->TPC(tpcid);
          for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
            CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
            auto mcpHits = PutMCPHitsInVector(ipart, inCTP);
            if(mcpHits.empty()) continue;
            // get the direction correct
            auto& fhit = evt.allHits[slc.slHits[mcpHits[0]]];
            float fwire = fhit->WireID().Wire;
            auto& lhit = evt.allHits[slc.slHits[mcpHits.size() - 1]];
            float lwire = lhit->WireID().Wire;
            if(std::abs(fwire - startWire[plane]) < std::abs(lwire - startWire[plane])) {
              myprt<<" "<<PrintHitShort(slc.slHits[mcpHits[0]])<<"-"<<PrintHitShort(slc.slHits[mcpHits[mcpHits.size() - 1]]);
            } else {
              myprt<<" "<<PrintHitShort(slc.slHits[mcpHits[mcpHits.size() - 1]])<<"-"<<PrintHitShort(slc.slHits[mcpHits[0]]);
            }
            myprt<<"_"<<mcpHits.size();
          } // plane
        } // tpcid
      } // ipart
    } // tcc.matchTruth[1] > 0
    
    // Match Tjs and PFParticles and accumulate statistics
    MatchAndSum(slc, hist, mcpSelect);
        
    // Tj histograms
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.mcpListIndex == UINT_MAX) continue;
      auto& mcp = slc.mcpList[tj.mcpListIndex];
      short truIndex = PDGCodeIndex(slc, mcp->PdgCode());
      if(truIndex == SHRT_MAX) continue;
      float frac = 0;
      float cnt = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Environment[kEnvNearTj]) ++frac;
        ++cnt;
      } // ipt
      if(cnt > 0) frac /= cnt;
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      hist.fNearTj[truIndex]->Fill(TMeV, frac);
    } // tj
    
    
    // PFParticle histograms
    constexpr double twopi = 2 * M_PI;
    std::array<int, 5> recoCodeList = {{0, 11, 13, 211, 2212}};
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      // ignore showers
      if(pfp.PDGCode == 1111) continue;
      // require match to MC
      if(pfp.mcpListIndex == UINT_MAX) continue;
      auto& mcp = slc.mcpList[pfp.mcpListIndex];
      short truIndex = PDGCodeIndex(slc, mcp->PdgCode());
      if(truIndex == SHRT_MAX) continue;
      short recIndex = 0;
      for(recIndex = 0; recIndex < 5; ++recIndex) if(pfp.PDGCode == recoCodeList[recIndex]) break;
      if(recIndex > 4) continue;
      hist.PDGCode_reco_true->Fill((float)truIndex, (float)recIndex);
      // get the true start position and shift it by the SCE offset
      Point3_t truStart {{mcp->Vx(), mcp->Vy(), mcp->Vz()}};
      Point3_t truEnd {{mcp->EndX(), mcp->EndY(), mcp->EndZ()}};
      float truLen = PosSep(truStart, truEnd);
      if(truLen < 2) continue;
      posOffsets = SCE->GetPosOffsets({truStart[0], truStart[1], truStart[2]});
      posOffsets.SetX(-posOffsets.X());
      truStart[0] += posOffsets.X();
      truStart[1] += posOffsets.Y();
      truStart[2] += posOffsets.Z();
      auto recDir = pfp.Dir[0];
      short startEnd = 0;
      if(PosSep(pfp.XYZ[1], truStart) < PosSep(pfp.XYZ[0], truStart)) {
        startEnd = 1;
        for(unsigned short xyz = 0; xyz < 3; ++xyz) recDir[xyz] *= -1;
      }
      hist.fPFPStartEnd->Fill((float)startEnd);
      hist.fPFPStartdX[truIndex]->Fill(pfp.XYZ[startEnd][0] - truStart[0]);
      hist.fPFPStartdY[truIndex]->Fill(pfp.XYZ[startEnd][1] - truStart[1]);
      hist.fPFPStartdZ[truIndex]->Fill(pfp.XYZ[startEnd][2] - truStart[2]);
        Vector3_t truDir {{mcp->Px(), mcp->Py(), mcp->Pz()}};
      SetMag(truDir, 1);
      float dang = DeltaAngle(truDir, pfp.Dir[startEnd]);
      while(dang >  M_PI) dang -= twopi;
      while(dang < -M_PI) dang += twopi;
      hist.fPFPStartAngDiff[truIndex]->Fill(dang);
    } // pfp

    StudyElectrons(slc, hist);

  } // MatchTruth

  ////////////////////////////////////////////////
  void TruthMatcher::MatchAndSum(TCSlice& slc, const HistStuff& hist, const std::vector<unsigned int>& mcpSelect)
  {
    // Match Tjs and PFParticles and accumulate performance statistics
    if(mcpSelect.empty()) return;
    
    unsigned int tpc = slc.TPCID.TPC;
    unsigned int cstat = slc.TPCID.Cryostat;
    
    // get the hits associated with all MCParticles in mcpSelect
    std::vector<std::vector<unsigned int>> mcpHits(mcpSelect.size());
    for(unsigned short isel = 0; isel < mcpSelect.size(); ++isel) {
      unsigned int mcpIndex = mcpSelect[isel];
      for(unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
        if(slc.slHits[iht].mcpListIndex != mcpIndex) continue;
        if(slc.slHits[iht].ArtPtr->WireID().TPC != tpc) continue;
        if(slc.slHits[iht].ArtPtr->WireID().Cryostat != cstat) continue;
        mcpHits[isel].push_back(iht);
      } // iht
    } // mcpIndex
    
    // get the hits in all Tjs in this TPCID
    std::vector<std::vector<unsigned int>> tjHits(slc.tjs.size());
    for(unsigned short itj = 0; itj < slc.tjs.size(); ++itj) {
      auto& tj = slc.tjs[itj];
      if(tj.AlgMod[kKilled]) continue;
      if(DecodeCTP(tj.CTP).TPC != tpc) continue;
      tjHits[itj] = PutTrajHitsInVector(tj, kUsedHits);
      tj.mcpListIndex = UINT_MAX;
    } // itj
    
    // match them up
    for(unsigned short isel = 0; isel < mcpSelect.size(); ++isel) {
      if(mcpHits[isel].empty()) continue;
      unsigned int mcpIndex = mcpSelect[isel];
      auto& mcp = slc.mcpList[mcpIndex];
      unsigned short pdgIndex = PDGCodeIndex(slc, mcp->PdgCode());
      if(pdgIndex > 4) continue;
      std::string particleName = "Other";
      int pdg = abs(mcp->PdgCode());
      if(pdg == 11) particleName = "Electron";
      if(pdg == 22) particleName = "Photon";
      if(pdg == 13) particleName = "Muon";
      if(pdg == 211) particleName = "Pion";
      if(pdg == 321) particleName = "Kaon";
      if(pdg == 2212) particleName = "Proton";
      if(particleName == "Other") particleName = "PDG_" + std::to_string(pdg);
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        // put the mcpHits (which are already in this TPCID) in this plane in a vector
        std::vector<unsigned int> mcpPlnHits;
        for(auto iht : mcpHits[isel]) {
          auto& hit = slc.slHits[iht];
          if(hit.ArtPtr->WireID().Plane == plane) mcpPlnHits.push_back(iht);
        } // iht
        // require 2 truth-matched hits
        if(mcpPlnHits.size() < 2) continue;
        if((float)mcpPlnHits.size() >= tcc.matchTruth[3]) ++nLongInPln;
        TSums[pdgIndex] += TMeV;
        ++EPCnts[pdgIndex];
        CTP_t inCTP = EncodeCTP(cstat, tpc, plane);
        unsigned short mtj = USHRT_MAX;
        float maxEP = 0;
        for(unsigned short itj = 0; itj < slc.tjs.size(); ++itj) {
          // No hits in this TPC and plane
          if(tjHits[itj].empty()) continue;
          auto& tj = slc.tjs[itj];
          // wrong CTP?
          if(tj.CTP != inCTP) continue;
          // make a list of hits that are common
          auto shared = SetIntersection(mcpPlnHits, tjHits[itj]);
          if(shared.empty()) continue;
          float eff = (float)shared.size() / (float)mcpPlnHits.size();
          float pur = (float)shared.size() / (float)tjHits[itj].size();
          float ep = eff * pur;
          // temp for checking
          if(pdg != 11) {
            hist.fEff->Fill(eff);
            hist.fPur->Fill(pur);
          }
          // Replace a previously made poorer match with a better one?
          if(tj.mcpListIndex != UINT_MAX && ep > tj.EffPur) {
            tj.EffPur = ep;
            tj.mcpListIndex = mcpIndex;
          }
          if(ep > maxEP) {
            maxEP = ep;
            mtj = itj;
          }
        } // itj
        if(mtj == USHRT_MAX) {
          // failed to match MCParticle to a trajectory
          // Enter 0 in the profile histogram
          hist.fEP_T[pdgIndex]->Fill(TMeV, 0);
          if((float)mcpPlnHits.size() > tcc.matchTruth[3]) {
            ++nBadEP;
            mf::LogVerbatim myprt("TC");
            myprt<<particleName<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to isel "<<isel;
            myprt<<" nTrue hits "<<mcpPlnHits.size();
            myprt<<" extent "<<PrintHit(slc.slHits[mcpPlnHits[0]])<<"-"<<PrintHit(slc.slHits[mcpPlnHits[mcpPlnHits.size() - 1]]);
            myprt<<" events processed "<<evt.eventsProcessed;
          }
          continue;
        } // match failed
        auto& tj = slc.tjs[mtj];
        // don't clobber a better match
        if(maxEP < tj.EffPur) continue;
        tj.EffPur = maxEP;
        tj.mcpListIndex = mcpIndex;
        EPTSums[pdgIndex] += TMeV * tj.EffPur;
        hist.fEP_T[pdgIndex]->Fill(TMeV, tj.EffPur);
        // ignore electrons
        if(tj.EffPur < tcc.matchTruth[2] && (float)mcpPlnHits.size() >= tcc.matchTruth[3] && pdgIndex > 0) {
          ++nBadEP;
          mf::LogVerbatim myprt("TC");
          myprt<<particleName<<" BadEP: "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
          myprt<<" mcpIndex "<<mcpIndex;
          myprt<<" TMeV "<<(int)TMeV<<" MCP hits "<<mcpPlnHits.size();
          myprt<<" extent "<<PrintHit(slc.slHits[mcpPlnHits[0]])<<"-"<<PrintHit(slc.slHits[mcpPlnHits[mcpPlnHits.size() - 1]]);
          myprt<<" T"<<tj.ID;
          myprt<<" Algs";
          for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
          myprt<<" events processed "<<evt.eventsProcessed;
        } // print BadEP
        // badep
      } // plane
    } // isel

    // Calculate PFParticle efficiency and purity
    std::vector<std::vector<unsigned int>> pfpHits;
    if(!slc.pfps.empty()) {
      pfpHits.resize(slc.pfps.size());
      // get the hits in all pfparticles in this TPCID
      for(unsigned short ipfp = 0; ipfp < slc.pfps.size(); ++ipfp) {
        auto& pfp = slc.pfps[ipfp];
        if(pfp.ID == 0) continue;
        // in the right TPCID?
        if(pfp.TPCID != inTPCID) continue;
        // ignore the neutrino PFParticle and true photons
        if(pfp.PDGCode == 14 || pfp.PDGCode == 12 || pfp.PDGCode == 22) continue;
        pfp.mcpListIndex = UINT_MAX;
        for(auto& tjid : pfp.TjIDs) {
          unsigned short itj = tjid - 1;
          pfpHits[ipfp].insert(pfpHits[ipfp].end(), tjHits[itj].begin(), tjHits[itj].end());
        } // tj
      } // ipfp
    } // pfps exist
    
    // match them up
    for(unsigned short isel = 0; isel < mcpSelect.size(); ++isel) {
      unsigned short mpfp = USHRT_MAX;
      float maxEP = 0;
      unsigned int mcpIndex = mcpSelect[isel];
      if(mcpHits[isel].empty()) continue;
      auto& mcp = slc.mcpList[mcpIndex];
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      MCP_TSum += TMeV;
      // Performance reconstructing long muons, pions, kaons and protons
      unsigned short pdgIndex = PDGCodeIndex(slc, mcp->PdgCode());
      bool longMCP = (pdgIndex > 0 && pdgIndex < 5 && (float)mcpHits[isel].size() >= 2 * tcc.matchTruth[3]);
      if(longMCP) ++nLongMCP;
      for(unsigned short ipfp = 0; ipfp < slc.pfps.size(); ++ipfp) {
        auto& pfp = slc.pfps[ipfp];
        if(pfp.ID == 0) continue;
        // in the right TPCID?
        if(pfp.TPCID != inTPCID) continue;
        // not enough hits?
        if(pfpHits[ipfp].empty()) continue;
        auto shared = SetIntersection(mcpHits[isel], pfpHits[ipfp]);
        if(shared.empty()) continue;
        float eff = (float)shared.size() / (float)mcpHits[isel].size();
        float pur = (float)shared.size() / (float)pfpHits[ipfp].size();
        float ep = eff * pur;
        // Replace a previously made poorer match with a better one?
        if(pfp.mcpListIndex != UINT_MAX && ep > pfp.EffPur) {
          pfp.EffPur = ep;
          pfp.mcpListIndex = mcpIndex;
        }
        if(ep > maxEP) {
          maxEP = ep;
          mpfp = ipfp;
        }
      } // ipfp
      if(mpfp == USHRT_MAX) {
        if(TMeV > 30) {
          mf::LogVerbatim myprt("TC");
          myprt<<"BadPFP: MCParticle "<<mcpSelect[isel]<<" w PDGCode "<<mcp->PdgCode()<<" T "<<(int)TMeV<<" not reconstructed.";
          myprt<<" matched Tjs:";
          for(auto& tj : slc.tjs) {
            if(tj.AlgMod[kKilled]) continue;
            if(tj.mcpListIndex == mcpIndex) myprt<<" "<<tj.ID<<" EP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
          } // tj
          myprt<<" events processed "<<evt.eventsProcessed;
        } // TMeV > 30
        continue;
      }
      auto& pfp = slc.pfps[mpfp];
      if(maxEP < pfp.EffPur) continue;
      pfp.EffPur = maxEP;
      pfp.mcpListIndex = mcpIndex;
      MCP_EPTSum += TMeV * maxEP;
      ++MCP_PFP_Cnt;
      if(longMCP && maxEP > 0.8) ++nGoodLongMCP;
    } // isel
    
    MCP_Cnt += mcpSelect.size();

  } // MatchAndSum
    
  ////////////////////////////////////////////////
  void TruthMatcher::PrintResults(int eventNum) const
  {
    // Print performance metrics for each selected event
    
    mf::LogVerbatim myprt("TC");
    myprt<<"Evt "<<eventNum;
    float sum = 0;
    float sumt = 0;
    for(unsigned short pdgIndex = 0; pdgIndex < TSums.size(); ++pdgIndex) {
      if(TSums[pdgIndex] == 0) continue;
      if(pdgIndex == 0) myprt<<" El";
      if(pdgIndex == 1) myprt<<" Mu";
      if(pdgIndex == 2) myprt<<" Pi";
      if(pdgIndex == 3) myprt<<" K";
      if(pdgIndex == 4) myprt<<" P";
      float ave = EPTSums[pdgIndex] / (float)TSums[pdgIndex];
      myprt<<" "<<std::fixed<<std::setprecision(2)<<ave;
//      myprt<<" "<<EPCnts[pdgIndex];
      if(pdgIndex > 0) {
        sum  += TSums[pdgIndex];
        sumt += EPTSums[pdgIndex];
      }
    } // pdgIndex
    if(sum > 0) myprt<<" MuPiKP "<<std::fixed<<std::setprecision(2)<<sumt / sum;
    myprt<<" BadEP "<<nBadEP;
    if(nLongInPln > 0) {
      float longGood = 1 - (float)nBadEP / (float)nLongInPln;
      myprt<<" longGood "<<std::fixed<<std::setprecision(2)<<longGood;
    }
    if(MCP_TSum > 0) {
      // PFParticle statistics
      float ep = MCP_EPTSum / MCP_TSum;
      myprt<<" MCP cnt "<<(int)MCP_Cnt<<" PFP "<<std::fixed<<std::setprecision(2)<<ep;
    }
    if(Prim_TSum > 0) {
      float ep = Prim_EPTSum / Prim_TSum;
      myprt<<" PrimPFP "<<std::fixed<<std::setprecision(2)<<ep;
    }
    if(nLongMCP > 0) {
      float longGood = (float)nGoodLongMCP / (float)nLongMCP;
      myprt<<" longGood "<<std::fixed<<std::setprecision(2)<<longGood;
    }
    if(TruVxCounts[1] > 0) {
      // True vertex is reconstructable
      float frac = (float)TruVxCounts[2] / (float)TruVxCounts[1];
      myprt<<" NuVx correct "<<std::fixed<<std::setprecision(2)<<frac;
    }

  } // PrintResults
  
  ////////////////////////////////////////////////
  bool TruthMatcher::CanReconstruct(unsigned int mcpIndex, unsigned short nDimensions, const geo::TPCID& inTPCID)
  {
    // returns true if the MCParticle can be reconstructed in nDimensions
    if(mcpIndex > slc.mcpList.size() - 1) return false;
    if(nDimensions < 2 || nDimensions > 3) return false;
    
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for(auto& hit : slc.slHits) {
      if(hit.ArtPtr->WireID().TPC != inTPCID.TPC) continue;
      if(hit.ArtPtr->WireID().Cryostat != inTPCID.Cryostat) continue;
      if(hit.mcpListIndex == mcpIndex) ++cntInPln[hit.ArtPtr->WireID().Plane];
    } // hit
    unsigned short nPlnOK = 0;
    // Require at least 2 truth-matched hits in a plane
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] > 1) ++nPlnOK;
    if(nPlnOK < nDimensions - 1) return false;
    return true;
  } // CanReconstruct
  
  
  ////////////////////////////////////////////////
  std::vector<unsigned int> TruthMatcher::PutMCPHitsInVector(unsigned int mcpIndex, CTP_t inCTP)
  {
    // put the hits matched to the MCParticle into a vector in the requested CTP
    std::vector<unsigned int> hitVec;
    if(mcpIndex > slc.mcpList.size() - 1) return hitVec;
    geo::PlaneID planeID = DecodeCTP(inCTP);
    for(unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      auto& hit = slc.slHits[iht];
      if(hit.mcpListIndex != mcpIndex) continue;
      if(hit.ArtPtr->WireID().Plane != planeID.Plane) continue;
      if(hit.ArtPtr->WireID().TPC != planeID.TPC) continue;
      if(hit.ArtPtr->WireID().Cryostat != planeID.Cryostat) continue;
      hitVec.push_back(iht);
    } // iht
    return hitVec;
  } // PutMCPHitsInVector

  ////////////////////////////////////////////////
  void MCParticleListUtils::MakeTruTrajPoint(unsigned int MCParticleListIndex, TrajPoint& tp)
  {
    // Creates a trajectory point at the start of the MCParticle with index MCParticleListIndex. The
    // projected length of the MCParticle in the plane coordinate system is stored in TruTp.Delta.
    // The calling function should specify the CTP in which this TP resides.
    
    if(MCParticleListIndex > slc.mcpList.size() - 1) return;
    
    const simb::MCParticle* mcp = slc.mcpList[MCParticleListIndex];
    
    Point3_t pos {{mcp->Vx(), mcp->Vy(), mcp->Vz()}};
    Vector3_t dir {{mcp->Px(), mcp->Py(), mcp->Pz()}};
    SetMag(dir, 1);
    tp = MakeBareTP(slc, pos, dir, tp.CTP);
  } // MakeTruTrajPoint

  ////////////////////////////////////////////////
  ShowerStruct3D MCParticleListUtils::MakeCheatShower(TCSlice& slc, unsigned int mcpIndex, Point3_t primVx, int& truParentPFP)
  {
    // Make a 3D shower for an electron or photon MCParticle in the current TPC. 
    // The shower ID is set to 0 if there is a failure. The ID of the most likely parent pfp is returned
    // as well - the one that is closest to the primary vertex position
    auto ss3 = CreateSS3(slc, slc.TPCID);
    truParentPFP = 0;
    // save the ID 
    int goodID = ss3.ID;
    // set it to the failure state
    ss3.ID = 0;
    ss3.Cheat = true;
    if(mcpIndex > slc.mcpList.size() - 1) return ss3;
    auto& mcp = slc.mcpList[mcpIndex];
    int pdg = abs(mcp->PdgCode());
    if(!(pdg == 11 || pdg == 111)) return ss3;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    int eveID = mcp->TrackId();
    float showerEnergy = 1000 * mcp->E();
    float shMaxAlong = 7.0 * log(showerEnergy / 15);
    std::cout<<"MCS: showerEnergy "<<(int)showerEnergy<<" MeV. shMaxAlong "<<shMaxAlong<<" cm\n";
    
    // put the shower hits in two vectors. One for the InShower hits
    // and one for the shower parent
    std::vector<std::vector<unsigned int>> showerHits(slc.nPlanes);
    std::vector<std::vector<unsigned int>> parentHits(slc.nPlanes);
    unsigned short nsh = 0;
    unsigned short npar = 0;
    unsigned int cstat = slc.TPCID.Cryostat;
    unsigned int tpc = slc.TPCID.TPC;
    for(unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      auto& hit = slc.slHits[iht];
      if(hit.mcpListIndex > slc.mcpList.size()-1) continue;
      if(hit.ArtPtr->WireID().TPC != tpc) continue;
      if(hit.ArtPtr->WireID().Cryostat != cstat) continue;
      if(hit.Integral <= 0) continue;
      // require that it be used in a tj to ignore hits that are far
      // from the shower
      if(hit.InTraj <= 0) continue;
      auto& shmcp = slc.mcpList[hit.mcpListIndex];
      // look for an electron
      if(abs(shmcp->PdgCode()) != 11) continue;
      // with eve ID = the primary electron being considered
      int eid = pi_serv->ParticleList().EveId(shmcp->TrackId());
      if(eid != eveID) continue;
      if(shmcp->TrackId() == eid) {
        // store the shower parent hit
        parentHits[hit.ArtPtr->WireID().Plane].push_back(iht);
        ++npar;
      } else {
        // store the InShower hit
        showerHits[hit.ArtPtr->WireID().Plane].push_back(iht);
        ++nsh;
      }
    } // iht
    // require more hits in the shower than in the parent
    if(npar > nsh) return ss3;

 // create 2D cheat showers in each plane
    std::vector<int> dummy_tjlist;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(cstat, tpc, plane);
      auto ss = CreateSS(slc, inCTP, dummy_tjlist);
      // fill the ShPts
      ss.ShPts.resize(showerHits[plane].size());
      for(unsigned short iht = 0; iht < showerHits[plane].size(); ++iht) {
        auto& hit = slc.slHits[showerHits[plane][iht]];
        ss.ShPts[iht].HitIndex = showerHits[plane][iht];
        ss.ShPts[iht].TID = 0;
        ss.ShPts[iht].Chg = hit.Integral;
        ss.ShPts[iht].Pos[0] = hit.ArtPtr->WireID().Wire;
        ss.ShPts[iht].Pos[1] = hit.PeakTime * slc.unitsPerTick;
      } // iht
      ss.SS3ID = goodID;
      if(!UpdateShower("MCS", slc, ss, true)) {
        std::cout<<"Failed to update 2S"<<ss.ID<<"\n";
        return ss3;
      } // UpdateShower failed
      // remove shower points that are far away
      std::vector<ShowerPoint> spts;
      for(auto& spt : ss.ShPts) {
        // don't make a tight cut right now. The shower parameterization isn't great
        double along = slc.WirePitch * spt.RotPos[0];
        float tau = along / shMaxAlong;
//        auto& hit = slc.slHits[spt.HitIndex];
//        if(ss.ID == 2) std::cout<<"chk "<<PrintHit(hit)<<" along "<<(int)along<<" tau "<<std::fixed<<std::setprecision(1)<<tau<<"\n";
        if(tau > -1 && tau < 2) {
          spts.push_back(spt);
        } else {
//          std::cout<<"Skip "<<PrintHit(hit)<<" along "<<(int)along<<"\n";
        }
      } // spt
      std::cout<<" 2S"<<ss.ID<<" shpts "<<ss.ShPts.size()<<" new "<<spts.size()<<"\n";
      ss.ShPts = spts;
      ss.NeedsUpdate = true;
      if(!UpdateShower("MCS", slc, ss, true)) {
        std::cout<<"Failed to update 2S"<<ss.ID<<"\n";
        return ss3;
      } // UpdateShower failed
      if(!StoreShower("MCS", slc, ss)) {
        std::cout<<"Failed to store 2S"<<ss.ID<<"\n";
        return ss3;
      } // UpdateShower failed
      ss3.CotIDs.push_back(ss.ID);
    } // plane
    ss3.ID = goodID;
    if(!UpdateShower("MCS", slc, ss3, true)) {
      std::cout<<"SS3 Failed...\n";
      ss3.ID = 0;
      return ss3;
    }
    // success
    
    // now look for a parent pfp
    std::vector<std::pair<int, int>> pcnt;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      for(auto iht : parentHits[plane]) {
        auto& hit = slc.slHits[iht];
        if(hit.InTraj <= 0) continue;
        auto& tj = slc.tjs[hit.InTraj - 1];
        if(!tj.AlgMod[kMat3D]) continue;
        auto TInP = GetAssns(slc, "T", tj.ID, "P");
        if(TInP.empty()) continue;
        auto& pfp = slc.pfps[TInP[0] - 1];
        unsigned short indx = 0;
        for(indx = 0; indx < pcnt.size(); ++indx) if(pfp.ID == pcnt[indx].first) break;
        if(indx == pcnt.size()) {
          pcnt.push_back(std::make_pair(pfp.ID, 1));
        } else {
          ++pcnt[indx].second;
        }
      } // iht
    } // plane
    float close = 5;
    for(auto pcn : pcnt) {
      if(pcn.second < 6) continue;
      auto& pfp = slc.pfps[pcn.first - 1];
      for(unsigned short end = 0; end < 2; ++end) {
        float sep = PosSep(pfp.XYZ[end], primVx);
        if(sep > close) continue;
        close = sep;
        truParentPFP = pfp.ID;
      }
    } // pcn
    std::cout<<"truParent P"<<truParentPFP<<" close "<<std::fixed<<std::setprecision(2)<<close<<"\n";
    
    return ss3;
  } // MakeCheatShower
  
  /////////////////////////////////////////
  bool MCParticleListUtils::PrimaryElectronStart(Point3_t& start, Vector3_t& dir, float& energy)
  {
    // returns the SCE corrected start position of a primary electron
    if(slc.mcpList.empty()) return false;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    for(unsigned int part = 0; part < slc.mcpList.size(); ++part) {
      auto& mcp = slc.mcpList[part];
      // require electron
      if(abs(mcp->PdgCode()) != 11) continue;
      int eveID = pi_serv->ParticleList().EveId(mcp->TrackId());
      if(mcp->TrackId() != eveID) continue;
      start = {{mcp->Vx(), mcp->Vy(), mcp->Vz()}};
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      geo::Vector_t posOffsets = SCE->GetPosOffsets({start[0], start[1], start[2]});
      posOffsets.SetX(-posOffsets.X());
      start[0] += posOffsets.X();
      start[1] += posOffsets.Y();
      start[2] += posOffsets.Z();
      dir = {{mcp->Px(), mcp->Py(), mcp->Pz()}};
      SetMag(dir, 1);
      energy = 1000 * (mcp->E() - mcp->Mass());
      return true;
    } // part
    return false;
  } // PrimaryElectronStart
  
  
  /////////////////////////////////////////
  int MCParticleListUtils::PrimaryElectronPFPID(TCSlice& slc)
  {
    // Returns the ID of a pfp that has hits whose eve ID is a primary electron
    // and is the closest (< 5 WSE units) to the primary electron start
    if(slc.mcpList.empty()) return 0;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    TruthMatcher tm{slc};
    for(unsigned int part = 0; part < slc.mcpList.size(); ++part) {
      auto& mcp = slc.mcpList[part];
      // require electron
      if(abs(mcp->PdgCode()) != 11) continue;
      int eveID = pi_serv->ParticleList().EveId(mcp->TrackId());
      if(mcp->TrackId() != eveID) continue;
      Point3_t start = {{mcp->Vx(), mcp->Vy(), mcp->Vz()}};
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      geo::Vector_t posOffsets = SCE->GetPosOffsets({start[0], start[1], start[2]});
      posOffsets.SetX(-posOffsets.X());
      start[0] += posOffsets.X();
      start[1] += posOffsets.Y();
      start[2] += posOffsets.Z();
      // make a list of pfps that use tjs that use these hits
      std::vector<int> pfplist;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
        std::vector<unsigned int> mcphits = tm.PutMCPHitsInVector(part, inCTP);
        if(mcphits.empty()) continue;
        for(auto iht : mcphits) {
          auto& hit = slc.slHits[iht];
          if(hit.InTraj <= 0) continue;
          // require that the tj is 3D-matched
          auto& tj = slc.tjs[hit.InTraj - 1];
          if(!tj.AlgMod[kMat3D]) continue;
          // find the number of hits that are in mcphits
          auto tjhits = PutTrajHitsInVector(tj, kUsedHits);
          auto shared = SetIntersection(tjhits, mcphits);
          if(shared.size() < 10) continue;
          // find out what pfp it is used in.
          auto TInP = GetAssns(slc, "T", tj.ID, "P");
          if(TInP.size() != 1) continue;
          int pid = TInP[0];
          if(std::find(pfplist.begin(), pfplist.end(), pid) == pfplist.end()) pfplist.push_back(pid);
        } // iht
      } // plane
      if(pfplist.empty()) return 0;
      // Use the one that is closest to the true start position, not the
      // one that has the most matching hits. Electrons are likely to be
      // poorly reconstructed
      int pfpid = 0;
      float close = 1E6;
      for(auto pid : pfplist) {
        auto& pfp = slc.pfps[pid - 1];
        unsigned short nearEnd = 1 - FarEnd(slc, pfp, start);
        float sep = PosSep2(pfp.XYZ[nearEnd], start);
        if(sep > close) continue;
        close = sep;
        pfpid = pid;
      }
      return pfpid;
    } // part
    return 0;
  } // PrimaryElectronPFPID
  
  /////////////////////////////////////////
  int MCParticleListUtils::PrimaryElectronTjID(TCSlice& slc, CTP_t inCTP)
  {
    // returns the ID of a tj in inCTP that is closest to the start of a primary electron
    if(slc.mcpList.empty()) return 0;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    for(unsigned int part = 0; part < slc.mcpList.size(); ++part) {
      auto& mcp = slc.mcpList[part];
      // require electron
      if(abs(mcp->PdgCode()) != 11) continue;
      int eveID = pi_serv->ParticleList().EveId(mcp->TrackId());
      if(mcp->TrackId() != eveID) continue;
      return MCParticleStartTjID(part, inCTP);
    } // part
    return 0;
  } // PrimaryElectronTjID

  /////////////////////////////////////////
  int MCParticleListUtils::MCParticleStartTjID(unsigned int mcpIndex, CTP_t inCTP)
  {
    // Finds the trajectory that has hits matched to the MC Particle and is the closest to the
    // MCParticle start vertex
    
    if(mcpIndex > slc.mcpList.size() - 1) return 0;
    
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    // tj ID and occurrence count
    std::vector<std::array<int, 2>> t_cnt;
    for(auto hit : slc.slHits) {
      if(hit.InTraj <= 0) continue;
      if(hit.mcpListIndex != mcpIndex) continue;
      if(hit.ArtPtr->WireID().TPC != planeID.TPC) continue;
      if(hit.ArtPtr->WireID().Plane != planeID.Plane) continue;
      if(hit.ArtPtr->WireID().Cryostat != planeID.Cryostat) continue;
      unsigned short indx = 0;
      for(indx = 0; indx < t_cnt.size(); ++indx) if(t_cnt[indx][0] == hit.InTraj) break;
      if(indx == t_cnt.size()) {
        // didn't find the traj ID in t_cnt so add it
        std::array<int, 2> tmp;
        tmp[0] = hit.InTraj;
        tmp[1] = 1;
        t_cnt.push_back(tmp);
      } else {
        // count occurrences of this tj ID
        ++t_cnt[indx][1];
      }
    } // hit
    
    if(t_cnt.empty()) return 0;
    unsigned short occMax = 0;
    unsigned short occIndx = 0;
    for(unsigned short ii = 0; ii < t_cnt.size(); ++ii) {
      auto& tcnt = t_cnt[ii];
      if(tcnt[1] < occMax) continue;
      occMax = tcnt[1];
      occIndx = ii;
    } // tcnt
    return t_cnt[occIndx][0];
    
  } // MCParticleStartTj
  
  /////////////////////////////////////////
  unsigned int MCParticleListUtils::GetMCPListIndex(const TrajPoint& tp)
  {
    // Returns the MCParticle index that best matches the hits used in the Tp
    if(slc.mcpList.empty()) return UINT_MAX;
    if(tp.Chg <= 0) return UINT_MAX;
    std::vector<unsigned int> mcpIndex;
    std::vector<unsigned short> mcpCnt;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      unsigned int mcpi = slc.slHits[tp.Hits[ii]].mcpListIndex;
      if(mcpi == UINT_MAX) continue;
      unsigned short indx = 0;
      for(indx = 0; indx < mcpIndex.size(); ++indx) if(mcpi == mcpIndex[indx]) break;
      if(indx == mcpIndex.size()) {
        // not in the list so add it
        mcpIndex.push_back(mcpi);
        mcpCnt.push_back(1);
      } else {
        ++mcpCnt[indx];
      }
    } // ii
    if(mcpIndex.empty()) return UINT_MAX;
    if(mcpIndex.size() == 1) return mcpIndex[0];
    unsigned int indx = 0;
    unsigned short maxCnt = 0;
    for(unsigned short ii = 0; ii < mcpIndex.size(); ++ii) {
      if(mcpCnt[ii] > maxCnt) {
        maxCnt = mcpCnt[ii];
        indx = mcpIndex[ii];
      }
    } // ii
    return indx;
  } // GetMCPListIndex
  
  /////////////////////////////////////////
  unsigned int MCParticleListUtils::GetMCPListIndex(const ShowerStruct& ss, unsigned short& nTruHits)
  {
    // Returns the index of the MCParticle that has the most number of matches
    // to the hits in this shower
    
    if(slc.mcpList.empty()) return UINT_MAX;
    if(ss.TjIDs.empty()) return UINT_MAX;
    
    std::vector<unsigned int> pListCnt(slc.mcpList.size());
    
    for(auto& tjid : ss.TjIDs) {
      Trajectory& tj = slc.tjs[tjid - 1];
      for(auto& tp : tj.Pts) {
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          // ignore unmatched hits
          if(slc.slHits[iht].mcpListIndex > slc.mcpList.size() - 1) continue;
          ++pListCnt[slc.slHits[iht].mcpListIndex];
        } // ii
      } // pt
    } // tjid
    
    unsigned int pIndex = UINT_MAX;
    nTruHits = 0;
    for(unsigned short ii = 0; ii < pListCnt.size(); ++ii) {
      if(pListCnt[ii] > nTruHits) {
        nTruHits = pListCnt[ii];
        pIndex = ii;
      }
    } // ii
    
    return pIndex;
    
  } // GetMCPListIndex
  
  /////////////////////////////////////////
  unsigned int MCParticleListUtils::GetMCPListIndex(const Trajectory& tj, unsigned short& nTruHits)
  {
    // Returns the index of the MCParticle that has the most number of matches
    // to the hits in this trajectory
    
    if(slc.mcpList.empty()) return UINT_MAX;
    
    // Check all hits associated with this Tj
    std::vector<unsigned int> pListCnt(slc.mcpList.size());
    
    for(auto& tp : tj.Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        // ignore unmatched hits
        if(slc.slHits[iht].mcpListIndex > slc.mcpList.size() - 1) continue;
        ++pListCnt[slc.slHits[iht].mcpListIndex];
      } // ii
    } // pt
    
    unsigned int pIndex = UINT_MAX;
    nTruHits = 0;
    for(unsigned short ii = 0; ii < pListCnt.size(); ++ii) {
      if(pListCnt[ii] > nTruHits) {
        nTruHits = pListCnt[ii];
        pIndex = ii;
      }
    } // ii
    
    return pIndex;
    
  } // GetMCPListIndex
*/
} // namespace tca
