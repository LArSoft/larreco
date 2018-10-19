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
  
/* This code was used to develop the TMVA showerParentReader. The MakeCheatShower function needs
   to be re-written if this function is used in the future
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
*/
  //////////////////////////////////////////
  void TruthMatcher::MatchTruth(std::vector<simb::MCParticle*> const& mcpList, std::vector<unsigned int> const& mcpListIndex)
  {
    // Match trajectories, PFParticles, etc to the MC truth matched hits then
    // calculate reconstruction efficiency and purity. This function should only be
    // called once per event after reconstruction has been done in all slices
    
    // mcpList is a vector of all MCParticles that have been selected in the module
    if(mcpList.empty()) return;
    // mcpListIndex points to the MCParticle to which each hit is matched
    if(mcpListIndex.size() != (*evt.allHits).size()) return;
    
    
/* TODO: fix this later
    // Form a list of mother-daughter pairs that should be considered as a single particle
    std::vector<std::pair<unsigned int, unsigned int>> moda;
    for(unsigned int part = 0; part < mcpList.size(); ++part) {
      auto& mcp = mcpList[part];
      if(mcp->NumberDaughters() == 0) continue;
      unsigned int ndtr = 0;
      unsigned int dtrSelIndex = 0;
      for(int idtr = 0; idtr < mcp->NumberDaughters(); ++idtr) {
        int dtrTrackId = mcp->Daughter(idtr);
        // ignore it if it's not in the list
        bool ignore = true;
        for(unsigned short dpart = part + 1; dpart < mcpList.size(); ++dpart) {
          auto& dmcp = mcpList[dpart];
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
      auto& dtr = mcpList[dtrSelIndex];
      if(dtr->PdgCode() != mcp->PdgCode()) continue;
      if(tcc.matchTruth[1] > 1) mf::LogVerbatim("TC")<<"Daughter MCP "<<dtrSelIndex<<" -> Mother MCP "<<part;
      moda.push_back(std::make_pair(mcpListIndex[part], mcpListIndex[dtrSelIndex]));
    } // part
    if(!moda.empty()) {
      // over-write the hit -> daughter MCParticle association to the mother MCParticle.
      // Note that mother-daughter pairs are ordered by increasing generation. Reverse
      // moda so that the grand-daughters come first. Grand-daughters will then be
      // over-written by daughters in the previous generation
      if(moda.size() > 1) std::reverse(moda.begin(), moda.end());
      for(unsigned int iht = 0; iht < (*evt.allHits).size(); ++iht) {
        for(auto& md : moda) if(md.second == mcpListIndex[iht]) mcpListIndex[iht] = md.first;
      } // iht
    } // moda not empty
*/
    // decide if electrons inside showers should be associated with the eve electron
//    bool showerRecoMode = (tcc.showerTag[0] == 2) || (tcc.showerTag[0] == 4);
    
    MatchAndSum(mcpList, mcpListIndex);

  } // MatchTruth
/*
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
*/
  ////////////////////////////////////////////////
  void TruthMatcher::MatchAndSum(std::vector<simb::MCParticle*> const& mcpList, std::vector<unsigned int> const& mcpListIndex)
  {
    // Match Tjs and PFParticles and accumulate performance statistics
    
    // A MCParticle may span more than one TPC but trajectories and PFParticles are
    // reconstructed in only one TPC so we need to consider them separately
    for(const geo::TPCID& tpcid : tcc.geom->IterateTPCIDs()) {
      unsigned int tpc = tpcid.TPC;
      unsigned int cstat = tpcid.Cryostat;
      // select the MCParticles with matched hits in this TPC
      // mcpList   list of hits in this TPC
      std::vector<std::vector<unsigned int>> mcpHits(mcpList.size());
      bool hitsInTPC = false;
      for(unsigned int iht = 0; iht < (*evt.allHits).size(); ++iht) {
        if(mcpListIndex[iht] > mcpList.size() - 1) continue;
        auto& hit = (*evt.allHits)[iht];
        if(hit.WireID().Cryostat != cstat) continue;
        if(hit.WireID().TPC != tpc) continue;
        mcpHits[mcpListIndex[iht]].push_back(iht);
        hitsInTPC = true;
      } // iht
      // no sense continuing if there are no selected MCParticles that have hits
      // in this TPC
      if(!hitsInTPC) continue;
      // get the location of a tj in terms of (slice index, tj index) 
      std::vector<std::pair<unsigned short, unsigned short>> tjLocs;
      // and the hits
      std::vector<std::vector<unsigned int>> tjHits;
      // do the same for pfps
      std::vector<std::pair<unsigned short, unsigned short>> pfpLocs;
      std::vector<std::vector<unsigned int>> pfpHits;
      for(unsigned short isl = 0; isl < slices.size(); ++isl) {
        auto& slc = slices[isl];
        for(auto& tj : slc.tjs) {
          if(tj.AlgMod[kKilled]) continue;
          if(DecodeCTP(tj.CTP).TPC != tpc) continue;
          tj.mcpListIndex = UINT_MAX;
          tj.EffPur = 0;
          tjLocs.push_back(std::make_pair(isl, (unsigned short)(tj.ID-1)));
          tjHits.push_back(PutTrajHitsInVector(tj, kUsedHits));
        } // tj
        for(auto& pfp : slc.pfps) {
          if(pfp.ID <= 0) continue;
          if(pfp.TPCID != tpcid) continue;
          // ignore neutrino PFParticles
          if(pfp.PDGCode == 14 || pfp.PDGCode == 12) continue;
          pfp.mcpListIndex = UINT_MAX;
          pfp.EffPur = 0;
          pfpLocs.push_back(std::make_pair(isl, (unsigned short)(pfp.ID-1)));
          std::vector<unsigned int> tmp;
          for(auto tjid : pfp.TjIDs) {
            auto& tj = slc.tjs[tjid - 1];
            auto thits = PutTrajHitsInVector(tj, kUsedHits);
            tmp.insert(tmp.end(), thits.begin(), thits.end());
          } // tjid
          pfpHits.push_back(tmp);
        } // pfp
      } // slc
      unsigned short nplanes = tcc.geom->Nplanes(tpc, cstat);
      // match them
      for(unsigned int imcp = 0; imcp < mcpList.size(); ++imcp) {
        if(mcpHits[imcp].empty()) continue;
        // ignore if it isn't reconstructable in 3D
        if(!CanReconstruct(mcpHits[imcp], 3, tpcid)) continue;
        auto& mcp = mcpList[imcp];
        unsigned short pdgIndex = PDGCodeIndex(mcp->PdgCode());
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
        ++MCP_Cnt;
        MCP_TSum += TMeV;
        for(unsigned short plane = 0; plane < nplanes; ++plane) {
          // get the MCP hits in this plane 
          std::vector<unsigned int> mcpPlnHits;
          for(auto iht : mcpHits[imcp]) {
            auto& hit = (*evt.allHits)[iht];
            if(hit.WireID().Plane == plane) mcpPlnHits.push_back(iht);
          } // iht
          // require 2 truth-matched hits
          if(mcpPlnHits.size() < 2) continue;
          if((float)mcpPlnHits.size() >= tcc.matchTruth[3]) ++nLongInPln;
          TSums[pdgIndex] += TMeV;
          ++EPCnts[pdgIndex];
          // the tjSlIDs index of the tj that has the highest EP
          unsigned short mtjLoc = USHRT_MAX;
          float maxEP = 0;
          for(unsigned short iit = 0; iit < tjLocs.size(); ++iit) {
            // check for hits in this TPC and plane
            if(tjHits[iit].empty()) continue;
            // find hits that are common
            auto shared = SetIntersection(mcpPlnHits, tjHits[iit]);
            if(shared.empty()) continue;
            float eff = (float)shared.size() / (float)mcpPlnHits.size();
            float pur = (float)shared.size() / (float)tjHits[iit].size();
            float ep = eff * pur;
            if(ep > maxEP) {
              maxEP = ep;
              mtjLoc = iit;
            } // ep > tj.EffPur
          } // iit
          if(mtjLoc == USHRT_MAX) {
            if((float)mcpPlnHits.size() > tcc.matchTruth[3]) {
              ++nBadEP;
              mf::LogVerbatim myprt("TC");
              myprt<<particleName<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to imcp "<<imcp;
              myprt<<" nTrue hits "<<mcpPlnHits.size();
//              myprt<<" extent "<<PrintHit(slc.slHits[mcpPlnHits[0]])<<"-"<<PrintHit(slc.slHits[mcpPlnHits[mcpPlnHits.size() - 1]]);
              myprt<<" events processed "<<evt.eventsProcessed;
            } // BadEP
          } else {
            // set EffPur for the best matching tj
            auto& tj = slices[tjLocs[mtjLoc].first].tjs[tjLocs[mtjLoc].second];
            if(maxEP > tj.EffPur) {
              tj.EffPur = maxEP;
              tj.mcpListIndex = imcp;
              EPTSums[pdgIndex] += TMeV * tj.EffPur;
            }
            // print BadEP ignoring electrons
            if(tj.EffPur < tcc.matchTruth[2] && (float)mcpPlnHits.size() >= tcc.matchTruth[3] && pdgIndex > 0) {
              ++nBadEP;
              mf::LogVerbatim myprt("TC");
              myprt<<particleName<<" BadEP: "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
              myprt<<" imcp "<<imcp;
              myprt<<" TMeV "<<(int)TMeV<<" MCP hits "<<mcpPlnHits.size();
              //            myprt<<" extent "<<PrintHit(slc.slHits[mcpPlnHits[0]])<<"-"<<PrintHit(slc.slHits[mcpPlnHits[mcpPlnHits.size() - 1]]);
              myprt<<" T"<<tj.ID;
              myprt<<" Algs";
              for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
              myprt<<" events processed "<<evt.eventsProcessed;
            } // print BadEP
          } // matched tj in this plane
        } // plane
        // pfp matching
        bool longMCP = (pdgIndex > 0 && pdgIndex < 5 && (float)mcpHits[imcp].size() >= 2 * tcc.matchTruth[3]);
        float maxEP = 0;
        unsigned short mpfpLoc = USHRT_MAX;
        for(unsigned short iit = 0; iit < pfpLocs.size(); ++iit) {
          // check for hits in this TPC and plane
          if(pfpHits[iit].empty()) continue;
          // find hits that are common
          auto shared = SetIntersection(mcpHits[imcp], pfpHits[iit]);
          if(shared.empty()) continue;
          float eff = (float)shared.size() / (float)mcpHits[imcp].size();
          float pur = (float)shared.size() / (float)pfpHits[iit].size();
          float ep = eff * pur;
          if(ep > maxEP) {
            maxEP = ep;
            mpfpLoc = iit;
          } // ep > tj.EffPur
        } // iit
        if(mpfpLoc == USHRT_MAX) {
          // no matching pfp
          if(TMeV > 30) {
            mf::LogVerbatim myprt("TC");
            myprt<<"BadPFP: MCParticle "<<imcp<<" w PDGCode "<<mcp->PdgCode()<<" T "<<(int)TMeV<<" not reconstructed.";
            myprt<<" events processed "<<evt.eventsProcessed;
          } // TMeV > 30
        } else {
          // matched pfp
          auto& pfp = slices[pfpLocs[mpfpLoc].first].pfps[pfpLocs[mpfpLoc].second];
          if(maxEP > pfp.EffPur) {
            pfp.EffPur = maxEP;
            pfp.mcpListIndex = imcp;
            MCP_EPTSum += TMeV * maxEP;
            ++MCP_PFP_Cnt;
            if(longMCP && maxEP > 0.8) ++nGoodLongMCP;
          } // maxEP > pfp.EffPur
        } // matched pfp
      } // imcp
    } // tpcid
    
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
  bool TruthMatcher::CanReconstruct(std::vector<unsigned int> mcpHits, unsigned short nDimensions, const geo::TPCID& inTPCID)
  {
    // returns true if the MCParticle that is matched to the hits in mcpHits can be reconstructed
    // in nDimensions in inTPCID
    if(mcpHits.empty()) return false;
    if(nDimensions < 2 || nDimensions > 3) return false;
    unsigned short tpc = inTPCID.TPC;
    unsigned short cstat = inTPCID.Cryostat;
    unsigned short nplanes = tcc.geom->Nplanes(tpc, cstat);
    std::vector<unsigned short> cntInPln(nplanes);
    for(auto iht : mcpHits) {
      if(iht > (*evt.allHits).size()) return false;
      auto& hit = (*evt.allHits)[iht];
      if(hit.WireID().TPC != tpc) continue;
      if(hit.WireID().Cryostat != cstat) continue;
      ++cntInPln[hit.WireID().Plane];
    } // hit
    unsigned short nPlnOK = 0;
    // Require at least 2 truth-matched hits in a plane
    for(unsigned short plane = 0; plane < nplanes; ++plane) if(cntInPln[plane] > 1) ++nPlnOK;
    return (nPlnOK >= 2);
  } // CanReconstruct
} // namespace tca
