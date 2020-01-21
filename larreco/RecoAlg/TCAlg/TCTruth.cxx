#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h"
#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <stdlib.h>
#include <string>
#include <utility>

namespace tca {

  //////////////////////////////////////////
  void TruthMatcher::Initialize()
  {
     // Initialize the variables used to calculate Efficiency * Purity (aka EP) for matching to truth
    EPCnts.fill(0);
    TSums.fill(0.0);
    EPTSums.fill(0.0);
    nBadT = 0;
    nBadP = 0;
    art::ServiceHandle<art::TFileService const> tfs;
    hist.CreateHists(tfs);
  } // Initialize

  //////////////////////////////////////////
  void TruthMatcher::MatchTruth()
  {
    // Match trajectories, PFParticles, etc to the MC truth matched hits then
    // calculate reconstruction efficiency and purity. This function should only be
    // called once per event after reconstruction has been done in all slices

    // check for a serious error
    if(!evt.mcpHandle) return;
    // and no MCParticles
    if((*evt.mcpHandle).empty()) return;
    if(slices.empty()) return;
    if(evt.allHitsMCPIndex.size() != (*evt.allHits).size()) return;

    MatchTAndSum();
    
    // Set the PFP mcpIndex and EP
    for(auto& slc : slices) {
      for(auto& pfp : slc.pfps) {
        if(pfp.ID <= 0) continue;
        float hitCnt = 0;
        float EPSum = 0;
        unsigned int mcpIndex = UINT_MAX;
        for(auto tid : pfp.TjIDs) {
          if(tid <= 0) continue;
          auto& tj = slc.tjs[tid - 1];
          if(tj.mcpIndex == UINT_MAX) continue;
          float npwc = NumPtsWithCharge(slc, tj, false);
          hitCnt += npwc;
          EPSum += npwc * tj.EffPur;
          if(mcpIndex == UINT_MAX) mcpIndex = tj.mcpIndex;
        } // tid
        if(hitCnt < 3) continue;
        pfp.mcpIndex = mcpIndex;
        pfp.EffPur = EPSum / hitCnt;
      } // pfp
    } // slc
    
  } // MatchTruth

  ////////////////////////////////////////////////
  void TruthMatcher::MatchTAndSum()
  {
    // match tjs in all tpcs that were reconstructed and sum
    
    // create a list of TPCs that were reconstructed
    std::vector<unsigned int> tpcList;
    for(auto& slc : slices) {
      unsigned int tpc = slc.TPCID.TPC;
      if(std::find(tpcList.end(), tpcList.end(), tpc) == tpcList.end()) tpcList.push_back(tpc);
    } // slc
    if(tpcList.empty()) return;
    
    // Hit -> T unique ID in all slices
    std::vector<int> inTUID((*evt.allHits).size(), 0);
    for(auto& slc : slices) {
      if(std::find(tpcList.begin(), tpcList.end(), slc.TPCID.TPC) == tpcList.end()) continue;
      for(auto& slh : slc.slHits) if(slh.InTraj > 0) {
        auto& tj = slc.tjs[slh.InTraj - 1];
        inTUID[slh.allHitsIndex] = tj.UID;
      }
    } // slc

    for(const geo::TPCID& tpcid : tcc.geom->IterateTPCIDs()) {
      // ignore protoDUNE dummy TPCs
      if(tcc.geom->TPC(tpcid).DriftDistance() < 25.0) continue;
      unsigned int tpc = tpcid.TPC;
      if(std::find(tpcList.begin(), tpcList.end(), tpc) == tpcList.end()) continue;
      // iterate over planes
      for(unsigned short plane = 0; plane < tcc.geom->Nplanes(); ++plane) {
        // form a list of MCParticles in this TPC and plane and the hit count
        std::vector<std::pair<unsigned int, float>> mCnt;
        // and lists of tj UIDs that use these hits
        std::vector<std::vector<std::pair<int, float>>> mtCnt;
        for(unsigned int iht = 0; iht < (*evt.allHits).size(); ++iht) {
          // require that it is MC-matched
          if(evt.allHitsMCPIndex[iht] == UINT_MAX) continue;
          // require that it resides in this tpc and plane
          auto& hit = (*evt.allHits)[iht];
          if(hit.WireID().TPC != tpc) continue;
          if(hit.WireID().Plane != plane) continue;
          unsigned int mcpi = evt.allHitsMCPIndex[iht];
          // find the mCnt entry
          unsigned short indx = 0;
          for(indx = 0; indx < mCnt.size(); ++indx) if(mcpi == mCnt[indx].first) break;
          if(indx == mCnt.size()) {
            mCnt.push_back(std::make_pair(mcpi, 0));
            mtCnt.resize(mCnt.size());
          }
          ++mCnt[indx].second;
          // see if it is used in a tj
          if(inTUID[iht] <= 0) continue;
          // find the mtCnt entry
          unsigned short tindx = 0;
          for(tindx = 0; tindx < mtCnt[indx].size(); ++tindx) if(mtCnt[indx][tindx].first == inTUID[iht]) break;
          if(tindx == mtCnt[indx].size()) mtCnt[indx].push_back(std::make_pair(inTUID[iht], 0));
          ++mtCnt[indx][tindx].second;
        } // iht
        if(mCnt.empty()) continue;
        for(unsigned short indx = 0; indx < mCnt.size(); ++indx) {
          // require at least 3 hits per plane to reconstruct
          if(mCnt[indx].second < 3) continue;
          // get a reference to the MCParticle and ensure that it is one we want to track
          auto& mcp = (*evt.mcpHandle)[mCnt[indx].first];
          float TMeV = 1000 * (mcp.E() - mcp.Mass());
          // don't weight entries by T
          float wght = 1;
          unsigned short pdgIndex = PDGCodeIndex(mcp.PdgCode());
          if(pdgIndex > 4) continue;
          hist.fTruT[pdgIndex]->Fill(TMeV);
          TSums[pdgIndex] += wght;
          ++EPCnts[pdgIndex];
          int pdg = abs(mcp.PdgCode());
          // find the tj with the highest match count
          std::pair<int, float> big = std::make_pair(0, 0);
          for(unsigned short tindx = 0; tindx < mtCnt[indx].size(); ++tindx) {
            if(mtCnt[indx][tindx].second > big.second) big = mtCnt[indx][tindx];
          } // tindx
          if(big.first == 0) continue;
          auto slcIndex = GetSliceIndex("T", big.first);
          if(slcIndex.first == USHRT_MAX) continue;
          auto& slc = slices[slcIndex.first];
          auto& tj = slc.tjs[slcIndex.second];
          auto tHits = PutTrajHitsInVector(tj, kUsedHits);
          float npwc = tHits.size();
          float eff = big.second / mCnt[indx].second;
          hist.fEff_T[pdgIndex]->Fill(TMeV, eff);
          float pur = big.second / npwc;
          hist.fPur_T[pdgIndex]->Fill(TMeV, pur);
          tj.EffPur = eff * pur;
          tj.mcpIndex = mCnt[indx].first;
          EPTSums[pdgIndex] += wght * tj.EffPur;
          // print BadEP ignoring electrons
          if(tj.EffPur < tcc.matchTruth[2] && (float)mCnt[indx].second > tcc.matchTruth[3] && pdgIndex > 0) {
            ++nBadT;
            std::string particleName = "Other";
            if(pdg == 11) particleName = "Electron";
            if(pdg == 22) particleName = "Photon";
            if(pdg == 13) particleName = "Muon";
            if(pdg == 211) particleName = "Pion";
            if(pdg == 321) particleName = "Kaon";
            if(pdg == 2212) particleName = "Proton";
            mf::LogVerbatim myprt("TC");
            myprt<<"badT"<<tj.ID<<" TU"<<tj.UID<<" tpc "<<tpc;
            myprt<<" slice index "<<slcIndex.first;
            myprt<<" -> mcp "<<tj.mcpIndex<<" with mCnt = "<<(int)mCnt[indx].second;
            myprt<<" "<<particleName<<" T = "<<(int)TMeV<<" MeV";
            myprt<<" EP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
            // print the first and last hit
            unsigned int firstHit = UINT_MAX;
            unsigned int lastHit = 0;
            for(unsigned int iht = 0; iht < (*evt.allHits).size(); ++iht) {
              // require that it is matched to this MCP
              if(evt.allHitsMCPIndex[iht] != mCnt[indx].first) continue;
              // require that it resides in this tpc and plane
              auto& hit = (*evt.allHits)[iht];
              if(hit.WireID().TPC != tpc) continue;
              if(hit.WireID().Plane != plane) continue;
              if(firstHit == UINT_MAX) firstHit = iht;
              lastHit = iht;
            } // iht
            auto& fhit = (*evt.allHits)[firstHit];
            myprt<<" Hit range "<<fhit.WireID().Plane<<":"<<fhit.WireID().Wire<<":"<<(int)fhit.PeakTime();
            auto& lhit = (*evt.allHits)[lastHit];
            myprt<<" - "<<lhit.WireID().Plane<<":"<<lhit.WireID().Wire<<":"<<(int)lhit.PeakTime();
            myprt<<" evts Processed "<<evt.eventsProcessed;
            myprt<<" Evt "<<evt.event;
          } // Poor EP
        } // indx
      } // plane
    } // tpcid
  } // MatchTAndSum

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
      myprt<<" "<<std::fixed<<std::setprecision(3)<<ave;
      if(pdgIndex > 0) {
        sum  += TSums[pdgIndex];
        sumt += EPTSums[pdgIndex];
      }
    } // pdgIndex
    if(sum > 0) myprt<<" MuPiKP "<<std::fixed<<std::setprecision(3)<<sumt / sum;
    myprt<<" nBadT "<<(int)nBadT;
    if(MCP_TSum > 0) {
      // PFParticle statistics
      float ep = MCP_EPTSum / MCP_TSum;
      myprt<<" MCP cnt "<<(int)MCP_Cnt<<" PFP EP "<<std::fixed<<std::setprecision(3)<<ep;
    }
    if(tcc.match3DCuts[0] > 0) { myprt<<" +Mat3D"; } else { myprt<<" -Mat3D"; }
    myprt<<" MCP Cnt:";
    for(unsigned short pdgIndex = 0; pdgIndex < TSums.size(); ++pdgIndex) {
      if(pdgIndex == 0) myprt<<" El";
      if(pdgIndex == 1) myprt<<" Mu";
      if(pdgIndex == 2) myprt<<" Pi";
      if(pdgIndex == 3) myprt<<" K";
      if(pdgIndex == 4) myprt<<" P";
      myprt<<" "<<EPCnts[pdgIndex];
    } // pdgIndex
  } // PrintResults

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
   art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
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
} // namespace tca
