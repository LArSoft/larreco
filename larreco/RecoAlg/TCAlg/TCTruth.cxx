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
    std::cout<<"TruthMatcher::Initialize isn't functional anymore\n";
/*
    EPCnts.fill(0);
    TSums.fill(0.0);
    EPTSums.fill(0.0);
    nBadT = 0;
    nBadP = 0;
    art::ServiceHandle<art::TFileService const> tfs;
    hist.CreateHists(tfs);
*/
  } // Initialize

  //////////////////////////////////////////
  void TruthMatcher::MatchTruth()
  {
    // Match trajectories, PFParticles, etc to the MC truth matched hits then
    // calculate reconstruction efficiency and purity. This function should only be
    // called once per event after reconstruction has been done in all slices

    std::cout<<"TruthMatcher::MatchTruth isn't functional anymore\n";
/*

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
*/
  } // MatchTruth
/*
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
*/
} // namespace tca
