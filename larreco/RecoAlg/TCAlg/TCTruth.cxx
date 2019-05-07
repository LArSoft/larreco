#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"


#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
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

    MatchAndSum();
/*
    // print electron likelihood to output to create an ntuple
    if(tcc.modes[kStudy2]) {
      for(auto& slc : slices) {
        for(auto& tj : slc.tjs) {
          if(tj.AlgMod[kKilled]) continue;
          if(tj.mcpListIndex != 0) continue;
          auto& mcp = mcpList[tj.mcpListIndex];
          int pdg = abs(mcp->PdgCode());
          short TMeV = 1000 * (mcpList[0]->E() - mcpList[0]->Mass());
          float eLike = ElectronLikelihood(slc, tj);
          mf::LogVerbatim myprt("TC");
          myprt<<"ntp "<<pdg<<" "<<TMeV;
          myprt<<" "<<tj.MCSMom;
          myprt<<" "<<tj.PDGCode;
          myprt<<" "<<std::fixed<<std::setprecision(1);
          myprt<<" "<<TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
          myprt<<" "<<std::fixed<<std::setprecision(2);
          myprt<<" "<<std::setprecision(3)<<tj.ChgRMS;
          myprt<<" "<<eLike;
          myprt<<" "<<tj.EffPur;
          if(pdg == 13 && eLike > 0.5) mf::LogVerbatim("TC")<<"Bad mu "<<eLike<<" "<<evt.eventsProcessed;
          if(pdg == 11 && eLike < 0.5) mf::LogVerbatim("TC")<<"Bad el "<<eLike<<" "<<evt.eventsProcessed;
        } // tj
      } // slc
    } // study2
*/
  } // MatchTruth

  ////////////////////////////////////////////////
  void TruthMatcher::MatchAndSum()
  {
    // Match Tjs and PFParticles and accumulate performance statistics

    if(!evt.mcpHandle) return;
    if(evt.allHitsMCPIndex.size() != (*evt.allHits).size()) return;

    // A MCParticle may span more than one TPC but trajectories and PFParticles are
    // reconstructed in only one TPC so we need to consider them separately
    for(const geo::TPCID& tpcid : tcc.geom->IterateTPCIDs()) {
      unsigned int tpc = tpcid.TPC;
      unsigned int cstat = tpcid.Cryostat;
      // find the MCParticles with matched hits in this TPC
      std::vector<unsigned int> mcpIndex;
      // mcpList   list of hits in this TPC
      std::vector<std::vector<unsigned int>> mcpHits;
      for(unsigned int iht = 0; iht < (*evt.allHits).size(); ++iht) {
        if(evt.allHitsMCPIndex[iht] == UINT_MAX) continue;
        auto& hit = (*evt.allHits)[iht];
        if(hit.WireID().Cryostat != cstat) continue;
        if(hit.WireID().TPC != tpc) continue;
        unsigned int indx = 0;
        for(indx = 0; indx < mcpIndex.size(); ++indx) if(evt.allHitsMCPIndex[iht] == mcpIndex[indx]) break;
        if(indx == mcpIndex.size()) {
          mcpIndex.push_back(evt.allHitsMCPIndex[iht]);
          mcpHits.resize(mcpIndex.size());
        }
        mcpHits[indx].push_back(iht);
      } // iht
      // no sense continuing if there are no selected MCParticles that have hits
      // in this TPC
      if(mcpIndex.empty()) continue;
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
          tj.mcpIndex = UINT_MAX;
          tj.EffPur = 0;
          tjLocs.push_back(std::make_pair(isl, (unsigned short)(tj.ID-1)));
          // get the tj hits which are indexed into slHits
          auto thits = PutTrajHitsInVector(tj, kUsedHits);
          // convert to the index into allHits
          for(unsigned short ii = 0; ii < thits.size(); ++ii) thits[ii] = slc.slHits[thits[ii]].allHitsIndex;
          tjHits.push_back(thits);
        } // tj
        for(auto& pfp : slc.pfps) {
          if(pfp.ID <= 0) continue;
          if(pfp.TPCID != tpcid) continue;
          // ignore neutrino PFParticles
          if(pfp.PDGCode == 14 || pfp.PDGCode == 12) continue;
          pfp.mcpIndex = UINT_MAX;
          pfp.EffPur = 0;
          pfpLocs.push_back(std::make_pair(isl, (unsigned short)(pfp.ID-1)));
          std::vector<unsigned int> tmp;
          for(auto tjid : pfp.TjIDs) {
            auto& tj = slc.tjs[tjid - 1];
            auto thits = PutTrajHitsInVector(tj, kUsedHits);
            for(unsigned short ii = 0; ii < thits.size(); ++ii) thits[ii] = slc.slHits[thits[ii]].allHitsIndex;
            tmp.insert(tmp.end(), thits.begin(), thits.end());
          } // tjid
          pfpHits.push_back(tmp);
        } // pfp
      } // slc
      if(pfpHits.size() != pfpLocs.size()) {
        std::cout<<"TCTruth coding error. pfpHits size "<<pfpHits.size()<<" != pfpLoocs size "<<pfpLocs.size()<<"\n";
        exit(1);
      }
      unsigned short nplanes = tcc.geom->Nplanes(tpc, cstat);
      // match them
      for(unsigned int imcp = 0; imcp < mcpIndex.size(); ++imcp) {
        if(mcpHits[imcp].empty()) continue;
        // ignore if it isn't reconstructable in 3D
        if(!CanReconstruct(mcpHits[imcp], 3, tpcid)) continue;
        auto& mcp = (*evt.mcpHandle)[mcpIndex[imcp]];
        unsigned short pdgIndex = PDGCodeIndex(mcp.PdgCode());
        if(pdgIndex > 4) continue;
        std::string particleName = "Other";
        int pdg = abs(mcp.PdgCode());
        if(pdg == 11) particleName = "Electron";
        if(pdg == 22) particleName = "Photon";
        if(pdg == 13) particleName = "Muon";
        if(pdg == 211) particleName = "Pion";
        if(pdg == 321) particleName = "Kaon";
        if(pdg == 2212) particleName = "Proton";
        if(particleName == "Other") particleName = "PDG_" + std::to_string(pdg);
        float TMeV = 1000 * (mcp.E() - mcp.Mass());
        ++MCP_Cnt;
        MCP_TSum += TMeV;
        for(unsigned short plane = 0; plane < nplanes; ++plane) {
          // get the MCP hits in this plane
          std::vector<unsigned int> mcpPlnHits;
          unsigned int firstHit = 0;
          unsigned int firstWire = USHRT_MAX;
          unsigned int lastHit = 0;
          unsigned int lastWire = 0;
          for(auto iht : mcpHits[imcp]) {
            auto& hit = (*evt.allHits)[iht];
            if(hit.WireID().Plane != plane) continue;
            mcpPlnHits.push_back(iht);
            if(hit.WireID().Wire < firstWire) {
              firstWire = hit.WireID().Wire;
              firstHit = iht;
            }
            if(hit.WireID().Wire > lastWire) {
              lastWire = hit.WireID().Wire;
              lastHit = iht;
            }
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
              myprt<<" in pln "<<plane;
              myprt<<" nTrue hits "<<mcpPlnHits.size();
              auto& fhit = (*evt.allHits)[firstHit];
              myprt<<" extent "<<fhit.WireID().Plane<<":"<<fhit.WireID().Wire<<":"<<(int)(fhit.PeakTime());
              auto& lhit = (*evt.allHits)[lastHit];
              myprt<<" - "<<lhit.WireID().Plane<<":"<<lhit.WireID().Wire<<":"<<(int)(lhit.PeakTime());
              myprt<<" events processed "<<evt.eventsProcessed;
            } // BadEP
          } else {
            // set EffPur for the best matching tj
            auto& tj = slices[tjLocs[mtjLoc].first].tjs[tjLocs[mtjLoc].second];
            if(maxEP > tj.EffPur) {
              tj.EffPur = maxEP;
              tj.mcpIndex = mcpIndex[imcp];
              EPTSums[pdgIndex] += TMeV * tj.EffPur;
            }
            // print BadEP ignoring electrons
            if(tj.EffPur < tcc.matchTruth[2] && (float)mcpPlnHits.size() >= tcc.matchTruth[3] && pdgIndex > 0) {
              ++nBadEP;
              mf::LogVerbatim myprt("TC");
              myprt<<particleName<<" BadEP: "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
              myprt<<" imcp "<<imcp;
              myprt<<" in pln "<<plane;
              myprt<<" TMeV "<<(int)TMeV<<" MCP hits "<<mcpPlnHits.size();
              auto& fhit = (*evt.allHits)[firstHit];
              myprt<<" extent "<<fhit.WireID().Plane<<":"<<fhit.WireID().Wire<<":"<<(int)(fhit.PeakTime());
              auto& lhit = (*evt.allHits)[lastHit];
              myprt<<" - "<<lhit.WireID().Plane<<":"<<lhit.WireID().Wire<<":"<<(int)(lhit.PeakTime());
              myprt<<" T"<<tj.ID;
              myprt<<" Algs";
              for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
              myprt<<" events processed "<<evt.eventsProcessed;
            } // print BadEP
          } // matched tj in this plane
        } // plane
        // pfp matching
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
        if(mpfpLoc >= pfpLocs.size()) continue;
        // matched pfp
        auto& pfp = slices[pfpLocs[mpfpLoc].first].pfps[pfpLocs[mpfpLoc].second];
        if(maxEP > pfp.EffPur) {
          pfp.EffPur = maxEP;
          pfp.mcpIndex = mcpIndex[imcp];
        } // maxEP > pfp.EffPur
      } // imcp
      
      // accumulate
      for(unsigned int imcp = 0; imcp < mcpIndex.size(); ++imcp) {
        if(mcpHits[imcp].empty()) continue;
        // ignore if it isn't reconstructable in 3D
        if(!CanReconstruct(mcpHits[imcp], 3, tpcid)) continue;
        auto& mcp = (*evt.mcpHandle)[mcpIndex[imcp]];
        unsigned short pdgIndex = PDGCodeIndex(mcp.PdgCode());
        if(pdgIndex == 0) continue;
        int pdg = abs(mcp.PdgCode());
        // ignore electrons
        if(pdg == 11) continue;
        bool longMCP = (pdgIndex > 0 && pdgIndex < 5 && (float)mcpHits[imcp].size() >= 2 * tcc.matchTruth[3]);
        float TMeV = 1000 * (mcp.E() - mcp.Mass());
        float maxEP = 0;
        unsigned short mpfpLoc = USHRT_MAX;
        for(unsigned short iit = 0; iit < pfpLocs.size(); ++iit) {
          if(pfpHits[iit].empty()) continue;
          auto& pfp = slices[pfpLocs[iit].first].pfps[pfpLocs[iit].second];
          if(pfp.mcpIndex != mcpIndex[imcp]) continue;
          if(maxEP > 0) std::cout<<"hmmm\n";
          maxEP = pfp.EffPur;
          mpfpLoc = iit;
        } // iit
        MCP_EPTSum += TMeV * maxEP;
        ++MCP_PFP_Cnt;
        if(longMCP && maxEP > 0.8) ++nGoodLongMCP;
        if(mpfpLoc >= pfpLocs.size()) {
          mf::LogVerbatim myprt("TC");
          myprt<<"BadP: MCParticle "<<mcpIndex[imcp]<<" PDG "<<mcp.PdgCode()<<" T "<<(int)TMeV<<" MeV not reconstructed.";
          myprt<<" "<<mcpHits.size()<<" mcp Hits";
          myprt<<" events processed "<<evt.eventsProcessed;
          std::vector<std::pair<int, int>> tlist;
          // find all hits that are matched to this mcp, find hit -> Tj assns and put the Tj IDs into
          // the vector for printing
          for(auto& slc : slices) {
            for(auto& sht : slc.slHits) {
              if(sht.InTraj <= 0) continue;
              auto ahi = sht.allHitsIndex;
              if(evt.allHitsMCPIndex[ahi] != mcpIndex[imcp]) continue;
              // We found a hit that is matched to this MCParticle and is used in a Tj. Add it to the list
              // if it isn't already in it
              unsigned short ii = 0;
              for(ii = 0; ii < tlist.size(); ++ii) {
                if(tlist[ii].first == slc.ID && tlist[ii].second == sht.InTraj) break;
              } // ii
              if(ii == tlist.size()) tlist.push_back(std::make_pair(slc.ID, sht.InTraj));
            } // sht
          } // slc
          // print a list of trajectories
          myprt<<" Hits in";
          for(auto tid : tlist) myprt<<" "<<tid.first<<":T"<<tid.second;
          continue;
        } // no match
        if(maxEP < 0.8) {
          auto& pfp = slices[pfpLocs[mpfpLoc].first].pfps[pfpLocs[mpfpLoc].second];
          mf::LogVerbatim myprt("TC");
          myprt<<"BadP"<<pfp.ID<<" -> mcp "<<pfp.mcpIndex<<" T "<<(int)TMeV;
          myprt<<" PDG "<<pdg;
          myprt<<" EP "<<std::setprecision(2)<<pfp.EffPur;
          for(auto tid : pfp.TjIDs) myprt<<" T"<<tid;
        } // Poor EP
      } // imcp
      
      
      // debug primary electron reconstruction
/*
      if(tcc.modes[kStudy3] && !mcpIndex.empty() && mcpHits[0].size() > 10) {
        short TMeV = 1000 * (mcpIndex[0]->E() - mcpIndex[0]->Mass());
        std::cout<<"Study3: Find Tjs matched to primary w PDGCode "<<mcpList[0]->PdgCode()<<". T = "<<TMeV<<"\n";
        std::array<bool, 3> inPln {{false}};
        for(auto& slc : slices) {
          for(auto& tj : slc.tjs) {
            if(tj.AlgMod[kKilled]) continue;
            if(tj.mcpListIndex != 0) continue;
            unsigned short plane = DecodeCTP(tj.CTP).Plane;
            inPln[plane] = true;
            std::cout<<"T"<<tj.UID<<" start ";
            auto& tp = tj.Pts[tj.EndPt[0]];
            std::cout<<PrintPos(slc, tp);
            std::cout<<" PDGCode "<<tj.PDGCode;
            std::cout<<" len "<<tj.EndPt[1] - tj.EndPt[0] + 1;
            std::cout<<" MCSMom "<<tj.MCSMom;
            std::cout<<" ChgRMS "<<std::fixed<<std::setprecision(2)<<tj.ChgRMS;
            std::cout<<" BraggPeak? "<<tj.StopFlag[1][kBragg];
            std::cout<<" eLike "<<ElectronLikelihood(slc, tj);
            std::cout<<"\n";
          } // tj
          if(!slc.pfps.empty()) {
            auto& pfp = slc.pfps[0];
            std::cout<<"P"<<pfp.UID<<" PDGCode "<<pfp.PDGCode<<" dtrs";
            for(auto dtruid : pfp.DtrUIDs) std::cout<<" P"<<dtruid;
            std::cout<<"\n";
          } // pfps exist
        } // slc
        for(unsigned short plane = 0; plane < 3; ++plane) if(!inPln[plane]) std::cout<<"No match in plane "<<plane<<"\n";
      } // kStudy2
*/
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
