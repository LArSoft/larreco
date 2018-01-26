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
  void TruthMatcher::MatchTrueHits()
  {
    // Matches reco hits to MC true tracks and puts the association into
    // TCHit TruTrkID. This code is almost identical to the first part of MatchTruth.
    
    tjs.MCPartList.clear();

    if(tjs.MatchTruth[0] < 0) return;
    if(tjs.fHits.empty()) return;
    
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
      if(tjs.MatchTruth[0] == 1) {
        // Look for beam neutrino or single particle
        if(theTruth->Origin() == simb::kBeamNeutrino) {
          sourcePtclTrackID = trackID;
          sourceOrigin = simb::kBeamNeutrino;
          if(tjs.MatchTruth[1] > 0) {
            Vector3_t dir;
            dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
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
          if(tjs.MatchTruth[1] > 0) {
            Vector3_t dir;
            dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
            SetMag(dir, 1);
            std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode();
            std::cout<<" Vx "<<std::fixed<<std::setprecision(1)<<mcp->Vx()<<" Vy "<<mcp->Vy()<<" Vz "<<mcp->Vz();
            std::cout<<" dir "<<std::setprecision(3)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            break;
          }
        } // single particle
      } else if(tjs.MatchTruth[0] == 2 && theTruth->Origin() == simb::kCosmicRay) {
        sourcePtclTrackID = trackID;
        sourceOrigin = simb::kCosmicRay;
      }
    } // ipart
    
    if(sourcePtclTrackID == -1) return;
      
    // flag MCParticles to select for measuring performance. 
    std::vector<bool> select(plist.size(), false);
    // ensure that we get all of the primary particles
    unsigned int indx = 0;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* mcp = (*ipart).second;
      // primaries come first
      if(mcp->Process() != "primary") break;
      select[indx] = true;
      ++indx;
    }
    // make a temp vector of hit -> geant trackID
    std::vector<int> gtid(tjs.fHits.size(), 0);
    // find hits that match to the source particle
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      std::vector<sim::TrackIDE> ides;
      auto& tcHit = tjs.fHits[iht];
      geo::PlaneID planeID = geo::PlaneID(tcHit.ArtPtr->WireID().Cryostat, tcHit.ArtPtr->WireID().TPC, tcHit.ArtPtr->WireID().Plane);
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)tcHit.ArtPtr->WireID().Plane, (int)tcHit.ArtPtr->WireID().Wire, (int)tcHit.ArtPtr->WireID().TPC, (int)tcHit.ArtPtr->WireID().Cryostat);
      auto rhit = recob::Hit(channel,
                             tcHit.StartTick, tcHit.EndTick,
                             tcHit.PeakTime, tcHit.SigmaPeakTime,
                             tcHit.RMS,
                             tcHit.PeakAmplitude, tcHit.SigmaPeakAmp,
                             tcHit.Integral, tcHit.Integral, tcHit.SigmaIntegral,
                             tcHit.Multiplicity, tcHit.LocalIndex,
                             tcHit.GoodnessOfFit, tcHit.NDOF,
                             tjs.geom->View(channel),
                             tjs.geom->SignalType(planeID),
                             tcHit.ArtPtr->WireID());
      try {
        ides = bt_serv->HitToTrackIDEs(rhit);
      }
      catch(...) {}
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

    // save the selected MCParticles in tjs
    indx = 0;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      if(select[indx]) {
        simb::MCParticle* mcp = (*ipart).second;
        tjs.MCPartList.push_back(mcp);
      }
      ++indx;
    } // ipart
    
    if(tjs.MCPartList.size() > UINT_MAX) {
      std::cout<<"MatchTrueHits: Crazy large number of MCParticles "<<tjs.MCPartList.size()<<". Ignoring this event\n";
      tjs.MCPartList.clear();
      return;
    }
    
    // define MCPartListIndex for the hits
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(gtid[iht] == 0) continue;
      auto& hit = tjs.fHits[iht];
      for(unsigned int indx = 0; indx < tjs.MCPartList.size(); ++indx) {
        auto& mcp = tjs.MCPartList[indx];
        if(mcp->TrackId() != gtid[iht]) continue;
        hit.MCPartListIndex = indx;
      } // indx
    } // iht

  } // MatchTrueHits

  //////////////////////////////////////////
  void TruthMatcher::MatchTruth(const HistStuff& hist, bool fStudyMode)
  {
    // The hits have already been matched to the truth in MatchTrueHits. Here we match reconstructed objects
    // to the truth-matched hits to measure performance    
    
    if(tjs.MatchTruth[0] < 0) return;
    if(tjs.MCPartList.empty()) return;
    
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // list of all true particles
    sim::ParticleList const& plist = pi_serv->ParticleList();
    if(plist.empty()) return;

    Point3_t PrimVtx;
    // these are only applicable to neutrinos
    bool neutrinoVxInFiducialVolume = false;
    bool neutrinoVxReconstructable = false;
    bool neutrinoVxReconstructed = false;
    bool neutrinoVxCorrect = false;

    // Look for the MC truth process that should be considered (beam neutrino,
    // single particle, cosmic rays), then form a list of selected MCParticles 
    // that will be used to measure performance
    std::vector<unsigned int> mcpSelect;
    int sourceParticleIndex = -1;
    geo::TPCID inTPCID = tjs.TPCID;
    for(unsigned int part = 0; part < tjs.MCPartList.size(); ++part) {
      auto& mcp = tjs.MCPartList[part];
      int trackID = mcp->TrackId();
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
      bool originNeutrino = theTruth->Origin() == simb::kBeamNeutrino;
      int pdg = abs(mcp->PdgCode());
      bool isNeutrino = (pdg == 12) || (pdg == 14) || (pdg == 16);
      // ignore neutrinos that aren't primary
      if(isNeutrino && sourceParticleIndex >= 0) continue;
      if(isNeutrino && originNeutrino && sourceParticleIndex < 0) {
        sourceParticleIndex = part;
        continue;
      } // primary neutrino
      bool selectMCP = ((tjs.MatchTruth[0] == 1 && originNeutrino) || 
                        (tjs.MatchTruth[0] == 1 && theTruth->Origin() == simb::kSingleParticle) ||
                        (tjs.MatchTruth[0] == 2 && theTruth->Origin() == simb::kCosmicRay));
      if(!selectMCP) continue;
      // First occurrence of something the user wants to track other than neutrino origin handled above.
      // Note that the sourceParticleIndex isn't useful if there are multiple single particles or cosmic rays.
      if(sourceParticleIndex < 0) {
        sourceParticleIndex = part;
        PrimVtx[0] = mcp->Vx();
        PrimVtx[1] = mcp->Vy();
        PrimVtx[2] = mcp->Vz();
        // select the primary whether it is charged or not
        mcpSelect.push_back(part);
        // Determine which TPC this is in. Ignore this event if the primary start is
        // outside the fiducial volume
        if(!InsideTPC(tjs, PrimVtx, inTPCID)) return;
        // print out?
        if(tjs.MatchTruth[1] > 0) {
          Vector3_t dir;
          dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
          SetMag(dir, 1);
          std::cout<<"Found primary MCParticle "<<trackID<<" PDG code "<<mcp->PdgCode();
          std::cout<<" start";
          for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) std::cout<<" "<<std::setprecision(1)<<PrimVtx[ixyz];
          std::cout<<" dir";
          for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) std::cout<<" "<<std::setprecision(1)<<dir[ixyz];
          std::cout<<"\n";
        } // prt
        continue;
      } // first useMCP
      // only select charged particles
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      // cut on kinetic energy
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      if(TMeV < 10) continue;
      // require that it can be reconstructed in 3D
      if(!CanReconstruct(part, 3, inTPCID)) continue;
      mcpSelect.push_back(part);
    } // ipart
    
    if(mcpSelect.empty()) return;
    tjs.SelectEvent = true;
//    mf::LogVerbatim("TC")<<"SelectEvent "<<tjs.Run<<" "<<tjs.SubRun<<" "<<tjs.Event;

    if(neutrinoVxInFiducialVolume) ++TruVxCounts[0];
    
    // Form a list of mother-daughter pairs that should be considered as a single particle
    std::vector<std::pair<unsigned int, unsigned int>> moda;
    for(unsigned int part = 0; part < mcpSelect.size(); ++part) {
      auto& mcp = tjs.MCPartList[mcpSelect[part]];
      if(mcp->NumberDaughters() == 0) continue;
      unsigned int ndtr = 0;
      unsigned int dtrSelIndex = 0;
      for(int idtr = 0; idtr < mcp->NumberDaughters(); ++idtr) {
        int dtrTrackId = mcp->Daughter(idtr);
        // ignore it if it's not in the list
        bool ignore = true;
        for(unsigned short dpart = part + 1; dpart < mcpSelect.size(); ++dpart) {
          auto& dmcp = tjs.MCPartList[mcpSelect[dpart]];
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
      auto& dtr = tjs.MCPartList[mcpSelect[dtrSelIndex]];
      if(dtr->PdgCode() != mcp->PdgCode()) continue;
      if(tjs.MatchTruth[1] > 1) mf::LogVerbatim("TC")<<"Daughter MCP "<<dtrSelIndex<<" -> Mother MCP "<<part;
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
      for(auto& hit : tjs.fHits) {
        // see if this is a daughter
        for(auto& md : moda) if(md.second == hit.MCPartListIndex) hit.MCPartListIndex = md.first;
      } // hit
    } // !moda.empty

    // Match tjs to MC particles. Declare it a match to the MC particle that has the most
    // hits in the tj
    std::vector<std::vector<int>> tjid(tjs.NumPlanes);
    std::vector<std::vector<unsigned short>> ntht(tjs.NumPlanes);
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
      tjid[plane].resize(tjs.MCPartList.size());
      ntht[plane].resize(tjs.MCPartList.size());
    }
    
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.Cryostat != inTPCID.Cryostat) continue;
      if(planeID.TPC != inTPCID.TPC) continue;
      tj.MCPartListIndex = UINT_MAX;
      std::vector<unsigned int> mcpIndex, cnt;
      auto tjhits = PutTrajHitsInVector(tj, kUsedHits);
      for(auto iht : tjhits) {
        auto& hit = tjs.fHits[iht];
        if(hit.MCPartListIndex > tjs.MCPartList.size() - 1) continue;
        // ignore Tj hits that aren't in mcpSelect
        if(std::find(mcpSelect.begin(), mcpSelect.end(), hit.MCPartListIndex) == mcpSelect.end()) continue;
        unsigned int indx = 0;
        for(indx = 0; indx < mcpIndex.size(); ++indx) if(hit.MCPartListIndex == mcpIndex[indx]) break;
        if(indx == mcpIndex.size()) {
          mcpIndex.push_back(hit.MCPartListIndex);
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
          auto& mtj = tjs.allTraj[tjid[plane][tmpIndex] - 1];
          mtj.EffPur = 0;
          mtj.MCPartListIndex = UINT_MAX;
          // use the new one
          tj.EffPur = ep;
          tj.MCPartListIndex = tmpIndex;
          tjid[plane][tmpIndex] = tj.ID;
          ntht[plane][tmpIndex] = maxcnt;
        } else {
          // use the old one
          continue;
        }
      } else {
        // no previous match
        tj.EffPur = eff * pur;
        tj.MCPartListIndex = tmpIndex;
        tjid[plane][tmpIndex] = tj.ID;
        ntht[plane][tmpIndex] = maxcnt;
      }
    } // tj
    
    // count the number of primary tracks that can be reconstructed
    unsigned short nTruPrimary = 0;
    // vector of reconstructable primary particles
    std::vector<unsigned int> primMCPs;
    for(auto mcpIndex : mcpSelect) {
      auto& mcp = tjs.MCPartList[mcpIndex];
      if(mcp->Mother() != 0) continue;
      ++nTruPrimary;
      if(CanReconstruct(mcpIndex, 3, inTPCID)) primMCPs.push_back(mcpIndex);
    } // mcpIndex
    
    neutrinoVxReconstructable = (neutrinoVxInFiducialVolume && primMCPs.size() > 1);
    
    if(neutrinoVxReconstructable) {
      ++TruVxCounts[1];
      // Find the closest reconstructed vertex to the true vertex
      float closest = 1;
      unsigned short imTheOne = 0;
      for(auto& aVtx3 : tjs.vtx3) {
        float dx = aVtx3.X - PrimVtx[0];
        float dy = aVtx3.Y - PrimVtx[1];
        float dz = aVtx3.Z - PrimVtx[2];
        hist.fNuVtx_dx->Fill(dx);
        hist.fNuVtx_dy->Fill(dy);
        hist.fNuVtx_dz->Fill(dz);
        float sep = dx * dx + dy * dy + dz * dz;
        if(sep < closest) {
          closest = sep;
          imTheOne = aVtx3.ID;
        }
      } // aVtx3
      if(imTheOne > 0) {
        neutrinoVxReconstructed = true;
        ++TruVxCounts[2];
        auto& vx3 = tjs.vtx3[imTheOne - 1];
        hist.fNuVx3Score->Fill(vx3.Score);
        // Histogram the score of 2D vertices
        for(auto vx2id : vx3.Vx2ID) {
          if(vx2id == 0) continue;
          auto& vx2 = tjs.vtx[vx2id - 1];
          hist.fNuVx2Score->Fill(vx2.Score);
        } // vx2id
        // histogram the relative score of other vertices
        float maxScore = 0;
        for(auto& ovx3 : tjs.vtx3) {
          if(ovx3.ID == 0) continue;
          if(ovx3.Score > maxScore) maxScore = ovx3.Score;
          if(ovx3.ID == vx3.ID) continue;
          float dScore = ovx3.Score - vx3.Score;
          hist.fNuVx3ScoreDiff->Fill(dScore);
        } // ovx3
        neutrinoVxCorrect = (maxScore == vx3.Score);
        if(neutrinoVxCorrect) ++TruVxCounts[3];
        // find the most common Topo of the 2D vertices that were matched to the
        // primary vertex. This might be a useful to tag neutrino interaction vertices
        float vx3Topo = Vx3Topo(tjs, vx3);
        hist.fVxTopoMat->Fill(vx3Topo);
      }
    } // neutrinoVxInFiducialVolume
    
    if(neutrinoVxInFiducialVolume && neutrinoVxReconstructable && !neutrinoVxCorrect) {
      // More than one reconstructable primaries so there must a reconstructable neutrino vertex
      mf::LogVerbatim("TC")<<"BadVtx Reconstructed? "<<neutrinoVxReconstructed<<" and not correct. events processed "<<tjs.EventsProcessed;
    }

    if(tjs.MatchTruth[1] > 0) {
      // print out
      mf::LogVerbatim myprt("TC");
      myprt<<"Number of primary particles "<<nTruPrimary<<" Vtx";
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) myprt<<" "<<std::fixed<<std::setprecision(1)<<PrimVtx[ixyz];
      myprt<<" Reconstructable? "<<neutrinoVxReconstructable<<" Reconstructed? "<<neutrinoVxReconstructed<<" Correct? "<<neutrinoVxCorrect<<"\n";
      myprt<<"mcpIndex   PDG  momIndex    KE _________Dir___________       Process         TrajectoryExtentInPlane_nTruHits \n";
      for(auto mcpIndex : mcpSelect) {
        auto& mcp = tjs.MCPartList[mcpIndex];
        // Kinetic energy in MeV
        int TMeV = 1000 * (mcp->E() - mcp->Mass());
        myprt<<std::setw(8)<<mcpIndex;
        myprt<<std::setw(6)<<mcp->PdgCode();
        unsigned int momIndex = 0;
        if(mcp->Mother() > 0) {
          for(unsigned int part = 0; part < tjs.MCPartList.size(); ++part) {
            if(tjs.MCPartList[part]->TrackId() == mcp->Mother()) momIndex = part;
          }
        } // Mother() > 0
        myprt<<std::setw(10)<<momIndex;
        myprt<<std::setw(6)<<TMeV;
        Vector3_t dir;
        dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
        SetMag(dir, 1);
        myprt<<std::setprecision(3);
        for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) myprt<<std::setw(8)<<std::setprecision(3)<<dir[ixyz];
        myprt<<std::setw(22)<<mcp->Process();
        // print the extent of the particle in each TPC plane
        for(const geo::TPCID& tpcid : tjs.geom->IterateTPCIDs()) {
          geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
          for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
            CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
            auto mcpHits = PutMCPHitsInVector(mcpIndex, inCTP);
            if(mcpHits.empty()) continue;
//            if(mcpHits.size() < 3) continue;
            myprt<<" "<<PrintHitShort(tjs.fHits[mcpHits[0]])<<"-"<<PrintHitShort(tjs.fHits[mcpHits[mcpHits.size() - 1]]);
            myprt<<"_"<<mcpHits.size();
          } // plane
        } // tpcid
        myprt<<"\n";
      } // ipart
    } // tjs.MatchTruth[1] > 0
    
    // Match Tjs and PFParticles and accumulate statistics
    MatchAndSum(hist, mcpSelect, inTPCID);
/* This needs work...
    // match 2D vertices (crudely)
    for(auto& tj : tjs.allTraj) {
      // obsolete vertex
      if(tj.AlgMod[kKilled]) continue;
      // require a truth match
      if(tj.MCPartListIndex == UINT_MAX) continue;
      if (tj.MCPartListIndex < tjs.MCPartList.size()) { 
        // ignore electrons unless it is a primary electron    
        auto& mcp = tjs.MCPartList[tj.MCPartListIndex];
        int pdg = abs(mcp->PdgCode());
        if(pdg == 11 && mcp->Mother() != 0) continue;
      }
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == 0) continue;
        VtxStore& vx2 = tjs.vtx[tj.VtxID[end]-1];
        vx2.Stat[kVtxTruMatch] = true;
      } // end
    } // tj
*/
    for(auto& vx2 : tjs.vtx) if(vx2.ID > 0) ++RecoVx2Count;
    
    if (fStudyMode) {
      // nomatch
      for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
        auto& pfp = tjs.pfps[ipfp];
        if(pfp.ID == 0) continue;
        if(pfp.MCPartListIndex == UINT_MAX) {
          // unmatched PFPs
          hist.hasfit_nomatch->Fill(pfp.Track.ID()>=0);
          if (pfp.Track.ID()<0) continue;
          hist.nvalidpoints_nomatch->Fill(pfp.Track.CountValidPoints());
          hist.nrejectpoints_nomatch->Fill(pfp.Track.NPoints()-pfp.Track.CountValidPoints());
          hist.fracreject_nomatch->Fill( float(pfp.Track.NPoints()-pfp.Track.CountValidPoints())/float(pfp.Track.NPoints()) );
          if (pfp.Track.CountValidPoints()>1) {
            hist.nchi2_nomatch->Fill( pfp.Track.Chi2PerNdof() );
            auto cv = pfp.Track.VertexCovarianceLocal5D();
            hist.covtrace_nomatch->Fill( cv(0,0)+cv(1,1)+cv(2,2)+cv(3,3) );
          }
        } else {
          // matched PFPs
          auto& mcp = tjs.MCPartList[pfp.MCPartListIndex];
          if (abs(mcp->PdgCode())==11) continue;
          bool wrongid = ( std::abs(mcp->PdgCode())!=std::abs(pfp.Track.ParticleId()) );
          //
          hist.hasfit_match->Fill(pfp.Track.ID()>=0);
          if (wrongid) {
            hist.hasfit_wrongid->Fill(pfp.Track.ID()>=0);
          } else {
            hist.hasfit_okid->Fill(pfp.Track.ID()>=0);
          }
          //
          if (pfp.Track.ID()<0) continue;
          // std::cout << "pfp #" << ipfp << " pdg=" << pfp.PDGCode << " tj vtx=" << recob::tracking::Point_t(pfp.XYZ[0][0],pfp.XYZ[0][1],pfp.XYZ[0][2]) << " dir=" << recob::tracking::Vector_t(pfp.Dir[0][0],pfp.Dir[0][1],pfp.Dir[0][2]) << " fit vtx=" << pfp.Track.Start() << " dir=" << pfp.Track.StartDirection() << " match part #" << ipart << " vtx=" << recob::tracking::Point_t(mcp->Vx(), mcp->Vy(), mcp->Vz()) << " dir=" << recob::tracking::Vector_t(mcp->Px()/mcp->P(), mcp->Py()/mcp->P(), mcp->Pz()/mcp->P()) << " mom=" << mcp->P() << std::endl;
          hist.nvalidpoints_match->Fill(pfp.Track.CountValidPoints());
          hist.nrejectpoints_match->Fill(pfp.Track.NPoints()-pfp.Track.CountValidPoints());
          hist.fracreject_match->Fill( float(pfp.Track.NPoints()-pfp.Track.CountValidPoints())/float(pfp.Track.NPoints()) );
          if (pfp.Track.CountValidPoints()>1) {
            //
            detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
            double g4Ticks = detClocks->TPCG4Time2Tick(mcp->T())+tjs.detprop->GetXTicksOffset(0,0,0)-tjs.detprop->TriggerOffset();
            double xOffset = tjs.detprop->ConvertTicksToX(g4Ticks, 0, 0, 0);
            //
            trkf::TrackStatePropagator prop(1.0, 0.1, 10, 10., 0.01, false);
            recob::tracking::Plane mcplane(recob::tracking::Point_t(mcp->Vx()+xOffset, mcp->Vy(), mcp->Vz()),
                                           recob::tracking::Vector_t(mcp->Px()/mcp->P(), mcp->Py()/mcp->P(), mcp->Pz()/mcp->P()));
            recob::tracking::Plane tkplane(pfp.Track.Start(), pfp.Track.StartDirection());
            trkf::TrackState tkstart(pfp.Track.VertexParametersLocal5D(), pfp.Track.VertexCovarianceLocal5D(), tkplane, true, pfp.Track.ParticleId());
            //
            recob::tracking::Plane tjplane(recob::tracking::Point_t(pfp.XYZ[0][0],pfp.XYZ[0][1],pfp.XYZ[0][2]),
                                           recob::tracking::Vector_t(pfp.Dir[0][0],pfp.Dir[0][1],pfp.Dir[0][2]));
            auto tjpar5 = tjplane.Global6DToLocal5DParameters({pfp.XYZ[0][0],pfp.XYZ[0][1],pfp.XYZ[0][2],pfp.Dir[0][0],pfp.Dir[0][1],pfp.Dir[0][2]});
            trkf::TrackState tjstart(tjpar5, pfp.Track.VertexCovarianceLocal5D(), tjplane, true, pfp.Track.ParticleId());
            //
            bool propok = true;
            trkf::TrackState tkatmc = prop.propagateToPlane(propok, tkstart, mcplane, true, true, trkf::TrackStatePropagator::UNKNOWN);
            trkf::TrackState tjatmc = prop.propagateToPlane(propok, tjstart, mcplane, true, true, trkf::TrackStatePropagator::UNKNOWN);
            //
            if (!propok) continue;
            //
            auto recmom = tkatmc.momentum().R();
            auto c = tkatmc.plane().Local5DToGlobal6DCovariance(tkatmc.covariance(),false,tkatmc.momentum());//consider direction not momentum
            auto cv = pfp.Track.VertexCovarianceLocal5D();
            //
            hist.dXkf_match->Fill(tkatmc.position().X()-mcp->Vx()-xOffset);
            hist.dXtc_match->Fill(tjatmc.position().X()-mcp->Vx()-xOffset);
            hist.dYkf_match->Fill(tkatmc.position().Y()-mcp->Vy());
            hist.dYtc_match->Fill(tjatmc.position().Y()-mcp->Vy());
            hist.dZkf_match->Fill(tkatmc.position().Z()-mcp->Vz());
            hist.dZtc_match->Fill(tjatmc.position().Z()-mcp->Vz());
            hist.dUXkf_match->Fill(tkatmc.momentum().X()/recmom-mcp->Px()/mcp->P());
            hist.dUXtc_match->Fill(tjatmc.momentum().Unit().X()-mcp->Px()/mcp->P());
            hist.dUYkf_match->Fill(tkatmc.momentum().Y()/recmom-mcp->Py()/mcp->P());
            hist.dUYtc_match->Fill(tjatmc.momentum().Unit().Y()-mcp->Py()/mcp->P());
            hist.dUZkf_match->Fill(tkatmc.momentum().Z()/recmom-mcp->Pz()/mcp->P());
            hist.dUZtc_match->Fill(tjatmc.momentum().Unit().Z()-mcp->Pz()/mcp->P());
            hist.dXpull_match->Fill( (tkatmc.position().X()-mcp->Vx()-xOffset)/ sqrt(c(0,0)) );
            hist.dYpull_match->Fill( (tkatmc.position().Y()-mcp->Vy())/ sqrt(c(1,1)) );
            hist.dZpull_match->Fill( (tkatmc.position().Z()-mcp->Vz())/ sqrt(c(2,2)) );
            hist.dUXpull_match->Fill( (tkatmc.momentum().X()/recmom-mcp->Px()/mcp->P())/ sqrt(c(3,3)) );
            hist.dUYpull_match->Fill( (tkatmc.momentum().Y()/recmom-mcp->Py()/mcp->P())/ sqrt(c(4,4)) );
            hist.dUZpull_match->Fill( (tkatmc.momentum().Z()/recmom-mcp->Pz()/mcp->P())/ sqrt(c(5,5)) );
            hist.nchi2_match->Fill( pfp.Track.Chi2PerNdof() );
            hist.covtrace_match->Fill( cv(0,0)+cv(1,1)+cv(2,2)+cv(3,3) );
            if (wrongid) {
              int wid = 0;
              if (std::abs(mcp->PdgCode())==11) wid=1;
              if (std::abs(mcp->PdgCode())==13) wid=2;
              if (std::abs(mcp->PdgCode())==211) wid=3;
              if (std::abs(mcp->PdgCode())==2212) wid=4;
              hist.pdgid_wrongid->Fill( wid );
              hist.nchi2_wrongid->Fill( pfp.Track.Chi2PerNdof() );
              hist.dXkf_wrongid->Fill(tkatmc.position().X()-mcp->Vx()-xOffset);
              hist.dXtc_wrongid->Fill(tjatmc.position().X()-mcp->Vx()-xOffset);
              hist.dYkf_wrongid->Fill(tkatmc.position().Y()-mcp->Vy());
              hist.dYtc_wrongid->Fill(tjatmc.position().Y()-mcp->Vy());
              hist.dZkf_wrongid->Fill(tkatmc.position().Z()-mcp->Vz());
              hist.dZtc_wrongid->Fill(tjatmc.position().Z()-mcp->Vz());
              hist.dUXkf_wrongid->Fill(tkatmc.momentum().X()/recmom-mcp->Px()/mcp->P());
              hist.dUXtc_wrongid->Fill(tjatmc.momentum().Unit().X()-mcp->Px()/mcp->P());
              hist.dUYkf_wrongid->Fill(tkatmc.momentum().Y()/recmom-mcp->Py()/mcp->P());
              hist.dUYtc_wrongid->Fill(tjatmc.momentum().Unit().Y()-mcp->Py()/mcp->P());
              hist.dUZkf_wrongid->Fill(tkatmc.momentum().Z()/recmom-mcp->Pz()/mcp->P());
              hist.dUZtc_wrongid->Fill(tjatmc.momentum().Unit().Z()-mcp->Pz()/mcp->P());
            } else {
              hist.dXkf_okid->Fill(tkatmc.position().X()-mcp->Vx()-xOffset);
              hist.dXtc_okid->Fill(tjatmc.position().X()-mcp->Vx()-xOffset);
              hist.dYkf_okid->Fill(tkatmc.position().Y()-mcp->Vy());
              hist.dYtc_okid->Fill(tjatmc.position().Y()-mcp->Vy());
              hist.dZkf_okid->Fill(tkatmc.position().Z()-mcp->Vz());
              hist.dZtc_okid->Fill(tjatmc.position().Z()-mcp->Vz());
              hist.dUXkf_okid->Fill(tkatmc.momentum().X()/recmom-mcp->Px()/mcp->P());
              hist.dUXtc_okid->Fill(tjatmc.momentum().Unit().X()-mcp->Px()/mcp->P());
              hist.dUYkf_okid->Fill(tkatmc.momentum().Y()/recmom-mcp->Py()/mcp->P());
              hist.dUYtc_okid->Fill(tjatmc.momentum().Unit().Y()-mcp->Py()/mcp->P());
              hist.dUZkf_okid->Fill(tkatmc.momentum().Z()/recmom-mcp->Pz()/mcp->P());
              hist.dUZtc_okid->Fill(tjatmc.momentum().Unit().Z()-mcp->Pz()/mcp->P());
              hist.dXpull_okid->Fill( (tkatmc.position().X()-mcp->Vx()-xOffset)/ sqrt(c(0,0)) );
              hist.dYpull_okid->Fill( (tkatmc.position().Y()-mcp->Vy())/ sqrt(c(1,1)) );
              hist.dZpull_okid->Fill( (tkatmc.position().Z()-mcp->Vz())/ sqrt(c(2,2)) );
              hist.dUXpull_okid->Fill( (tkatmc.momentum().X()/recmom-mcp->Px()/mcp->P())/ sqrt(c(3,3)) );
              hist.dUYpull_okid->Fill( (tkatmc.momentum().Y()/recmom-mcp->Py()/mcp->P())/ sqrt(c(4,4)) );
              hist.dUZpull_okid->Fill( (tkatmc.momentum().Z()/recmom-mcp->Pz()/mcp->P())/ sqrt(c(5,5)) );
              hist.nchi2_okid->Fill( pfp.Track.Chi2PerNdof() );
            }
            //
            // if (std::abs(tkatmc.position().X()-mcp->Vx())>10. && std::abs(tjatmc.position().X()-mcp->Vx())<1.) {
            //   std::cout << "DEBUG ME" << std::endl;
            // }
            //
          }
        } // matched PFP
      } // ipfp
    } //fStudyMode
    
    
    // histogram reconstructed PDG code vs true PDG code
    std::array<int, 5> recoCodeList = {0, 11, 13, 211, 2212};
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      // ignore showers
      if(pfp.PDGCode == 1111) continue;
      // require match to MC
      if(pfp.MCPartListIndex == UINT_MAX) continue;
      short truIndex = PDGCodeIndex(tjs, tjs.MCPartList[pfp.MCPartListIndex]->PdgCode());
      if(truIndex == SHRT_MAX) continue;
      short recIndex = 0;
      for(recIndex = 0; recIndex < 5; ++recIndex) if(pfp.PDGCode == recoCodeList[recIndex]) break;
      if(recIndex == 5) {
        std::cout<<"MatchTruth: Found an unknown PDGCode "<<pfp.PDGCode<<" in PFParticle "<<pfp.ID<<"\n";
        continue;
      }
      hist.PDGCode_reco_true->Fill((float)truIndex, (float)recIndex);
    } // pfp


  } // MatchTruth
  
  ////////////////////////////////////////////////
  void TruthMatcher::MatchAndSum(const HistStuff& hist, const std::vector<unsigned int>& mcpSelect, const geo::TPCID& inTPCID)
  {
    // Match Tjs and PFParticles and accumulate performance statistics
    if(mcpSelect.empty()) return;
    
    unsigned int tpc = inTPCID.TPC;
    unsigned int cstat = inTPCID.Cryostat;
    
    // get the hits associated with all MCParticles in mcpSelect
    std::vector<std::vector<unsigned int>> mcpHits(mcpSelect.size());
    for(unsigned short isel = 0; isel < mcpSelect.size(); ++isel) {
      unsigned int mcpIndex = mcpSelect[isel];
      for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
        if(tjs.fHits[iht].MCPartListIndex != mcpIndex) continue;
        if(tjs.fHits[iht].ArtPtr->WireID().TPC != tpc) continue;
        if(tjs.fHits[iht].ArtPtr->WireID().Cryostat != cstat) continue;
        mcpHits[isel].push_back(iht);
      } // iht
    } // mcpIndex
    
    // get the hits in all Tjs in this TPCID
    std::vector<std::vector<unsigned int>> tjHits(tjs.allTraj.size());
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      auto& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      if(DecodeCTP(tj.CTP).TPC != tpc) continue;
      tjHits[itj] = PutTrajHitsInVector(tj, kUsedHits);
      tj.MCPartListIndex = UINT_MAX;
    } // itj
    
    // match them up
    for(unsigned short isel = 0; isel < mcpSelect.size(); ++isel) {
      if(mcpHits[isel].empty()) continue;
      unsigned int mcpIndex = mcpSelect[isel];
      auto& mcp = tjs.MCPartList[mcpIndex];
      unsigned short pdgIndex = PDGCodeIndex(tjs, mcp->PdgCode());
      if(pdgIndex > 4) continue;
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        // put the mcpHits (which are already in this TPCID) in this plane in a vector
        std::vector<unsigned int> mcpPlnHits;
        for(auto iht : mcpHits[isel]) {
          auto& hit = tjs.fHits[iht];
          if(hit.ArtPtr->WireID().Plane == plane) mcpPlnHits.push_back(iht);
        } // iht
        // require 2 truth-matched hits
        if(mcpPlnHits.size() < 2) continue;
        if((float)mcpPlnHits.size() >= tjs.MatchTruth[3]) ++nLongInPln;
        TSums[pdgIndex] += TMeV;
        ++EPCnts[pdgIndex];
        CTP_t inCTP = EncodeCTP(cstat, tpc, plane);
        unsigned short mtj = USHRT_MAX;
        float maxEP = 0;
        for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
          // No hits in this TPC and plane
          if(tjHits[itj].empty()) continue;
          auto& tj = tjs.allTraj[itj];
          // wrong CTP?
          if(tj.CTP != inCTP) continue;
          // make a list of hits that are common
          auto shared = SetIntersection(mcpPlnHits, tjHits[itj]);
          if(shared.empty()) continue;
          float eff = (float)shared.size() / (float)mcpPlnHits.size();
          float pur = (float)shared.size() / (float)tjHits[itj].size();
          float ep = eff * pur;
          // Replace a previously made poorer match with a better one?
          if(tj.MCPartListIndex != UINT_MAX && ep > tj.EffPur) {
            tj.EffPur = ep;
            tj.MCPartListIndex = mcpIndex;
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
          if((float)mcpPlnHits.size() > tjs.MatchTruth[3]) {
            ++nBadEP;
            mf::LogVerbatim myprt("TC");
            myprt<<"pdgIndex "<<pdgIndex<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to isel "<<isel;
            myprt<<" nTrue hits "<<mcpPlnHits.size();
            myprt<<" extent "<<PrintHit(tjs.fHits[mcpPlnHits[0]])<<"-"<<PrintHit(tjs.fHits[mcpPlnHits[mcpPlnHits.size() - 1]]);
            myprt<<" events processed "<<tjs.EventsProcessed;
          }
          continue;
        } // match failed
        auto& tj = tjs.allTraj[mtj];
        // don't clobber a better match
        if(maxEP < tj.EffPur) continue;
        tj.EffPur = maxEP;
        tj.MCPartListIndex = mcpIndex;
        EPTSums[pdgIndex] += TMeV * tj.EffPur;
        hist.fEP_T[pdgIndex]->Fill(TMeV, tj.EffPur);
        if(tj.EffPur < tjs.MatchTruth[2] && (float)mcpPlnHits.size() >= tjs.MatchTruth[3]) {
          ++nBadEP;
          mf::LogVerbatim myprt("TC");
          myprt<<"pdgIndex "<<pdgIndex<<" BadEP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
          myprt<<" mcpIndex "<<mcpIndex;
          myprt<<" TMeV "<<(int)TMeV<<" MCP hits "<<mcpPlnHits.size();
          myprt<<" extent "<<PrintHit(tjs.fHits[mcpPlnHits[0]])<<"-"<<PrintHit(tjs.fHits[mcpPlnHits[mcpPlnHits.size() - 1]]);
          myprt<<" events processed "<<tjs.EventsProcessed;
        }
        // badep
      } // plane
    } // isel

    // Calculate PFParticle efficiency and purity
    if(tjs.pfps.empty()) return;
    
    // get the hits in all pfparticles in this TPCID
    std::vector<std::vector<unsigned int>> pfpHits(tjs.pfps.size());
    for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
      auto& pfp = tjs.pfps[ipfp];
      if(pfp.ID == 0) continue;
      // in the right TPCID?
      if(pfp.TPCID != inTPCID) continue;
      // ignore the neutrino PFParticle
      if(pfp.PDGCode == 14 || pfp.PDGCode == 12) continue;
      pfp.MCPartListIndex = UINT_MAX;
      for(auto& tjid : pfp.TjIDs) {
        unsigned short itj = tjid - 1;
        pfpHits[ipfp].insert(pfpHits[ipfp].end(), tjHits[itj].begin(), tjHits[itj].end());
      } // tj
    } // ipfp
    
    // match them up
    for(unsigned short isel = 0; isel < mcpSelect.size(); ++isel) {
      unsigned short mpfp = USHRT_MAX;
      float maxEP = 0;
      unsigned int mcpIndex = mcpSelect[isel];
      if(mcpHits[isel].empty()) continue;
      auto& mcp = tjs.MCPartList[mcpIndex];
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      MCP_TSum += TMeV;
      // Performance reconstructing long muons, pions, kaons and protons
      unsigned short pdgIndex = PDGCodeIndex(tjs, mcp->PdgCode());
      bool longMCP = (pdgIndex > 0 && pdgIndex < 5 && (float)mcpHits[isel].size() >= 2 * tjs.MatchTruth[3]);
      if(longMCP) ++nLongMCP;
      for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
        auto& pfp = tjs.pfps[ipfp];
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
        if(pfp.MCPartListIndex != UINT_MAX && ep > pfp.EffPur) {
          pfp.EffPur = ep;
          pfp.MCPartListIndex = mcpIndex;
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
          for(auto& tj : tjs.allTraj) {
            if(tj.AlgMod[kKilled]) continue;
            if(tj.MCPartListIndex == mcpIndex) myprt<<" "<<tj.ID<<" EP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
          } // tj
          myprt<<" events processed "<<tjs.EventsProcessed;
        } // TMeV > 30
        continue;
      }
      auto& pfp = tjs.pfps[mpfp];
      if(maxEP < pfp.EffPur) continue;
      pfp.EffPur = maxEP;
      pfp.MCPartListIndex = mcpIndex;
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
    if(!tjs.SelectEvent) return;
    
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

  } // PrintResults
  
  ////////////////////////////////////////////////
  bool TruthMatcher::CanReconstruct(unsigned int mcpIndex, unsigned short nDimensions, const geo::TPCID& inTPCID)
  {
    // returns true if the MCParticle can be reconstructed in nDimensions
    if(mcpIndex > tjs.MCPartList.size() - 1) return false;
    if(nDimensions < 2 || nDimensions > 3) return false;
    
    std::vector<unsigned short> cntInPln(tjs.NumPlanes);
    for(auto& hit : tjs.fHits) {
      if(hit.ArtPtr->WireID().TPC != inTPCID.TPC) continue;
      if(hit.ArtPtr->WireID().Cryostat != inTPCID.Cryostat) continue;
      if(hit.MCPartListIndex == mcpIndex) ++cntInPln[hit.ArtPtr->WireID().Plane];
    } // hit
    unsigned short nPlnOK = 0;
    // Require at least 2 truth-matched hits in a plane
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) if(cntInPln[plane] > 1) ++nPlnOK;
    if(nPlnOK < nDimensions - 1) return false;
    return true;
  } // CanReconstruct
  
  
  ////////////////////////////////////////////////
  std::vector<unsigned int> TruthMatcher::PutMCPHitsInVector(unsigned int mcpIndex, CTP_t inCTP)
  {
    // put the hits matched to the MCParticle into a vector in the requested CTP
    std::vector<unsigned int> hitVec;
    if(mcpIndex > tjs.MCPartList.size() - 1) return hitVec;
    geo::PlaneID planeID = DecodeCTP(inCTP);
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      auto& hit = tjs.fHits[iht];
      if(hit.MCPartListIndex != mcpIndex) continue;
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
    
    if(MCParticleListIndex > tjs.MCPartList.size() - 1) return;
    
    const simb::MCParticle* mcp = tjs.MCPartList[MCParticleListIndex];
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    
    tp.Pos[0] = tjs.geom->WireCoordinate(mcp->Vy(), mcp->Vz(), planeID);
    tp.Pos[1] = tjs.detprop->ConvertXToTicks(mcp->Vx(), planeID) * tjs.UnitsPerTick;
    
    TVector3 dir;
    dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
    if(dir.Mag() == 0) return;
    dir.SetMag(1);
    TVector3 pos;
    pos[0] = mcp->Vx() + 100 * dir[0];
    pos[1] = mcp->Vy() + 100 * dir[1];
    pos[2] = mcp->Vz() + 100 * dir[2];
    // use HitPos as a work vector
    tp.HitPos[0] = tjs.geom->WireCoordinate(pos[1], pos[2], planeID);
    tp.HitPos[1] = tjs.detprop->ConvertXToTicks(pos[0], planeID) * tjs.UnitsPerTick;
    
    tp.Dir[0] = tp.HitPos[0] - tp.Pos[0];
    tp.Dir[1] = tp.HitPos[1] - tp.Pos[1];
    double norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    tp.Delta = norm / 100;
    
    // The Orth vectors are not unit normalized so we need to correct for this
    double w0 = tjs.geom->WireCoordinate(0, 0, planeID);
    // cosine-like component
    double cs = tjs.geom->WireCoordinate(1, 0, planeID) - w0;
    // sine-like component
    double sn = tjs.geom->WireCoordinate(0, 1, planeID) - w0;
    norm = sqrt(cs * cs + sn * sn);
    tp.Delta /= norm;
    
  } // MakeTruTrajPoint
  
  /////////////////////////////////////////
  unsigned short MCParticleListUtils::MCParticleStartTjID(unsigned int MCParticleListIndex, CTP_t inCTP)
  {
    // Finds the trajectory that has hits matched to the MC Particle and is the closest to the
    // MCParticle start vertex
    
    if(MCParticleListIndex > tjs.MCPartList.size() - 1) return 0;
    
    const simb::MCParticle* mcp = tjs.MCPartList[MCParticleListIndex];
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    TrajPoint truTp;
    truTp.Pos[0] = tjs.geom->WireCoordinate(mcp->Vy(), mcp->Vz(), planeID);
    truTp.Pos[1] = tjs.detprop->ConvertXToTicks(mcp->Vx(), planeID) * tjs.UnitsPerTick;
    
    unsigned short imTheOne = 0;
    unsigned short length = 5;
    unsigned short nTruHits;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled] && !tj.AlgMod[kInShower]) continue;
      if(tj.CTP != inCTP) continue;
      if(tj.Pts.size() < length) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        unsigned short ept = tj.EndPt[end];
        float sep2 = PosSep2(tj.Pts[ept].Pos, truTp.Pos);
        if(sep2 > 20) continue;
        // found a close trajectory point. See if this is the right one
        if(GetMCPartListIndex(tj, nTruHits) != MCParticleListIndex) continue;
        imTheOne = tj.ID;
        length = tj.Pts.size();
      } // end
    } // tj
    
    return imTheOne;
    
  } // MCParticleStartTj
  
  /////////////////////////////////////////
  unsigned int MCParticleListUtils::GetMCPartListIndex(const TrajPoint& tp)
  {
    // Returns the MCParticle index that best matches the hits used in the Tp
    if(tjs.MCPartList.empty()) return UINT_MAX;
    if(tp.Chg <= 0) return UINT_MAX;
    std::vector<unsigned int> mcpIndex;
    std::vector<unsigned short> mcpCnt;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      unsigned int mcpi = tjs.fHits[tp.Hits[ii]].MCPartListIndex;
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
  } // GetMCPartListIndex
  
  /////////////////////////////////////////
  unsigned int MCParticleListUtils::GetMCPartListIndex(const ShowerStruct& ss, unsigned short& nTruHits)
  {
    // Returns the index of the MCParticle that has the most number of matches
    // to the hits in this shower
    
    if(tjs.MCPartList.empty()) return UINT_MAX;
    if(ss.TjIDs.empty()) return UINT_MAX;
    
    std::vector<unsigned int> pListCnt(tjs.MCPartList.size());
    
    for(auto& tjid : ss.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjid - 1];
      for(auto& tp : tj.Pts) {
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          // ignore unmatched hits
          if(tjs.fHits[iht].MCPartListIndex > tjs.MCPartList.size() - 1) continue;
          ++pListCnt[tjs.fHits[iht].MCPartListIndex];
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
    
  } // GetMCPartListIndex
  
  /////////////////////////////////////////
  unsigned int MCParticleListUtils::GetMCPartListIndex(const Trajectory& tj, unsigned short& nTruHits)
  {
    // Returns the index of the MCParticle that has the most number of matches
    // to the hits in this trajectory
    
    if(tjs.MCPartList.empty()) return UINT_MAX;
    
    // Check all hits associated with this Tj
    std::vector<unsigned int> pListCnt(tjs.MCPartList.size());
    
    for(auto& tp : tj.Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        // ignore unmatched hits
        if(tjs.fHits[iht].MCPartListIndex > tjs.MCPartList.size() - 1) continue;
        ++pListCnt[tjs.fHits[iht].MCPartListIndex];
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
    
  } // GetMCPartListIndex

} // namespace tca
