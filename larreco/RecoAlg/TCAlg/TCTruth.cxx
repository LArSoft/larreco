#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"


#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "larsim/MCCheater/BackTracker.h"

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
  void TruthMatcher::MatchTrueHits(const HistStuff& hist)
  {
    // Matches reco hits to MC true tracks and puts the association into
    // TCHit TruTrkID. This code is almost identical to the first part of MatchTruth.
    
    tjs.MCPartList.clear();

    if(tjs.MatchTruth[0] < 0) return;
    if(tjs.fHits.empty()) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;

    // first find the source particle 
    simb::Origin_t sourceOrigin = simb::kUnknown;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* mcp = (*ipart).second;
      int trackID = mcp->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(tjs.MatchTruth[0] == 1) {
        // Look for beam neutrino or single particle
        if(theTruth->Origin() == simb::kBeamNeutrino) {
          fSourceParticleEnergy = 1000 * mcp->E();
          sourcePtclTrackID = trackID;
          sourceOrigin = simb::kBeamNeutrino;
          if(tjs.MatchTruth[1] > 2) std::cout<<"Found beam neutrino sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode()<<"\n";
          break;
        } // beam neutrino
        if(theTruth->Origin() == simb::kSingleParticle) {
          fSourceParticleEnergy = 1000 * mcp->E();
          sourcePtclTrackID = trackID;
          sourceOrigin = simb::kSingleParticle;
          if(tjs.MatchTruth[1] > 0) {
            TVector3 dir;
            dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
            dir.SetMag(1);
            std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode()<<" Vx "<<(int)mcp->Vx()<<" Vy "<<(int)mcp->Vy()<<" Vz "<<(int)mcp->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
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
    // make a temp vector of hit -> geant trackID
    std::vector<int> gtid(tjs.fHits.size(), 0);
    // find hits that match to the source particle
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      std::vector<sim::IDE> ides;
      auto& tcHit = tjs.fHits[iht];
      geo::PlaneID planeID = geo::PlaneID(tcHit.WireID.Cryostat, tcHit.WireID.TPC, tcHit.WireID.Plane);
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)tcHit.WireID.Plane, (int)tcHit.WireID.Wire, (int)tcHit.WireID.TPC, (int)tcHit.WireID.Cryostat);
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
                             tcHit.WireID);
      try {
        bt->HitToSimIDEs(rhit, ides);
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
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(hitTruTrkID);
      if(theTruth->Origin() != sourceOrigin) continue;
      unsigned short indx = 0;
      for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
        ++indx;
        simb::MCParticle* mcp = (*ipart).second;
        if(mcp->TrackId() != hitTruTrkID) continue;
        select[indx - 1] = true;
        gtid[iht] = mcp->TrackId();
        // find the eve particle and select it as well
        const simb::MCParticle* momMCP = bt->TrackIDToMotherParticle(mcp->TrackId());
        unsigned short mindx = 0;
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
    unsigned short indx = 0;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      if(select[indx]) {
        simb::MCParticle* mcp = (*ipart).second;
        tjs.MCPartList.push_back(mcp);
      }
      ++indx;
    } // ipart
    
    if(tjs.MCPartList.size() > USHRT_MAX) {
      std::cout<<"MatchTrueHits: Crazy large number of MCParticles "<<tjs.MCPartList.size()<<". Ignoring this event\n";
      tjs.MCPartList.clear();
      return;
    }
    
    // define MCPartListIndex for the hits
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(gtid[iht] == 0) continue;
      auto& hit = tjs.fHits[iht];
      for(unsigned short indx = 0; indx < tjs.MCPartList.size(); ++indx) {
        auto& mcp = tjs.MCPartList[indx];
        if(mcp->TrackId() != gtid[iht]) continue;
        hit.MCPartListIndex = indx;
      } // indx
    } // iht
      

  } // MatchTrueHits

/*
  //////////////////////////////////////////
  void TruthMatcher::MatchTrueHits(const HistStuff& hist)
  {
    // Matches reco hits to MC true tracks and puts the association into
    // TCHit TruTrkID. This code is almost identical to the first part of MatchTruth.
    
    if(tjs.MatchTruth[0] < 0) return;
    if(tjs.fHits.empty()) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    if(plist.size() > USHRT_MAX) {
      std::cout<<"MatchTrueHits: plist too large for this code\n";
      return;
    }
    
    tjs.MCPartList.clear();
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;
    
    simb::Origin_t sourceOrigin = simb::kUnknown;
    // partList is the vector of MC particles that we want to use
    tjs.MCPartList.reserve(plist.size());
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* mcp = (*ipart).second;
      int trackID = mcp->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(sourcePtclTrackID < 0) {
        if(tjs.MatchTruth[0] == 1) {
          // Look for beam neutrino or single particle
          if(theTruth->Origin() == simb::kBeamNeutrino) {
            fSourceParticleEnergy = 1000 * mcp->E();
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
            if(tjs.MatchTruth[1] > 2) std::cout<<"Found beam neutrino sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode()<<"\n";
          } // beam neutrino
          if(theTruth->Origin() == simb::kSingleParticle) {
            fSourceParticleEnergy = 1000 * mcp->E();
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
            if(tjs.MatchTruth[1] > 0) {
              TVector3 dir;
              dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
              dir.SetMag(1);
              std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode()<<" Vx "<<(int)mcp->Vx()<<" Vy "<<(int)mcp->Vy()<<" Vz "<<(int)mcp->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            }
          } // single particle
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            if(tjs.MatchTruth[1] > 2) std::cout<<"True vertex position "<<(int)mcp->Vx()<<" "<<(int)mcp->Vy()<<" "<<(int)mcp->Vz()<<"\n";
          } // sourceOrigin != simb::kUnknown
        } else {
          // look for cosmic rays
          if(theTruth->Origin() == simb::kCosmicRay) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kCosmicRay;
          }
        }
      }
      if(trackID == 2474363) std::cout<<"Found it "<<mcp->PdgCode()<<" origin "<<theTruth->Origin()<<" "<<1000 * mcp->E()<<"\n";
      // ignore anything that has the incorrect origin
      if(theTruth->Origin() != sourceOrigin) continue;
      // ignore processes that aren't a stable final state particle
      if(mcp->Process() == "neutronInelastic") continue;
      if(mcp->Process() == "hadElastic") continue;
      // ignore anything that isn't charged
      unsigned short pdg = abs(mcp->PdgCode());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      tjs.MCPartList.push_back(mcp);
    } // ipart
    
    tjs.MCPartList.shrink_to_fit();
    
    // vector of (mother, daughter) pairs of partList indices
    std::vector<std::pair<unsigned short, unsigned short>> moda;
    
    if(sourcePtclTrackID > 0) {
      // enter grandmother-daughter pairs for primary electrons if we are using shower finding code
      if(tjs.ShowerTag[0] > 1) {
        for(unsigned short ii = 0; ii < tjs.MCPartList.size(); ++ii) {
          // check for the end of the primary particles
          if(tjs.MCPartList[ii]->Mother() != 0) break;
          // check for an electron
          if(abs(tjs.MCPartList[ii]->PdgCode()) != 11) continue;
          int primElectronTrackID = tjs.MCPartList[ii]->TrackId();
          for(unsigned short jj = ii + 1; jj < tjs.MCPartList.size(); ++jj) {
            int trackID = tjs.MCPartList[jj]->TrackId();
            const simb::MCParticle* gmom = bt->TrackIDToMotherParticle(trackID);
            if(gmom == 0 || gmom->TrackId() != primElectronTrackID) continue;
            moda.push_back(std::make_pair(ii, jj));
          } // jj
        } // ii
      } // Shower finding mode
      // Now enter mother-daughter pairs for soft interactions.
      // Daughters appear later in the list so reverse iterate
      for(unsigned short ii = 0; ii < tjs.MCPartList.size(); ++ii) {
        unsigned short dpl = tjs.MCPartList.size() - 1 - ii;
        if(tjs.MCPartList[dpl]->Mother() == 0) continue;
        // ignore previous entries
        int trackID = tjs.MCPartList[ii]->TrackId();
        bool skipit = false;
        for(auto& md : moda) if(md.second == trackID) skipit = true;
        if(skipit) continue;
        int motherID = tjs.MCPartList[dpl]->Mother() + sourcePtclTrackID - 1;
        // count the number of daughters
        unsigned short ndtr = 0;
        for(unsigned short jj = 0; jj < tjs.MCPartList.size(); ++jj) {
          // some processes to ignore. Delta-rays
          if(tjs.MCPartList[jj]->Process() == "hIoni") continue;
          if(tjs.MCPartList[jj]->Process() == "eIoni") continue;
          if(tjs.MCPartList[jj]->Mother() == tjs.MCPartList[dpl]->Mother()) ++ndtr;
        } // jj
        // require only one daughter
        if(ndtr != 1) continue;
        // Then find the mother index
        unsigned short momIndex = USHRT_MAX;
        for(unsigned short jj = 0; jj < tjs.MCPartList.size(); ++jj) {
          if(tjs.MCPartList[jj]->TrackId() == motherID) {
            momIndex = jj;
            break;
          }
        } // jj
        // Mother not found for some reason
        if(momIndex == USHRT_MAX) continue;
        // ensure that mother and daughter have the same PDG code
        if(tjs.MCPartList[momIndex]->PdgCode() != tjs.MCPartList[dpl]->PdgCode()) continue;
        moda.push_back(std::make_pair(ii, dpl));
      } // ii
    } // sourcePtclTrackID >= 0
    
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      std::vector<sim::IDE> ides;
      auto& tcHit = tjs.fHits[iht];
      geo::PlaneID planeID = geo::PlaneID(tcHit.WireID.Cryostat, tcHit.WireID.TPC, tcHit.WireID.Plane);
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)tcHit.WireID.Plane, (int)tcHit.WireID.Wire, (int)tcHit.WireID.TPC, (int)tcHit.WireID.Cryostat);
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
                 tcHit.WireID);
      try {
        bt->HitToSimIDEs(rhit, ides);
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
      unsigned int partListIndex;
      for(partListIndex = 0; partListIndex < tjs.MCPartList.size(); ++partListIndex) if(hitTruTrkID == tjs.MCPartList[partListIndex]->TrackId()) break;
      if(partListIndex == tjs.MCPartList.size()) continue;
      
      // Try to re-assign it to a mother. Note that the mother-daughter pairs
      // are in reverse order, so this loop will transfer all-generation daughters
      // to the (grand) mother
      for(auto& md : moda) if(md.second == partListIndex) partListIndex = md.first;
      tcHit.MCPartListIndex = partListIndex;
    } // iht
    
    // Dump all MC particle info 
    if(tjs.MatchTruth[1] > 1) {
      mf::LogVerbatim myprt("TC");
      myprt<<"part   PDG     TrkID     MomID momIndex KE(MeV) ____Process______   ________Start________  ___Direction____\n";
      for(unsigned short ipart = 0; ipart < tjs.MCPartList.size(); ++ipart) {
        auto& mcp = tjs.MCPartList[ipart];
        unsigned short pdg = abs(mcp->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (mcp->E() - mcp->Mass());
        int motherID = mcp->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipart;
        myprt<<std::setw(6)<<mcp->PdgCode();
        myprt<<std::setw(10)<<mcp->TrackId();
        myprt<<std::setw(10)<<motherID;
        const simb::MCParticle* momMCP = bt->TrackIDToMotherParticle(mcp->TrackId());
        myprt<<std::setw(9)<<momMCP->TrackId();
        myprt<<std::setw(9)<<momMCP->PdgCode();
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<mcp->Process();
        myprt<<std::fixed<<std::setprecision(1);
        myprt<<std::setw(8)<<mcp->Vx();
        myprt<<std::setw(8)<<mcp->Vy();
        myprt<<std::setw(8)<<mcp->Vz();
        TVector3 dir;
        dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
        dir.SetMag(1);
        myprt<<std::fixed<<std::setprecision(2);
        myprt<<std::setw(6)<<dir[0];
        myprt<<std::setw(6)<<dir[1];
        myprt<<std::setw(6)<<dir[2];
        myprt<<" hits: ";
        std::vector<unsigned int> mhits;
        for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
          auto& hit = tjs.fHits[iht];
          if(hit.MCPartListIndex == ipart) mhits.push_back(iht);
        }
        if(!mhits.empty()) {
          if(mhits.size() < 5) {
            myprt<<" hits:";
            for(auto iht : mhits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
          } else {
            myprt<<" nhits = "<<mhits.size();
          }
        }
        myprt<<"\n";
      } // ipart
    }

    // Look for hits without an MC match 
    unsigned int nomat = 0;
    for(auto& hit : tjs.fHits) {
      if(hit.Multiplicity > 1) continue;
      if(hit.GoodnessOfFit <= 0) continue;
      if(hit.MCPartListIndex == USHRT_MAX) ++nomat;
    } // hit
    float noMatFrac = (float)nomat / (float)tjs.fHits.size();
    hist.fUnMatchedHitFrac->Fill(noMatFrac);
//    if(nomat > 0.1) std::cout<<"MTH: reco-true unmatched hitfraction  "<<std::fixed<<std::setprecision(3)<<noMatFrac<<"\n";

  } // MatchTrueHits
*/
  //////////////////////////////////////////
  void TruthMatcher::MatchTruth(const HistStuff& hist, bool fStudyMode)
  {
    // The hits have already been matched to the truth in MatchTrueHits. Here we match reconstructed objects
    // to the truth-matched hits to measure performance    
    
    if(tjs.MatchTruth[0] < 0) return;
    if(tjs.MCPartList.empty()) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;
    
    bool neutrinoVxInFiducialVolume = false;
    bool neutrinoVxReconstructable = false;
    bool neutrinoVxReconstructed = false;
    bool neutrinoVxCorrect = false;
    float truPrimVtxX = -100;
    float truPrimVtxY = -100;
    float truPrimVtxZ = -100;
    
    // See if the interaction is inside the TPC
    simb::Origin_t sourceOrigin = simb::kUnknown;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* mcp = (*ipart).second;
      int trackID = mcp->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(sourcePtclTrackID < 0) {
        if(tjs.MatchTruth[0] == 1) {
          // Look for beam neutrino or single particle
          if(theTruth->Origin() == simb::kBeamNeutrino) {
            fSourceParticleEnergy = 1000 * mcp->E(); // in MeV
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
            fNeutrinoEnergy = 1000 * theTruth->GetNeutrino().Nu().E();
            if(tjs.MatchTruth[1] > 2) std::cout<<"Found beam neutrino E = "<<fNeutrinoEnergy<<" sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode()<<"\n";
          }
          if(theTruth->Origin() == simb::kSingleParticle) {
            fSourceParticleEnergy = 1000 * mcp->E(); // in MeV
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
            truPrimVtxX = mcp->Vx();
            truPrimVtxY = mcp->Vy();
            truPrimVtxZ = mcp->Vz();
            if(tjs.MatchTruth[1] > 0) {
              TVector3 dir;
              dir[0] = mcp->Px(); dir[1] = mcp->Py(); dir[2] = mcp->Pz();
              dir.SetMag(1);
              std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<mcp->PdgCode()<<" Vx "<<(int)mcp->Vx()<<" Vy "<<(int)mcp->Vy()<<" Vz "<<(int)mcp->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            }
          }
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            if(tjs.MatchTruth[1] > 2) std::cout<<" True vertex position "<<(int)mcp->Vx()<<" "<<(int)mcp->Vy()<<" "<<(int)mcp->Vz()<<" energy "<<(int)(1000*mcp->E())<<"\n";
            truPrimVtxX = mcp->Vx();
            truPrimVtxY = mcp->Vy();
            truPrimVtxZ = mcp->Vz();
            neutrinoVxInFiducialVolume = (truPrimVtxX > tjs.XLo && truPrimVtxX < tjs.XHi &&
                                truPrimVtxY > tjs.YLo && truPrimVtxY < tjs.YHi &&
                                truPrimVtxZ > tjs.ZLo && truPrimVtxZ < tjs.ZHi);
          } // sourceOrigin != simb::kUnknown
        } else {
          // look for cosmic rays
          if(theTruth->Origin() == simb::kCosmicRay) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kCosmicRay;
          }
        }
      }
      if(sourcePtclTrackID != -1) break;
    } // ipart

    if(neutrinoVxInFiducialVolume) ++TruVxCounts[0];

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
      tj.MCPartListIndex = USHRT_MAX;
      std::vector<unsigned short> mcpIndex, cnt;
      auto tjhits = PutTrajHitsInVector(tj, kUsedHits);
      for(auto iht : tjhits) {
        auto& hit = tjs.fHits[iht];
        if(hit.MCPartListIndex > tjs.MCPartList.size() - 1) continue;
        unsigned short indx = 0;
        for(indx = 0; indx < mcpIndex.size(); ++indx) if(hit.MCPartListIndex == mcpIndex[indx]) break;
        if(indx == mcpIndex.size()) {
          mcpIndex.push_back(hit.MCPartListIndex);
          cnt.push_back(1);
        } else {
          ++cnt[indx];
        }
      } // iht
      unsigned short maxcnt = 0;
      unsigned short tmpIndex = USHRT_MAX;
      for(unsigned short ii = 0; ii < cnt.size(); ++ii) {
        if(cnt[ii] > maxcnt) {
          maxcnt = cnt[ii];
          tmpIndex = mcpIndex[ii];
        }
      } // ii
      if(tmpIndex > tjs.MCPartList.size() - 1) continue;
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
          mtj.MCPartListIndex = USHRT_MAX;
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
    std::vector<unsigned short> primMCPs;
    for(unsigned short ipart = 0; ipart < tjs.MCPartList.size(); ++ipart) {
      if(tjs.MCPartList[ipart]->Mother() != 0) continue;
      ++nTruPrimary;
      if(CanReconstruct(ipart, 3)) primMCPs.push_back(ipart);
    } // ipart
    
    neutrinoVxReconstructable = (neutrinoVxInFiducialVolume && primMCPs.size() > 1);
    
    if(neutrinoVxReconstructable) {
      ++TruVxCounts[1];
      // Find the closest reconstructed vertex to the true vertex
      float closest = 1;
      unsigned short imTheOne = 0;
      for(auto& aVtx3 : tjs.vtx3) {
        float dx = aVtx3.X - truPrimVtxX;
        float dy = aVtx3.Y - truPrimVtxY;
        float dz = aVtx3.Z - truPrimVtxZ;
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
    
    // Print MCParticles that can be reconstructed
    if(tjs.MatchTruth[1] > 0) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Number of primary particles "<<nTruPrimary<<" Vtx reconstructable? "<<neutrinoVxReconstructable<<" Reconstructed? "<<neutrinoVxReconstructed<<" Correct? "<<neutrinoVxCorrect<<"\n";
      myprt<<"part   PDG     TrkID     MomID KE(MeV)   Process         TrajectoryExtentInPlane_nTruHits \n";
      for(unsigned short ipart = 0; ipart < tjs.MCPartList.size(); ++ipart) {
        // Skip if it can't be reconstructed in 2D
        if(!CanReconstruct(ipart, 2)) continue;
        auto& mcp = tjs.MCPartList[ipart];
        unsigned short pdg = abs(mcp->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (mcp->E() - mcp->Mass());
        int motherID = mcp->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipart;
        myprt<<std::setw(6)<<mcp->PdgCode();
        myprt<<std::setw(10)<<mcp->TrackId();
        myprt<<std::setw(10)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<mcp->Process();
        // print the extent of the particle in each TPC plane
        for(const geo::TPCID& tpcid : tjs.geom->IterateTPCIDs()) {
          geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
          for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
            CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
            auto mcpHits = PutMCPHitsInVector(ipart, inCTP);
            if(mcpHits.size() < 3) continue;
            myprt<<" "<<PrintHitShort(tjs.fHits[mcpHits[0]])<<"-"<<PrintHitShort(tjs.fHits[mcpHits[mcpHits.size() - 1]]);
            myprt<<"_"<<mcpHits.size();
          } // plane
        } // tpcid
        myprt<<"\n";
      } // ipart
    } // printing
    
    // Fill counts for PFParticles
    AccumulatePFPSums();
    
    // Check primary particle reconstruction performance
    if(!primMCPs.empty()) {
      float tsum = 0;
      float eptsum = 0;
      for(auto primMCP : primMCPs) {
        auto& mcp = tjs.MCPartList[primMCP];
        float TMeV = 1000 * (mcp->E() - mcp->Mass());
        mf::LogVerbatim("TC")<<"Primary particle "<<primMCP<<" T "<<TMeV;
        tsum += TMeV;
        for(auto& pfp : tjs.pfps) {
          if(pfp.ID == 0) continue;
          if(pfp.MCPartListIndex != primMCP) continue;
          mf::LogVerbatim("TC")<<"Primary particle "<<primMCP<<" T "<<TMeV<<" pfp EffPur v"<<pfp.EffPur;
          eptsum += pfp.EffPur * TMeV;
        } // pfp
      } // primMCP
      if(tsum == 0) {
        mf::LogVerbatim("TC")<<"MatchTruth: tsum = 0 for all primary particles...";
      } else {
        float ep = eptsum / tsum;
        mf::LogVerbatim("TC")<<"Primary particles EP "<<ep<<" num Primaries "<<primMCPs.size();
      }
    } // primMCPs exist

    // Update the trajectory EP sums
    for(unsigned short ipart = 0; ipart < tjs.MCPartList.size(); ++ipart) {
      auto& mcp = tjs.MCPartList[ipart];
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      unsigned short pdgIndex = PDGCodeIndex(tjs, mcp->PdgCode());
      if(pdgIndex > 4) continue;
      for(const geo::TPCID& tpcid : tjs.geom->IterateTPCIDs()) {
        geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
        for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
          CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
          auto mcpHits = PutMCPHitsInVector(ipart, inCTP);
          if(mcpHits.size() < 3) continue;
          TSums[pdgIndex] += TMeV;
          ++EPCnts[pdgIndex];
          // look for a matched Tj
          int matchedTjID = 0;
          for(auto& tj : tjs.allTraj) {
            if(tj.AlgMod[kKilled]) continue;
            if(tj.MCPartListIndex == ipart) {
              matchedTjID = tj.ID;
              break;
            }
          } // tj
          if(matchedTjID == 0) {
            // Enter 0 in the profile histogram
            hist.fEP_T[pdgIndex]->Fill(TMeV, 0);
            if((float)mcpHits.size() > tjs.MatchTruth[3]) {
              ++nBadEP;
              mf::LogVerbatim myprt("TC");
              myprt<<"pdgIndex "<<pdgIndex<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to mcpIndex["<<ipart<<"]";
              myprt<<" nTrue hits "<<mcpHits.size();
              myprt<<" extent "<<PrintHit(tjs.fHits[mcpHits[0]])<<"-"<<PrintHit(tjs.fHits[mcpHits[mcpHits.size() - 1]]);
              myprt<<" events processed "<<tjs.EventsProcessed;
            }
            continue;
          }
          auto& tj = tjs.allTraj[matchedTjID - 1];
          EPTSums[pdgIndex] += TMeV * tj.EffPur;
          hist.fEP_T[pdgIndex]->Fill(TMeV, tj.EffPur);
          // print debug info?
          if(tj.EffPur < tjs.MatchTruth[2] && (float)mcpHits.size() > tjs.MatchTruth[3]) {
            ++nBadEP;
            mf::LogVerbatim myprt("TC");
            myprt<<"pdgIndex "<<pdgIndex<<" BadEP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
            myprt<<" mcpIndex "<<ipart;
            myprt<<" TMeV "<<(int)TMeV<<" MCP hits "<<mcpHits.size();
            myprt<<" extent "<<PrintHit(tjs.fHits[mcpHits[0]])<<"-"<<PrintHit(tjs.fHits[mcpHits[mcpHits.size() - 1]]);
            myprt<<" events processed "<<tjs.EventsProcessed;
          }
        } // plane
      } // tpcid
    } // ipart
    
    // match 2D vertices (crudely)
    for(auto& tj : tjs.allTraj) {
      // obsolete vertex
      if(tj.AlgMod[kKilled]) continue;
      // require a truth match

      if(tj.MCPartListIndex == USHRT_MAX) continue;
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
    
    for(auto& vx2 : tjs.vtx) if(vx2.ID > 0) ++RecoVx2Count;
    
    if (fStudyMode) {
      // nomatch
      for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
        auto& pfp = tjs.pfps[ipfp];
        if(pfp.ID == 0) continue;
        if(pfp.MCPartListIndex == USHRT_MAX) {
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
      if(pfp.MCPartListIndex == USHRT_MAX) continue;
      short truIndex = PDGCodeIndex(tjs, tjs.MCPartList[pfp.MCPartListIndex]->PdgCode());
      if(truIndex == USHRT_MAX) continue;
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
  void TruthMatcher::AccumulatePFPSums()
  {
    // accumulate performance statistics for PFParticles
    
    // Calculate PFParticle efficiency and purity
    
    std::vector<unsigned short> pfpid(tjs.MCPartList.size(), 0);
    std::vector<unsigned short> nthts(tjs.MCPartList.size(), 0);
    
    // initialize the PFP MCP indices
    for(auto& pfp : tjs.pfps) pfp.MCPartListIndex = USHRT_MAX;
    
    for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
      auto& pfp = tjs.pfps[ipfp];
      if(pfp.ID == 0) continue;
      ++PFP_Cnt;
      pfp.EffPur = 0;
      // put all of the hits into a vector
      std::vector<int> pfpHits;
      for(auto& tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        auto tmp = PutTrajHitsInVector(tj, kUsedHits);
        pfpHits.insert(pfpHits.end(), tmp.begin(), tmp.end());
      } // tj
      if(pfpHits.empty()) continue;
      // find the MC particle with the most hits
      std::vector<unsigned short> mcpIndex, cnt;
      for(auto iht : pfpHits) {
        auto& hit = tjs.fHits[iht];
        if(hit.MCPartListIndex > tjs.MCPartList.size() - 1) continue;
        unsigned short indx = 0;
        for(indx = 0; indx < mcpIndex.size(); ++indx) if(hit.MCPartListIndex == mcpIndex[indx]) break;
        if(indx == mcpIndex.size()) {
          mcpIndex.push_back(hit.MCPartListIndex);
          cnt.push_back(1);
        } else {
          ++cnt[indx];
        }
      } // iht
      unsigned short maxcnt = 0;
      unsigned short tmpIndex = USHRT_MAX;
      for(unsigned short ii = 0; ii < cnt.size(); ++ii) {
        if(cnt[ii] > maxcnt) {
          maxcnt = cnt[ii];
          tmpIndex = mcpIndex[ii];
        }
      } // ii
      if(tmpIndex == USHRT_MAX) continue;
      // count the number of truth-matched hits to this MCParticle in all planes in this tpcid
      float ntru = 0;
      for(auto& hit : tjs.fHits) {
        if(hit.MCPartListIndex != tmpIndex) continue; 
        if(hit.WireID.Cryostat != pfp.TPCID.Cryostat) continue;
        if(hit.WireID.TPC != pfp.TPCID.TPC) continue;
        ++ntru;
      } // hit
      if(ntru == 0) continue;
      float eff = (float)maxcnt / ntru;
      float pur = (float)maxcnt / (float)pfpHits.size();
      float ep = eff * pur;
//      std::cout<<"pfp "<<pfp.ID<<" tmpIndex "<<tmpIndex<<" ntru "<<ntru<<" pfpHits "<<pfpHits.size()<<" maxcnt "<<maxcnt<<" ep "<<ep<<"\n";
      // see if a match already exists and is better
      if(pfpid[tmpIndex] > 0) {
        if(maxcnt < nthts[tmpIndex]) continue;
        // remove the old one
        auto& mpfp = tjs.pfps[pfpid[tmpIndex] - 1];
        mpfp.EffPur = 0;
        mpfp.MCPartListIndex = USHRT_MAX;
        // use the new one
        pfp.EffPur = ep;
        pfp.MCPartListIndex = tmpIndex;
        pfpid[tmpIndex] = pfp.ID;
        nthts[tmpIndex] = maxcnt;
      } else {
        // no previous match
        pfp.EffPur = ep;
        pfp.MCPartListIndex = tmpIndex;
        pfpid[tmpIndex] = pfp.ID;
        nthts[tmpIndex] = maxcnt;
      }
    } // pfp
    
    // accumulate the sums
    for(unsigned short ipart = 0; ipart < tjs.MCPartList.size(); ++ipart) {
      auto& mcp = tjs.MCPartList[ipart];
      float TMeV = 1000 * (mcp->E() - mcp->Mass());
      // See which TPCs this PFParticle should be reconstructed in
      for(const geo::TPCID& tpcid : tjs.geom->IterateTPCIDs()) {
        std::vector<unsigned short> cntInPln(tjs.NumPlanes);
        for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
          auto& hit = tjs.fHits[iht];
          if(hit.MCPartListIndex != ipart) continue;
          if(hit.WireID.Cryostat != tpcid.Cryostat) continue;
          if(hit.WireID.TPC != tpcid.TPC) continue;
          unsigned short plane = hit.WireID.Plane;
          ++cntInPln[plane];
        } // iht
        unsigned short nOKInPln = 0;
        // Require at least 3 hits in at least 2 planes 
        for(unsigned short plane = 0; plane < cntInPln.size(); ++plane) if(cntInPln[plane] > 2) ++nOKInPln;
        if(nOKInPln < 2) continue;
        // This MC Particle should be reconstructed. See if it is
        ++MCP_Cnt;
        MCP_TSum += TMeV;
        unsigned short pfpID = 0;
        for(auto& pfp : tjs.pfps) {
          if(pfp.ID == 0) continue;
          if(pfp.MCPartListIndex == ipart) {
            pfpID = pfp.ID;
            break;
          }
        } // pfp
        if(pfpID == 0) {
          if(TMeV > 30) {
            mf::LogVerbatim myprt("TC");
            myprt<<"BadPFP PDGCode "<<mcp->PdgCode()<<" TMeV "<<(int)TMeV;
            myprt<<" Tjs matched to MCParticle: ";
            for(auto& tj : tjs.allTraj) {
              if(tj.AlgMod[kKilled]) continue;
              if(tj.MCPartListIndex == ipart) myprt<<" "<<tj.ID<<" EP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
            } // tj
            myprt<<" events processed "<<tjs.EventsProcessed;
          } // TMeV > 30
          continue;
        }
        auto& pfp = tjs.pfps[pfpID - 1];
        MCP_EPTSum += pfp.EffPur * TMeV;
        ++PFP_CntMat;
      } // tpcid
    } // ipart
      

  } // AccumulatePFPSums
  
  ////////////////////////////////////////////////
  void TruthMatcher::PrintResults(int eventNum) const
  {
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
    // Vertex reconstruction
    if(MCP_TSum > 0 && PFP_Cnt > 0) {
      // PFParticle statistics
      float ep = MCP_EPTSum / MCP_TSum;
      float nofrac = 1 - (PFP_CntMat / PFP_Cnt);
      myprt<<" PFP "<<ep<<" MCP cnt "<<(int)MCP_Cnt<<" PFP cnt "<<(int)PFP_Cnt<<" noMatFrac "<<nofrac;
    }
    myprt<<" VxCount";
    for(auto cnt : TruVxCounts) myprt<<" "<<cnt;
    float vx2Cnt = 0;
    if(tjs.EventsProcessed > 0) vx2Cnt = (float)RecoVx2Count / (float)tjs.EventsProcessed;
    myprt<<" RecoVx2Cnt/Evt "<<std::fixed<<std::setprecision(1)<<vx2Cnt;
  } // PrintResults
  
  ////////////////////////////////////////////////
  bool TruthMatcher::CanReconstruct(unsigned short mcpIndex, unsigned short nDimensions)
  {
    // returns true if the MCParticle can be reconstructed in nDimensions
    if(mcpIndex > tjs.MCPartList.size() - 1) return false;
    if(nDimensions < 2 || nDimensions > 3) return false;
    
    std::vector<unsigned short> cntInPln(tjs.NumPlanes);
    for(auto& hit : tjs.fHits) if(hit.MCPartListIndex == mcpIndex) ++cntInPln[hit.WireID.Plane];
    unsigned short nPlnOK = 0;
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) if(cntInPln[plane] > 2) ++nPlnOK;
    if(nPlnOK < nDimensions - 1) return false;
    return true;
  } // CanReconstruct
  
  ////////////////////////////////////////////////
  std::vector<unsigned int> TruthMatcher::PutMCPHitsInVector(unsigned short mcpIndex, CTP_t inCTP)
  {
    // put the hits matched to the MCParticle into a vector in the requested CTP
    std::vector<unsigned int> hitVec;
    if(mcpIndex > tjs.MCPartList.size() - 1) return hitVec;
    geo::PlaneID planeID = DecodeCTP(inCTP);
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      auto& hit = tjs.fHits[iht];
      if(hit.MCPartListIndex != mcpIndex) continue;
      if(hit.WireID.Plane != planeID.Plane) continue;
      if(hit.WireID.TPC != planeID.TPC) continue;
      if(hit.WireID.Cryostat != planeID.Cryostat) continue;
      hitVec.push_back(iht);
    } // iht
    return hitVec;
  } // PutMCPHitsInVector

  ////////////////////////////////////////////////
  void MCParticleListUtils::MakeTruTrajPoint(unsigned short MCParticleListIndex, TrajPoint& tp)
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
  unsigned short MCParticleListUtils::MCParticleStartTjID(unsigned short MCParticleListIndex, CTP_t inCTP)
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
  unsigned short MCParticleListUtils::GetMCPartListIndex(const ShowerStruct& ss, unsigned short& nTruHits)
  {
    // Returns the index of the MCParticle that has the most number of matches
    // to the hits in this shower
    
    if(tjs.MCPartList.empty()) return USHRT_MAX;
    if(ss.TjIDs.empty()) return USHRT_MAX;
    
    std::vector<unsigned short> pListCnt(tjs.MCPartList.size());
    
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
    
    unsigned short pIndex = USHRT_MAX;
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
  unsigned short MCParticleListUtils::GetMCPartListIndex(const Trajectory& tj, unsigned short& nTruHits)
  {
    // Returns the index of the MCParticle that has the most number of matches
    // to the hits in this trajectory
    
    if(tjs.MCPartList.empty()) return USHRT_MAX;
    
    // Check all hits associated with this Tj
    std::vector<unsigned short> pListCnt(tjs.MCPartList.size());
    
    for(auto& tp : tj.Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        // ignore unmatched hits
        if(tjs.fHits[iht].MCPartListIndex > tjs.MCPartList.size() - 1) continue;
        ++pListCnt[tjs.fHits[iht].MCPartListIndex];
      } // ii
    } // pt
    
    unsigned short pIndex = USHRT_MAX;
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
