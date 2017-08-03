#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"


#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "larsim/MCCheater/BackTracker.h"


namespace tca {

  //////////////////////////////////////////
  void TruthMatcher::Initialize()
  {
     // Initialize the variables used to calculate Efficiency * Purity (aka EP) for matching to truth
    for(unsigned short pdgIndex = 0; pdgIndex < 6; ++pdgIndex) {
      EPTSums[pdgIndex] = 0;
      EPSums[pdgIndex] = 0;
      EPCnts[pdgIndex] = 0;
    }
    nTruPrimaryVtxOK = 0;
    nTruPrimaryVtxReco = 0;
  } // Initialize

  //////////////////////////////////////////
  void TruthMatcher::MatchTrueHits()
  {
    // Matches reco hits to MC true tracks and puts the association into
    // TCHit TruTrkID. This code is almost identical to the first part of MatchTruth.
    
    if(tjs.MatchTruth[0] < 0) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    tjs.MCPartList.clear();
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;
    
    simb::Origin_t sourceOrigin = simb::kUnknown;
    // partList is the vector of MC particles that we want to use
    tjs.MCPartList.reserve(plist.size());
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* part = (*ipart).second;
      int trackID = part->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(sourcePtclTrackID < 0) {
        if(tjs.MatchTruth[0] == 1) {
          // Look for beam neutrino or single particle
          if(theTruth->Origin() == simb::kBeamNeutrino) {
            fSourceParticleEnergy = 1000 * part->E();
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
            if(tjs.MatchTruth[1] > 2) std::cout<<"Found beam neutrino sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<"\n";
          } // beam neutrino
          if(theTruth->Origin() == simb::kSingleParticle) {
            fSourceParticleEnergy = 1000 * part->E();
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
            if(tjs.MatchTruth[1] > 0) {
              TVector3 dir;
              dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
              dir.SetMag(1);
              std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<" Vx "<<(int)part->Vx()<<" Vy "<<(int)part->Vy()<<" Vz "<<(int)part->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            }
          } // single particle
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            if(tjs.MatchTruth[1] > 2) std::cout<<"True vertex position "<<(int)part->Vx()<<" "<<(int)part->Vy()<<" "<<(int)part->Vz()<<"\n";
          } // sourceOrigin != simb::kUnknown
        } else {
          // look for cosmic rays
          if(theTruth->Origin() == simb::kCosmicRay) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kCosmicRay;
          }
        }
      }
      // ignore anything that has the incorrect origin
      if(theTruth->Origin() != sourceOrigin) continue;
      // ignore processes that aren't a stable final state particle
      if(part->Process() == "neutronInelastic") continue;
      if(part->Process() == "hadElastic") continue;
      // ignore anything that isn't charged
      unsigned short pdg = abs(part->PdgCode());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      tjs.MCPartList.push_back(part);
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
      // Now enter mother-daughter pairs for soft interactions
      // daughters appear later in the list so reverse iterate
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
          // some processes to ignore
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
      TCHit& hit = tjs.fHits[iht];
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)hit.WireID.Plane, (int)hit.WireID.Wire, (int)hit.WireID.TPC, (int)hit.WireID.Cryostat);
      double startTick = hit.PeakTime - hit.RMS;
      double endTick = hit.PeakTime + hit.RMS;
      unsigned short hitTruTrkID = 0;
      // get a list of track IDEs that are close to this hit
      std::vector<sim::TrackIDE> tides;
      bt->ChannelToTrackIDEs(tides, channel, startTick, endTick);
      // Declare a match to the one which has an energy fraction > 0.5
      for(auto itide = tides.begin(); itide != tides.end(); ++itide) {
        if(itide->energyFrac > 0.5) {
          hitTruTrkID = itide->trackID;
          break;
        }
      } // itid
      // not matched (confidently) to a MC track
      if(hitTruTrkID == 0) continue;
      // find out which partList entry corresponds to this track ID
      unsigned short partListIndex;
      for(partListIndex = 0; partListIndex < tjs.MCPartList.size(); ++partListIndex) if(hitTruTrkID == tjs.MCPartList[partListIndex]->TrackId()) break;
      if(partListIndex == tjs.MCPartList.size()) {
        //        std::cout<<"MatchTrueHits: Didn't find partList entry for MC Track ID "<<hitTruTrkID<<"\n";
        continue;
      }
      
      // Try to re-assign it to a mother. Note that the mother-daughter pairs
      // are in reverse order, so this loop will transfer all-generation daughters
      // to the (grand) mother
      for(auto& md : moda) if(md.second == partListIndex) partListIndex = md.first;
      hit.MCPartListIndex = partListIndex;
    } // iht
    
    
    if(tjs.MatchTruth[1] > 1) {
      mf::LogVerbatim myprt("TC");
      myprt<<"part   PDG TrkID MomID KE(MeV)   Process         Trajectory_extent_in_plane \n";
      for(unsigned short ipl = 0; ipl < tjs.MCPartList.size(); ++ipl) {
        unsigned short pdg = abs(tjs.MCPartList[ipl]->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (tjs.MCPartList[ipl]->E() - tjs.MCPartList[ipl]->Mass());
        int motherID = tjs.MCPartList[ipl]->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipl;
        myprt<<std::setw(6)<<tjs.MCPartList[ipl]->PdgCode();
        myprt<<std::setw(6)<<tjs.MCPartList[ipl]->TrackId();
        myprt<<std::setw(6)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<tjs.MCPartList[ipl]->Process();
        // print the extent of the particle in each plane
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          unsigned int fht = UINT_MAX;
          unsigned int lht = 0;
          for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
            if(tjs.fHits[iht].WireID.Plane != plane) continue;
            unsigned short partListIndex = ipl;
            // Look for the real mother
            for(auto& md : moda) if(md.second == partListIndex) partListIndex = md.first;
            if(tjs.fHits[iht].MCPartListIndex != partListIndex) continue;
            if(fht == UINT_MAX) fht = iht;
            lht = iht;
          } // iht
          if(fht == UINT_MAX) continue;
          myprt<<" "<<PrintHitShort(tjs.fHits[fht])<<"-"<<PrintHitShort(tjs.fHits[lht]);
        } // plane
        myprt<<"\n";
      } // ipl
    }
    
  } // MatchTrueHits
  
  //////////////////////////////////////////
  void TruthMatcher::MatchTruth(const HistStuff& hist, unsigned int& fEventsProcessed)
  {
    
    if(tjs.MatchTruth[0] < 0) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    // set true if there is a reconstructed 3D vertex within 1 cm of the true vertex
    bool nuVtxRecoOK = false;
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;
    
    simb::Origin_t sourceOrigin = simb::kUnknown;
    std::vector<simb::MCParticle*> partList;
    // partList is the vector of MC particles that we want to use
    partList.reserve(plist.size());
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* part = (*ipart).second;
      int trackID = part->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(sourcePtclTrackID < 0) {
        if(tjs.MatchTruth[0] == 1) {
          // Look for beam neutrino or single particle
          if(theTruth->Origin() == simb::kBeamNeutrino) {
            fSourceParticleEnergy = 1000 * part->E(); // in MeV
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
            fNeutrinoEnergy = 1000 * theTruth->GetNeutrino().Nu().E();
            if(tjs.MatchTruth[1] > 2) std::cout<<"Found beam neutrino E = "<<fNeutrinoEnergy<<" sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<"\n";
          }
          if(theTruth->Origin() == simb::kSingleParticle) {
            fSourceParticleEnergy = 1000 * part->E(); // in MeV
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
            if(tjs.MatchTruth[1] > 0) {
              TVector3 dir;
              dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
              dir.SetMag(1);
              std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<" Vx "<<(int)part->Vx()<<" Vy "<<(int)part->Vy()<<" Vz "<<(int)part->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            }
          }
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            if(tjs.MatchTruth[1] > 2) std::cout<<" True vertex position "<<(int)part->Vx()<<" "<<(int)part->Vy()<<" "<<(int)part->Vz()<<" energy "<<(int)(1000*part->E())<<"\n";
            for(auto& aVtx3 : tjs.vtx3) {
              hist.fNuVtx_dx->Fill(part->Vx() - aVtx3.X);
              hist.fNuVtx_dy->Fill(part->Vy() - aVtx3.Y);
              hist.fNuVtx_dz->Fill(part->Vz() - aVtx3.Z);
              if(std::abs(part->Vx()-aVtx3.X) < 1 && std::abs(part->Vy()-aVtx3.Y) < 1 && std::abs(part->Vz()-aVtx3.Z) < 1) {
                nuVtxRecoOK = true;
                float score = 0;
                for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
                  if(aVtx3.Vx2ID[ipl] == 0) continue;
                  unsigned short iv2 = aVtx3.Vx2ID[ipl] - 1;
                  score += tjs.vtx[iv2].Score;
                } // ipl
                hist.fNuVtx_Score->Fill(score);
                hist.fNuVtx_Enu_Score_p->Fill(fNeutrinoEnergy, score);
              }
            } // aVtx3
          } // sourceOrigin != simb::kUnknown
        } else {
          // look for cosmic rays
          if(theTruth->Origin() == simb::kCosmicRay) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kCosmicRay;
          }
        }
      }
      // ignore anything that has the incorrect origin
      if(theTruth->Origin() != sourceOrigin) continue;
      // ignore processes that aren't a stable final state particle
      if(part->Process() == "neutronInelastic") continue;
      if(part->Process() == "hadElastic") continue;
      // ignore anything that isn't charged
      unsigned short pdg = abs(part->PdgCode());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      partList.push_back(part);
    } // ipart
    
    if(tjs.MatchTruth[1] > 2) {
      for(unsigned int ii = 0; ii < partList.size(); ++ii) {
        int trackID = partList[ii]->TrackId();
        const simb::MCParticle* gmom = bt->TrackIDToMotherParticle(trackID);
        std::cout<<ii<<" PDG Code  "<<partList[ii]->PdgCode()<<" TrackId "<<trackID<<" Mother  "<<partList[ii]->Mother()<<" Grandmother "<<gmom->TrackId()<<" Process "<<partList[ii]->Process()<<"\n";
      } // ii
    }
    
    // vector of (mother, daughter) pairs of TrackIds
    std::vector<std::pair<int, int>> moda;
    
    if(sourcePtclTrackID > 0) {
      // enter grandmother-daughter pairs for primary electrons if we are using shower finding code
      if(tjs.ShowerTag[0] > 1) {
        for(unsigned short ii = 0; ii < partList.size(); ++ii) {
          // check for the end of the primary particles
          if(partList[ii]->Mother() != 0) break;
          // check for an electron
          if(abs(partList[ii]->PdgCode()) != 11) continue;
          int primElectronTrackID = partList[ii]->TrackId();
          for(unsigned short jj = ii + 1; jj < partList.size(); ++jj) {
            int trackID = partList[jj]->TrackId();
            const simb::MCParticle* gmom = bt->TrackIDToMotherParticle(trackID);
            if(gmom == 0 || gmom->TrackId() != primElectronTrackID) continue;
            moda.push_back(std::make_pair(primElectronTrackID, trackID));
          } // jj
        } // ii
      } // Shower finding mode
      // Now enter mother-daughter pairs for soft interactions
      // daughters appear later in the list so reverse iterate
      for(unsigned short ii = 0; ii < partList.size(); ++ii) {
        unsigned short dpl = partList.size() - 1 - ii;
        if(partList[dpl]->Mother() == 0) continue;
        // ignore previous entries
        int trackID = partList[ii]->TrackId();
        bool skipit = false;
        for(auto& md : moda) if(md.second == trackID) skipit = true;
        if(skipit) continue;
        int motherID = partList[dpl]->Mother() + sourcePtclTrackID - 1;
        // count the number of daughters
        unsigned short ndtr = 0;
        for(unsigned short jj = 0; jj < partList.size(); ++jj) {
          // some processes to ignore
          if(partList[jj]->Process() == "hIoni") continue;
          if(partList[jj]->Process() == "eIoni") continue;
          if(partList[jj]->Mother() == partList[dpl]->Mother()) ++ndtr;
        } // jj
        // require only one daughter
        if(ndtr != 1) continue;
        // Then find the mother index
        unsigned short momIndex = USHRT_MAX;
        for(unsigned short jj = 0; jj < partList.size(); ++jj) {
          if(partList[jj]->TrackId() == motherID) {
            momIndex = jj;
            break;
          }
        } // jj
        // Mother not found for some reason
        if(momIndex == USHRT_MAX) continue;
        // ensure that mother and daughter have the same PDG code
        if(partList[momIndex]->PdgCode() != partList[dpl]->PdgCode()) continue;
        moda.push_back(std::make_pair(partList[momIndex]->TrackId(), partList[dpl]->TrackId()));
      } // ii
    } // sourcePtclTrackID >= 0
    
    if(tjs.MatchTruth[1] > 2 && !moda.empty()) {
      std::cout<<"Mother-Daughter track IDs\n";
      unsigned short cnt = 0;
      for(auto& md : moda) {
        std::cout<<" "<<md.first<<"-"<<md.second;
        ++cnt;
        if(!(cnt % 20)) std::cout<<"\n";
      } // md
      std::cout<<"\n";
    } // tjs.MatchTruth[1] > 2
    
    // Match all hits to the truth. Put the MC track ID in a temp vector
    std::vector<int> hitTruTrkID(tjs.fHits.size());
    // Prepare to count of the number of hits matched to each MC Track in each plane
    std::vector<std::vector<unsigned short>> nMatchedHitsInPartList(plist.size());
    for(unsigned short ipl = 0; ipl < plist.size(); ++ipl) nMatchedHitsInPartList[ipl].resize(tjs.NumPlanes);
    // and make a list of the TJs and hit count for each MC Track
    std::vector<std::vector<std::array<unsigned short, 2>>> nMatchedHitsInTj(partList.size());
    
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      TCHit& hit = tjs.fHits[iht];
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)hit.WireID.Plane, (int)hit.WireID.Wire, (int)hit.WireID.TPC, (int)hit.WireID.Cryostat);
      double startTick = hit.PeakTime - hit.RMS;
      double endTick = hit.PeakTime + hit.RMS;
      unsigned short plane = tjs.fHits[iht].WireID.Plane;
      // get a list of track IDEs that are close to this hit
      std::vector<sim::TrackIDE> tides;
      bt->ChannelToTrackIDEs(tides, channel, startTick, endTick);
      // Declare a match to the one which has an energy fraction > 0.5
      for(auto itide = tides.begin(); itide != tides.end(); ++itide) {
        if(itide->energyFrac > 0.5) {
          hitTruTrkID[iht] = itide->trackID;
          break;
        }
      } // itid
      // not matched (confidently) to a MC track
      if(hitTruTrkID[iht] == 0) continue;
      
      // Try to re-assign it to a mother. Note that the mother-daughter pairs
      // are in reverse order, so this loop will transfer all-generation daughters
      // to the (grand) mother
      for(auto& md : moda) if(md.second == hitTruTrkID[iht]) hitTruTrkID[iht] = md.first;
      
      // count the number of matched hits for each MC track in each plane
      for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
        if(hitTruTrkID[iht] == partList[ipl]->TrackId()) {
          ++nMatchedHitsInPartList[ipl][plane];
          if(tjs.fHits[iht].InTraj > 0) {
            unsigned short itj = tjs.fHits[iht].InTraj - 1;
            bool gotit = false;
            for(auto& hitInTj : nMatchedHitsInTj[ipl]) {
              if(hitInTj[0] == itj) {
                ++hitInTj[1];
                gotit = true;
              }
            } //  hitInTj
            if(!gotit) {
              std::array<unsigned short, 2> tmp {itj, 1};
              nMatchedHitsInTj[ipl].push_back(tmp);
            }
          } // inTraj > 0
        } // hit matched to partList
      } // ipl
    } // iht
    
    // remove partList elements that have no matched hits
    std::vector<simb::MCParticle*> newPartList;
    std::vector<std::vector<unsigned short>> newnMatchedHitsInPartList;
    std::vector<std::vector<std::array<unsigned short, 2>>> newnMatchedHitsInTj;
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      unsigned short nht = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) nht += nMatchedHitsInPartList[ipl][plane];
      if(nht == 0) continue;
      newPartList.push_back(partList[ipl]);
      newnMatchedHitsInPartList.push_back(nMatchedHitsInPartList[ipl]);
      newnMatchedHitsInTj.push_back(nMatchedHitsInTj[ipl]);
    } // ipl
    partList = newPartList;
    nMatchedHitsInPartList = newnMatchedHitsInPartList;
    nMatchedHitsInTj = newnMatchedHitsInTj;
    
    // count the number of primary tracks that have at least 3 hits in at least 2 planes
    nTruPrimary = 0;
    nTruPrimaryOK = 0;
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      if(partList[ipl]->Mother() != 0) continue;
      ++nTruPrimary;
      unsigned short nInPln = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(nMatchedHitsInPartList[ipl][plane] > 2) ++nInPln;
      } // plane
      if(nInPln > 1) ++nTruPrimaryOK;
    } // ipl
    
    if(nTruPrimaryOK > 1) {
      // More than one reconstructable primaries so there must a reconstructable neutrino vertex
      ++nTruPrimaryVtxOK;
      // was it reconstructed?
      if(nuVtxRecoOK) ++nTruPrimaryVtxReco;
      if(fSourceParticleEnergy > 0 && !nuVtxRecoOK) mf::LogVerbatim("TC")<<"BadVtx fSourceParticleEnergy "<<std::fixed<<std::setprecision(2)<<fSourceParticleEnergy<<" events processed "<<fEventsProcessed;
    }
    
    if(tjs.MatchTruth[1] > 1) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Number of primary particles "<<nTruPrimary<<" Number reconstructable "<<nTruPrimaryOK<<" Found neutrino vertex? "<<nuVtxRecoOK<<"\n";
      myprt<<"part   PDG TrkID MomID KE(MeV)   Process         Trajectory_extent_in_plane \n";
      for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
        unsigned short pdg = abs(partList[ipl]->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (partList[ipl]->E() - partList[ipl]->Mass());
        int motherID = partList[ipl]->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipl;
        myprt<<std::setw(6)<<partList[ipl]->PdgCode();
        myprt<<std::setw(6)<<partList[ipl]->TrackId();
        myprt<<std::setw(6)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<partList[ipl]->Process();
        // print the extent of the particle in each plane
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          unsigned int fht = UINT_MAX;
          unsigned int lht = 0;
          for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
            if(tjs.fHits[iht].WireID.Plane != plane) continue;
            unsigned short momTrackID = partList[ipl]->TrackId();
            // Look for the real mother
            for(auto& md : moda) if(md.second == momTrackID) momTrackID = md.first;
            if(hitTruTrkID[iht] != momTrackID) continue;
            if(fht == UINT_MAX) fht = iht;
            lht = iht;
          } // iht
          if(fht == UINT_MAX) continue;
          myprt<<" "<<PrintHitShort(tjs.fHits[fht])<<"-"<<PrintHitShort(tjs.fHits[lht]);
        } // plane
        myprt<<"\n";
      } // ipl
    }
    
    // Declare a TJ - partlist match for the trajectory which has the most true hits.
    // another temp vector for the one-to-one match
    std::vector<std::vector<unsigned short>> partListToTjID(partList.size());
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) partListToTjID[ipl].resize(tjs.NumPlanes);
    
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(nMatchedHitsInPartList[ipl][plane] < 2) continue;
        unsigned short mostHits = 0;
        unsigned short tjWithMostHits = USHRT_MAX;
        for(unsigned short ii = 0; ii < nMatchedHitsInTj[ipl].size(); ++ii) {
          unsigned short itj = nMatchedHitsInTj[ipl][ii][0];
          geo::PlaneID planeID = DecodeCTP(tjs.allTraj[itj].CTP);
          // ensure we only check Tjs in the correct plane
          if(planeID.Plane != plane) continue;
          unsigned short nMatHits = nMatchedHitsInTj[ipl][ii][1];
          if(nMatHits > mostHits) {
            mostHits = nMatHits;
            tjWithMostHits = itj;
          }
        } // ii
        if(tjWithMostHits == USHRT_MAX) continue;
        // Count the total number of used hits in the TJ
        auto tmp = PutTrajHitsInVector(tjs.allTraj[tjWithMostHits], kUsedHits);
        if(tjs.allTraj[tjWithMostHits].ParentTrajID > 0) {
          // This is a daughter trajectory that has more truth matched hits than the parent
          unsigned short ptj = tjs.allTraj[tjWithMostHits].ParentTrajID - 1;
          // add the parent hits to the vector of daughter hits
          auto ptmp = PutTrajHitsInVector(tjs.allTraj[ptj], kUsedHits);
          tmp.insert(tmp.end(), ptmp.begin(), ptmp.end());
          // revise the mostHits count. To do this we need to find the parent tj index in nMatchedHitsInTj
          for(unsigned short ii = 0; ii < nMatchedHitsInTj[ipl].size(); ++ii) {
            unsigned short itj = nMatchedHitsInTj[ipl][ii][0];
            if(itj == ptj) {
              mostHits += nMatchedHitsInTj[ipl][ii][1];
              // re-direct the calculation to the parent
              tjWithMostHits = ptj;
              break;
            } // found the parent tj
          } // ii
        } // deal with daughters
        // count the number matched to a true particle
        float nTjHits = 0;
        for(auto& iht : tmp) if(hitTruTrkID[iht] > 0) ++nTjHits;
        float nTruHits = nMatchedHitsInPartList[ipl][plane];
        float nTjTruRecHits = mostHits;
        float eff = nTjTruRecHits / nTruHits;
        float pur = nTjTruRecHits / nTjHits;
        float effpur = eff * pur;
        // This overwrites any previous match that has poorer efficiency * purity
        if(effpur > tjs.allTraj[tjWithMostHits].EffPur) {
          tjs.allTraj[tjWithMostHits].MCPartListIndex = ipl;
          tjs.allTraj[tjWithMostHits].EffPur = effpur;
          partListToTjID[ipl][plane] = tjs.allTraj[tjWithMostHits].ID;
        }
      } // plane
    } // ipl
    
    // Update the EP sums
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      float TMeV = 1000 * (partList[ipl]->E() - partList[ipl]->Mass());
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        // require at least 2 matched hits
        if(nMatchedHitsInPartList[ipl][plane] < 2) continue;
        unsigned short pdgIndex = PDGCodeIndex(tjs, partList[ipl]->PdgCode());
        // count the number of EP sums for this PDG code
        EPSums[pdgIndex] += TMeV;
        ++EPCnts[pdgIndex];
        // find the first and last matched hit in this plane
        unsigned int fht = UINT_MAX;
        unsigned int lht = 0;
        // find the first and last matched hit in this plane
        for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
          if(tjs.fHits[iht].WireID.Plane != plane) continue;
          unsigned short momTrackID = partList[ipl]->TrackId();
          // Look for the real mother
          for(auto& md : moda) if(md.second == momTrackID) momTrackID = md.first;
          if(hitTruTrkID[iht] != momTrackID) continue;
          if(fht == UINT_MAX) fht = iht;
          lht = iht;
        } // iht
        if(fht == UINT_MAX) continue;
        if(partListToTjID[ipl][plane] == 0) {
          // Enter 0 in the profile histogram
          hist.fEP_T[pdgIndex]->Fill(TMeV, 0);
          if(nMatchedHitsInPartList[ipl][plane] > tjs.MatchTruth[3]) {
            mf::LogVerbatim myprt("TC");
            myprt<<"pdgIndex "<<pdgIndex<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to partList["<<ipl<<"]";
            myprt<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipl][plane];
            myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht])<<" events processed "<<fEventsProcessed;
          }
          continue;
        }
        unsigned short itj = partListToTjID[ipl][plane] - 1;
        EPTSums[pdgIndex] += TMeV * tjs.allTraj[itj].EffPur;
        hist.fEP_T[pdgIndex]->Fill(TMeV, tjs.allTraj[itj].EffPur);
        // print out some debugging information if the EP was pitiful and the number of matched hits is large
        if(tjs.allTraj[itj].EffPur < tjs.MatchTruth[2] && nMatchedHitsInPartList[ipl][plane] > tjs.MatchTruth[3]) {
          mf::LogVerbatim myprt("TC");
          myprt<<"pdgIndex "<<pdgIndex<<" BadEP "<<std::fixed<<std::setprecision(2)<<tjs.allTraj[itj].EffPur;
          myprt<<" TMeV "<<(int)TMeV<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipl][plane];
          myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht])<<" events processed "<<fEventsProcessed;
          // print alg names
          for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tjs.allTraj[itj].AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        }
        // check for a bad match to a primary electron shower
        if(tjs.ShowerTag[0] > 1 && pdgIndex == 0 && partList[ipl]->Mother() == 0 && fSourceParticleEnergy > 50) {
          Trajectory& ptj = tjs.allTraj[itj];
          // determine if this is a parent of a shower Tj
          int dtrID = 0;
          for(auto& tj : tjs.allTraj) {
            if(!tj.AlgMod[kShowerTj]) continue;
            if(tj.ParentTrajID == ptj.ID) {
              dtrID = tj.ID;
              break;
            }
          } // tj
          if(dtrID == 0) mf::LogVerbatim("TC")<<"BadShower Primary electron -> Traj "<<ptj.ID<<"_"<<ptj.CTP<<" Wrong shower parent. Events processed "<<fEventsProcessed;
        } // check primary electron shower
        // histogram the MC-reco stopping point difference
        unsigned short endPt = tjs.allTraj[itj].EndPt[0];
        float recoWire0 = tjs.allTraj[itj].Pts[endPt].Pos[0];
        endPt = tjs.allTraj[itj].EndPt[1];
        float recoWire1 = tjs.allTraj[itj].Pts[endPt].Pos[0];
        float trueFirstWire = tjs.fHits[fht].WireID.Wire;
        float trueLastWire = tjs.fHits[lht].WireID.Wire;
        // decide which ends should be compared
        if(std::abs(recoWire0 - trueFirstWire) < std::abs(recoWire1 - trueFirstWire)) {
          hist.fdWire[pdgIndex]->Fill(recoWire0 - trueFirstWire);
          hist.fdWire[pdgIndex]->Fill(recoWire1 - trueLastWire);
        } else {
          hist.fdWire[pdgIndex]->Fill(recoWire1 - trueFirstWire);
          hist.fdWire[pdgIndex]->Fill(recoWire0 - trueLastWire);
        }
      } // plane
    } // ipl
    
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
    } // vx2
    
  } // MatchTruth
  
  ////////////////////////////////////////////////
  void TruthMatcher::PrintResults(int eventNum) const
  {
    mf::LogVerbatim myprt("TC");
    myprt<<"Event "<<eventNum;
    float sum = 0;
    float sumt = 0;
    for(unsigned short pdgIndex = 0; pdgIndex < EPSums.size(); ++pdgIndex) {
      if(EPSums[pdgIndex] == 0) continue;
      if(pdgIndex == 0) myprt<<" Electron";
      if(pdgIndex == 1) myprt<<" Muon";
      if(pdgIndex == 2) myprt<<" Pion";
      if(pdgIndex == 3) myprt<<" Kaon";
      if(pdgIndex == 4) myprt<<" Proton";
      float ave = EPTSums[pdgIndex] / (float)EPSums[pdgIndex];
      myprt<<" "<<std::fixed<<std::setprecision(2)<<ave;
      myprt<<" "<<EPCnts[pdgIndex];
      if(pdgIndex > 0) {
        sum  += EPSums[pdgIndex];
        sumt += EPTSums[pdgIndex];
      }
    } // pdgIndex
    if(sum > 0) myprt<<" MuPiKP "<<std::fixed<<std::setprecision(2)<<sumt / sum;
  } // PrintResults

  
  ////////////////////////////////////////////////
  void MCParticleListUtils::MakeTruTrajPoint(unsigned short MCParticleListIndex, TrajPoint& tp)
  {
    // Creates a trajectory point at the start of the MCParticle with index MCParticleListIndex. The
    // projected length of the MCParticle in the plane coordinate system is stored in TruTp.Delta.
    // The calling function should specify the CTP in which this TP resides.
    
    if(MCParticleListIndex > tjs.MCPartList.size() - 1) return;
    
    const simb::MCParticle* part = tjs.MCPartList[MCParticleListIndex];
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    
    tp.Pos[0] = tjs.geom->WireCoordinate(part->Vy(), part->Vz(), planeID);
    tp.Pos[1] = tjs.detprop->ConvertXToTicks(part->Vx(), planeID) * tjs.UnitsPerTick;
    
    TVector3 dir;
    dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
    if(dir.Mag() == 0) return;
    dir.SetMag(1);
    TVector3 pos;
    pos[0] = part->Vx() + 100 * dir[0];
    pos[1] = part->Vy() + 100 * dir[1];
    pos[2] = part->Vz() + 100 * dir[2];
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
    
    const simb::MCParticle* part = tjs.MCPartList[MCParticleListIndex];
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    TrajPoint truTp;
    truTp.Pos[0] = tjs.geom->WireCoordinate(part->Vy(), part->Vz(), planeID);
    truTp.Pos[1] = tjs.detprop->ConvertXToTicks(part->Vx(), planeID) * tjs.UnitsPerTick;
    
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
