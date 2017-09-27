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
    
    if(tjs.MatchTruth[0] < 0) return;
    if(tjs.fHits.empty()) return;
    
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
      myprt<<"part   PDG     TrkID     MomID KE(MeV) ____Process______   ________Start________  ___Direction____\n";
      for(unsigned short ipart = 0; ipart < tjs.MCPartList.size(); ++ipart) {
        auto& part = tjs.MCPartList[ipart];
        unsigned short pdg = abs(part->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (part->E() - part->Mass());
        int motherID = part->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipart;
        myprt<<std::setw(6)<<part->PdgCode();
        myprt<<std::setw(10)<<part->TrackId();
        myprt<<std::setw(10)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<part->Process();
        myprt<<std::fixed<<std::setprecision(1);
        myprt<<std::setw(8)<<part->Vx();
        myprt<<std::setw(8)<<part->Vy();
        myprt<<std::setw(8)<<part->Vz();
        TVector3 dir;
        dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
        dir.SetMag(1);
        myprt<<std::fixed<<std::setprecision(2);
        myprt<<std::setw(6)<<dir[0];
        myprt<<std::setw(6)<<dir[1];
        myprt<<std::setw(6)<<dir[2];
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
  
  //////////////////////////////////////////
  void TruthMatcher::MatchTruth(const HistStuff& hist, bool fStudyMode)
  {
    
    if(tjs.MatchTruth[0] < 0) return;
    
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
            truPrimVtxX = part->Vx();
            truPrimVtxY = part->Vy();
            truPrimVtxZ = part->Vz();
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
            truPrimVtxX = part->Vx();
            truPrimVtxY = part->Vy();
            truPrimVtxZ = part->Vz();
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

    if(neutrinoVxInFiducialVolume) ++TruVxCounts[0];

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
    for(unsigned short ipart = 0; ipart < plist.size(); ++ipart) nMatchedHitsInPartList[ipart].resize(tjs.NumPlanes);
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
      for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
        if(hitTruTrkID[iht] == partList[ipart]->TrackId()) {
          ++nMatchedHitsInPartList[ipart][plane];
          if(tjs.fHits[iht].InTraj > 0) {
            unsigned short itj = tjs.fHits[iht].InTraj - 1;
            bool gotit = false;
            for(auto& hitInTj : nMatchedHitsInTj[ipart]) {
              if(hitInTj[0] == itj) {
                ++hitInTj[1];
                gotit = true;
              }
            } //  hitInTj
            if(!gotit) {
              std::array<unsigned short, 2> tmp {itj, 1};
              nMatchedHitsInTj[ipart].push_back(tmp);
            }
          } // inTraj > 0
        } // hit matched to partList
      } // ipart
    } // iht
    
    // remove partList elements that have no matched hits
    std::vector<simb::MCParticle*> newPartList;
    std::vector<std::vector<unsigned short>> newnMatchedHitsInPartList;
    std::vector<std::vector<std::array<unsigned short, 2>>> newnMatchedHitsInTj;
    for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
      unsigned short nht = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) nht += nMatchedHitsInPartList[ipart][plane];
      if(nht == 0) continue;
      newPartList.push_back(partList[ipart]);
      newnMatchedHitsInPartList.push_back(nMatchedHitsInPartList[ipart]);
      newnMatchedHitsInTj.push_back(nMatchedHitsInTj[ipart]);
    } // ipart
    partList = newPartList;
    nMatchedHitsInPartList = newnMatchedHitsInPartList;
    nMatchedHitsInTj = newnMatchedHitsInTj;
    
    // count the number of primary tracks that have at least 3 hits in at least 2 planes
    unsigned short nTruPrimary = 0;
    unsigned short nTruPrimaryOK = 0;
    for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
      if(partList[ipart]->Mother() != 0) continue;
      ++nTruPrimary;
      unsigned short nInPln = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(nMatchedHitsInPartList[ipart][plane] > 2) ++nInPln;
      } // plane
      if(nInPln > 1) ++nTruPrimaryOK;
    } // ipart
    
    neutrinoVxReconstructable = (neutrinoVxInFiducialVolume && nTruPrimaryOK > 1);
    
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
    
    if(tjs.MatchTruth[1] > 0) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Number of primary particles "<<nTruPrimary<<" Vtx reconstructable? "<<neutrinoVxReconstructable<<" Reconstructed? "<<neutrinoVxReconstructed<<" Correct? "<<neutrinoVxCorrect<<"\n";
      myprt<<"part   PDG   TrkID   MomID KE(MeV)   Process         Trajectory_extent_in_plane \n";
      for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
        unsigned short pdg = abs(partList[ipart]->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (partList[ipart]->E() - partList[ipart]->Mass());
        int motherID = partList[ipart]->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipart;
        myprt<<std::setw(6)<<partList[ipart]->PdgCode();
        myprt<<std::setw(10)<<partList[ipart]->TrackId();
        myprt<<std::setw(10)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<partList[ipart]->Process();
        // print the extent of the particle in each plane
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          unsigned int fht = UINT_MAX;
          unsigned int lht = 0;
          for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
            if(tjs.fHits[iht].WireID.Plane != plane) continue;
            unsigned short momTrackID = partList[ipart]->TrackId();
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
      } // ipart
    }
    
    // Declare a TJ - partlist match for the trajectory which has the most true hits.
    // another temp vector for the one-to-one match
    std::vector<std::vector<unsigned short>> partListToTjID(partList.size());
    for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) partListToTjID[ipart].resize(tjs.NumPlanes);
    
    for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(nMatchedHitsInPartList[ipart][plane] < 2) continue;
        unsigned short mostHits = 0;
        unsigned short tjWithMostHits = USHRT_MAX;
        for(unsigned short ii = 0; ii < nMatchedHitsInTj[ipart].size(); ++ii) {
          unsigned short itj = nMatchedHitsInTj[ipart][ii][0];
          geo::PlaneID planeID = DecodeCTP(tjs.allTraj[itj].CTP);
          // ensure we only check Tjs in the correct plane
          if(planeID.Plane != plane) continue;
          unsigned short nMatHits = nMatchedHitsInTj[ipart][ii][1];
          if(nMatHits > mostHits) {
            mostHits = nMatHits;
            tjWithMostHits = itj;
          }
        } // ii
        if(tjWithMostHits == USHRT_MAX) continue;
        // Count the total number of used hits in the TJ
        auto tmp = PutTrajHitsInVector(tjs.allTraj[tjWithMostHits], kUsedHits);
/* Sep20: This is a bad thing to do
        if(tjs.allTraj[tjWithMostHits].ParentID > 0) {
          // This is a daughter trajectory that has more truth matched hits than the parent
          unsigned short ptj = tjs.allTraj[tjWithMostHits].ParentID - 1;
          // add the parent hits to the vector of daughter hits
          auto ptmp = PutTrajHitsInVector(tjs.allTraj[ptj], kUsedHits);
          tmp.insert(tmp.end(), ptmp.begin(), ptmp.end());
          // revise the mostHits count. To do this we need to find the parent tj index in nMatchedHitsInTj
          for(unsigned short ii = 0; ii < nMatchedHitsInTj[ipart].size(); ++ii) {
            unsigned short itj = nMatchedHitsInTj[ipart][ii][0];
            if(itj == ptj) {
              mostHits += nMatchedHitsInTj[ipart][ii][1];
              // re-direct the calculation to the parent
              tjWithMostHits = ptj;
              break;
            } // found the parent tj
          } // ii
        } // deal with daughters
*/
        // count the number matched to a true particle
        float nTjHits = 0;
        for(auto& iht : tmp) if(hitTruTrkID[iht] > 0) ++nTjHits;
        float nTruHits = nMatchedHitsInPartList[ipart][plane];
        float nTjTruRecHits = mostHits;
        float eff = nTjTruRecHits / nTruHits;
        float pur = nTjTruRecHits / nTjHits;
        float effpur = eff * pur;
        // This overwrites any previous match that has poorer efficiency * purity
        if(effpur > tjs.allTraj[tjWithMostHits].EffPur) {
          tjs.allTraj[tjWithMostHits].MCPartListIndex = ipart;
          tjs.allTraj[tjWithMostHits].EffPur = effpur;
          partListToTjID[ipart][plane] = tjs.allTraj[tjWithMostHits].ID;
        }
      } // plane
    } // ipart
    
    // Update the EP sums
    for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
      float TMeV = 1000 * (partList[ipart]->E() - partList[ipart]->Mass());
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        // require at least 2 matched hits
        if(nMatchedHitsInPartList[ipart][plane] < 2) continue;
        unsigned short pdgIndex = PDGCodeIndex(tjs, partList[ipart]->PdgCode());
        // count the number of EP sums for this PDG code
        TSums[pdgIndex] += TMeV;
        ++EPCnts[pdgIndex];
        // find the first and last matched hit in this plane
        unsigned int fht = UINT_MAX;
        unsigned int lht = 0;
        // find the first and last matched hit in this plane
        for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
          if(tjs.fHits[iht].WireID.Plane != plane) continue;
          unsigned short momTrackID = partList[ipart]->TrackId();
          // Look for the real mother
          for(auto& md : moda) if(md.second == momTrackID) momTrackID = md.first;
          if(hitTruTrkID[iht] != momTrackID) continue;
          if(fht == UINT_MAX) fht = iht;
          lht = iht;
        } // iht
        if(fht == UINT_MAX) continue;
        if(partListToTjID[ipart][plane] == 0) {
          // Enter 0 in the profile histogram
          hist.fEP_T[pdgIndex]->Fill(TMeV, 0);
          if(nMatchedHitsInPartList[ipart][plane] > tjs.MatchTruth[3]) {
            ++nBadEP;
            mf::LogVerbatim myprt("TC");
            myprt<<"pdgIndex "<<pdgIndex<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to partList["<<ipart<<"]";
            myprt<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipart][plane];
            myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht])<<" events processed "<<tjs.EventsProcessed;
          }
          continue;
        }
        unsigned short itj = partListToTjID[ipart][plane] - 1;
        EPTSums[pdgIndex] += TMeV * tjs.allTraj[itj].EffPur;
        hist.fEP_T[pdgIndex]->Fill(TMeV, tjs.allTraj[itj].EffPur);
        // print out some debugging information if the EP was pitiful and the number of matched hits is large
        if(tjs.allTraj[itj].EffPur < tjs.MatchTruth[2] && nMatchedHitsInPartList[ipart][plane] > tjs.MatchTruth[3]) {
          mf::LogVerbatim myprt("TC");
          myprt<<"pdgIndex "<<pdgIndex<<" BadEP "<<std::fixed<<std::setprecision(2)<<tjs.allTraj[itj].EffPur;
          myprt<<" TMeV "<<(int)TMeV<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipart][plane];
          myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht])<<" events processed "<<tjs.EventsProcessed;
          // print alg names
          for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tjs.allTraj[itj].AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        }
        // check for a bad match to a primary electron shower
        if(tjs.ShowerTag[0] > 1 && pdgIndex == 0 && partList[ipart]->Mother() == 0 && fSourceParticleEnergy > 50) {
          Trajectory& ptj = tjs.allTraj[itj];
          // determine if this is a parent of a shower Tj
          int dtrID = 0;
          for(auto& tj : tjs.allTraj) {
            if(!tj.AlgMod[kShowerTj]) continue;
            if(tj.ParentID == ptj.ID) {
              dtrID = tj.ID;
              break;
            }
          } // tj
          if(dtrID == 0) mf::LogVerbatim("TC")<<"BadShower Primary electron -> Traj "<<ptj.ID<<"_"<<ptj.CTP<<" Wrong shower parent. Events processed "<<tjs.EventsProcessed;
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
    } // ipart
    
    // 
    
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
    
    // match PFParticles
    // initialize everything
    for(auto& pfp : tjs.pfps) pfp.MCPartListIndex = USHRT_MAX;
    
    // PFParticle reconstruction efficiency calculation
    // 1) Find the average EP for all Tjs in a PFParticle weighted by the Tj length, aveEP
    // 2) Accumulate MCP_EPTSum += aveEP * TMeV and MCP_TSum += TMeV in this function
    // 3) Calculate average PFParticle EP = MCP_EPTSum / TMeV in PrintResults()

    for(unsigned short ipart = 0; ipart < partList.size(); ++ipart) {
      auto& part = partList[ipart];
      float TMeV = 1000 * (part->E() - part->Mass());
      // skip low energy electrons
      int pdg = abs(part->PdgCode());
      if(pdg == 11 && TMeV < 100) continue;
      unsigned short nInPln = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        // require at least 2 matched hits
        if(nMatchedHitsInPartList[ipart][plane] < 2) continue;
        ++nInPln;
      } // plane
      // require matched hits in at least two planes
      if(nInPln < 2) continue;
      MCP_TSum += TMeV;
      ++MCP_Cnt;
      bool gotit = false;
      for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
        auto& pfp = tjs.pfps[ipfp];
        if(pfp.ID == 0) continue;
        unsigned short cnt = 0;
        float sum = 0;
        float epsum = 0;
        for(auto& tjID : pfp.TjIDs) {
          unsigned short itj = tjID - 1;
          Trajectory& tj = tjs.allTraj[itj];
          if(tj.MCPartListIndex == ipart) {
            ++cnt;
            float npts = NumPtsWithCharge(tjs, tj, false);
            sum += npts;
            epsum += npts * tj.EffPur;
          }
        } // tjID
        // require at least 2 Tjs match this PFParticle
        if(cnt > 1 && sum > 0) {
          float pfpEP = epsum / sum;
          pfp.EffPur = pfpEP;
          MCP_EPTSum += TMeV * pfpEP;
          pfp.MCPartListIndex = ipart;
          ++PFP_CntGoodMat;
          gotit = true;

          if (fStudyMode) {
            //
            if (abs(part->PdgCode())==11) continue;
            //
            bool wrongid = part->PdgCode()!=pfp.PDGCode;
            //
            hist.hasfit_match->Fill(pfp.Track.ID()>=0);
            if (wrongid) {
              hist.hasfit_wrongid->Fill(pfp.Track.ID()>=0);
            } else {
              hist.hasfit_okid->Fill(pfp.Track.ID()>=0);
            }
            //
            if (pfp.Track.ID()<0) continue;
            // std::cout << "pfp #" << ipfp << " pdg=" << pfp.PDGCode << " tj vtx=" << recob::tracking::Point_t(pfp.XYZ[0][0],pfp.XYZ[0][1],pfp.XYZ[0][2]) << " dir=" << recob::tracking::Vector_t(pfp.Dir[0][0],pfp.Dir[0][1],pfp.Dir[0][2]) << " fit vtx=" << pfp.Track.Start() << " dir=" << pfp.Track.StartDirection() << " match part #" << ipart << " vtx=" << recob::tracking::Point_t(part->Vx(), part->Vy(), part->Vz()) << " dir=" << recob::tracking::Vector_t(part->Px()/part->P(), part->Py()/part->P(), part->Pz()/part->P()) << " mom=" << part->P() << std::endl;
            hist.nvalidpoints_match->Fill(pfp.Track.CountValidPoints());
            hist.nrejectpoints_match->Fill(pfp.Track.NPoints()-pfp.Track.CountValidPoints());
            if (pfp.Track.CountValidPoints()>1) {
              //
              detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
              double g4Ticks = detClocks->TPCG4Time2Tick(part->T())+tjs.detprop->GetXTicksOffset(0,0,0)-tjs.detprop->TriggerOffset();
              double xOffset = tjs.detprop->ConvertTicksToX(g4Ticks, 0, 0, 0);
              //
              trkf::TrackStatePropagator prop(1.0, 0.1, 10, 10., 0.01, false);
              recob::tracking::Plane mcplane(recob::tracking::Point_t(part->Vx()+xOffset, part->Vy(), part->Vz()),
                                             recob::tracking::Vector_t(part->Px()/part->P(), part->Py()/part->P(), part->Pz()/part->P()));
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
              hist.dXkf_match->Fill(tkatmc.position().X()-part->Vx()-xOffset);
              hist.dXtc_match->Fill(tjatmc.position().X()-part->Vx()-xOffset);
              hist.dYkf_match->Fill(tkatmc.position().Y()-part->Vy());
              hist.dYtc_match->Fill(tjatmc.position().Y()-part->Vy());
              hist.dZkf_match->Fill(tkatmc.position().Z()-part->Vz());
              hist.dZtc_match->Fill(tjatmc.position().Z()-part->Vz());
              hist.dUXkf_match->Fill(tkatmc.momentum().X()/recmom-part->Px()/part->P());
              hist.dUXtc_match->Fill(tjatmc.momentum().Unit().X()-part->Px()/part->P());
              hist.dUYkf_match->Fill(tkatmc.momentum().Y()/recmom-part->Py()/part->P());
              hist.dUYtc_match->Fill(tjatmc.momentum().Unit().Y()-part->Py()/part->P());
              hist.dUZkf_match->Fill(tkatmc.momentum().Z()/recmom-part->Pz()/part->P());
              hist.dUZtc_match->Fill(tjatmc.momentum().Unit().Z()-part->Pz()/part->P());
              hist.dXpull_match->Fill( (tkatmc.position().X()-part->Vx()-xOffset)/ sqrt(c(0,0)) );
              hist.dYpull_match->Fill( (tkatmc.position().Y()-part->Vy())/ sqrt(c(1,1)) );
              hist.dZpull_match->Fill( (tkatmc.position().Z()-part->Vz())/ sqrt(c(2,2)) );
              hist.dUXpull_match->Fill( (tkatmc.momentum().X()/recmom-part->Px()/part->P())/ sqrt(c(3,3)) );
              hist.dUYpull_match->Fill( (tkatmc.momentum().Y()/recmom-part->Py()/part->P())/ sqrt(c(4,4)) );
              hist.dUZpull_match->Fill( (tkatmc.momentum().Z()/recmom-part->Pz()/part->P())/ sqrt(c(5,5)) );
              hist.nchi2_match->Fill( pfp.Track.Chi2PerNdof() );
              hist.covtrace_match->Fill( cv(0,0)+cv(1,1)+cv(2,2)+cv(3,3) );
              if (wrongid) {
                hist.nchi2_wrongid->Fill( pfp.Track.Chi2PerNdof() );
                hist.dXkf_wrongid->Fill(tkatmc.position().X()-part->Vx()-xOffset);
                hist.dXtc_wrongid->Fill(tjatmc.position().X()-part->Vx()-xOffset);
                hist.dYkf_wrongid->Fill(tkatmc.position().Y()-part->Vy());
                hist.dYtc_wrongid->Fill(tjatmc.position().Y()-part->Vy());
                hist.dZkf_wrongid->Fill(tkatmc.position().Z()-part->Vz());
                hist.dZtc_wrongid->Fill(tjatmc.position().Z()-part->Vz());
                hist.dUXkf_wrongid->Fill(tkatmc.momentum().X()/recmom-part->Px()/part->P());
                hist.dUXtc_wrongid->Fill(tjatmc.momentum().Unit().X()-part->Px()/part->P());
                hist.dUYkf_wrongid->Fill(tkatmc.momentum().Y()/recmom-part->Py()/part->P());
                hist.dUYtc_wrongid->Fill(tjatmc.momentum().Unit().Y()-part->Py()/part->P());
                hist.dUZkf_wrongid->Fill(tkatmc.momentum().Z()/recmom-part->Pz()/part->P());
                hist.dUZtc_wrongid->Fill(tjatmc.momentum().Unit().Z()-part->Pz()/part->P());
              } else {
                hist.dXkf_okid->Fill(tkatmc.position().X()-part->Vx()-xOffset);
                hist.dXtc_okid->Fill(tjatmc.position().X()-part->Vx()-xOffset);
                hist.dYkf_okid->Fill(tkatmc.position().Y()-part->Vy());
                hist.dYtc_okid->Fill(tjatmc.position().Y()-part->Vy());
                hist.dZkf_okid->Fill(tkatmc.position().Z()-part->Vz());
                hist.dZtc_okid->Fill(tjatmc.position().Z()-part->Vz());
                hist.dUXkf_okid->Fill(tkatmc.momentum().X()/recmom-part->Px()/part->P());
                hist.dUXtc_okid->Fill(tjatmc.momentum().Unit().X()-part->Px()/part->P());
                hist.dUYkf_okid->Fill(tkatmc.momentum().Y()/recmom-part->Py()/part->P());
                hist.dUYtc_okid->Fill(tjatmc.momentum().Unit().Y()-part->Py()/part->P());
                hist.dUZkf_okid->Fill(tkatmc.momentum().Z()/recmom-part->Pz()/part->P());
                hist.dUZtc_okid->Fill(tjatmc.momentum().Unit().Z()-part->Pz()/part->P());
                hist.dXpull_okid->Fill( (tkatmc.position().X()-part->Vx()-xOffset)/ sqrt(c(0,0)) );
                hist.dYpull_okid->Fill( (tkatmc.position().Y()-part->Vy())/ sqrt(c(1,1)) );
                hist.dZpull_okid->Fill( (tkatmc.position().Z()-part->Vz())/ sqrt(c(2,2)) );
                hist.dUXpull_okid->Fill( (tkatmc.momentum().X()/recmom-part->Px()/part->P())/ sqrt(c(3,3)) );
                hist.dUYpull_okid->Fill( (tkatmc.momentum().Y()/recmom-part->Py()/part->P())/ sqrt(c(4,4)) );
                hist.dUZpull_okid->Fill( (tkatmc.momentum().Z()/recmom-part->Pz()/part->P())/ sqrt(c(5,5)) );
                hist.nchi2_okid->Fill( pfp.Track.Chi2PerNdof() );
              }
              //
              // if (std::abs(tkatmc.position().X()-part->Vx())>10. && std::abs(tjatmc.position().X()-part->Vx())<1.) {
              //   std::cout << "DEBUG ME" << std::endl;
              // }
              //
            }
          } //fStudyMode
          break;
        }
      } // ipfp
      
      if(!gotit && TMeV > 30) {
        mf::LogVerbatim myprt("TC");
        myprt<<"BadPFP PDGCode "<<part->PdgCode()<<" TMeV "<<(int)TMeV;
        myprt<<" nMatchedHitsInPartList ";
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) myprt<<" "<<nMatchedHitsInPartList[ipart][plane];
        myprt<<" matched Tjs ";
        for(auto& tj : tjs.allTraj) {
          if(tj.AlgMod[kKilled]) continue;
          if(tj.MCPartListIndex == ipart) myprt<<" "<<tj.ID<<" EP "<<std::fixed<<std::setprecision(2)<<tj.EffPur;
        } // tj
        myprt<<" events processed "<<tjs.EventsProcessed;
      }
    } // ipart
    
    // histogram reconstructed PDG code vs true PDG code
    std::array<int, 5> recoCodeList = {0, 11, 13, 211, 2212};
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      // require match to MC
      if(pfp.MCPartListIndex == USHRT_MAX) continue;
      short truIndex = PDGCodeIndex(tjs, partList[pfp.MCPartListIndex]->PdgCode());
      short recIndex = 0;
      for(recIndex = 0; recIndex < 5; ++recIndex) if(pfp.PDGCode == recoCodeList[recIndex]) break;
      if(recIndex == 5) {
        std::cout<<"MT: Found an unknown PDGCode "<<pfp.PDGCode<<" in PFParticle "<<pfp.ID<<"\n";
        continue;
      }
      std::cout<<"PFP "<<pfp.ID<<" truIndex "<<truIndex<<" PDGCode "<<partList[pfp.MCPartListIndex]->PdgCode()<<" recIndex "<<recIndex<<" recPDG "<<pfp.PDGCode<<"\n";
      hist.PDGCode_reco_true->Fill((float)truIndex, (float)recIndex);
    } // pfp


    if (fStudyMode) {
      // nomatch
      for(unsigned short ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
        auto& pfp = tjs.pfps[ipfp];
        if(pfp.ID == 0) continue;
        if (pfp.MCPartListIndex!=USHRT_MAX) continue;
        //
        hist.hasfit_nomatch->Fill(pfp.Track.ID()>=0);
        if (pfp.Track.ID()<0) continue;
        hist.nvalidpoints_nomatch->Fill(pfp.Track.CountValidPoints());
	hist.nrejectpoints_nomatch->Fill(pfp.Track.NPoints()-pfp.Track.CountValidPoints());
        if (pfp.Track.CountValidPoints()>1) {
          hist.nchi2_nomatch->Fill( pfp.Track.Chi2PerNdof() );
          auto cv = pfp.Track.VertexCovarianceLocal5D();
          hist.covtrace_nomatch->Fill( cv(0,0)+cv(1,1)+cv(2,2)+cv(3,3) );
        }
      }
    } //fStudyMode
    
    // update the total PFParticle count
    for(auto& pfp : tjs.pfps) if(pfp.ID > 0) ++PFP_Cnt;
    
  } // MatchTruth
  
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
      float nofrac = 1 - (PFP_CntGoodMat / PFP_Cnt);
      myprt<<" PFP "<<ep<<" MCP cnt "<<(int)MCP_Cnt<<" PFP cnt "<<(int)PFP_Cnt<<" noMatFrac "<<nofrac;
    }
    myprt<<" VxCount";
    for(auto cnt : TruVxCounts) myprt<<" "<<cnt;
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
