//############################################################################
//### Name:        ShowerTrackFinderCheater                                ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/TRACSAlg.h"
#include "larreco/RecoAlg/TRACSCheatingAlg.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//C++ Includes
#include <iostream>

//Root Includes
#include "TMath.h"
#include "TVector.h"

namespace ShowerRecoTools {

  class ShowerTrackFinderCheater:IShowerTool {

    public:

      ShowerTrackFinderCheater(const fhicl::ParameterSet& pset);

      ~ShowerTrackFinderCheater();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			   art::Event& Event, 
			   reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    private:

      //Algorithm functions
      shower::TRACSCheatingAlg fTRACSCheatingAlg;


      //fcl
      bool fDebugEVD;
      art::InputTag fPFParticleModuleLabel;
      art::InputTag fHitModuleLabel;

      std::string fTrueParticleIntputLabel;
      std::string fShowerStartPositionInputTag;
      std::string fShowerDirectionInputTag;
      std::string fInitialTrackHitsOutputLabel;
      std::string fInitialTrackSpacePointsOutputLabel;
  };


  ShowerTrackFinderCheater::ShowerTrackFinderCheater(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fTRACSCheatingAlg(pset.get<fhicl::ParameterSet>("TRACSCheatingAlg")),
    fDebugEVD(pset.get<bool>("DebugEVD")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel")),
    fTrueParticleIntputLabel(pset.get<std::string>("TrueParticleIntputLabel")),
    fShowerStartPositionInputTag(pset.get<std::string>("ShowerStartPositionInputTag")),
    fShowerDirectionInputTag(pset.get<std::string>("ShowerDirectionInputTag")),
    fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel")),
    fInitialTrackSpacePointsOutputLabel(pset.get<std::string>("InitialTrackSpacePointsOutputLabel"))
  {
  }

  ShowerTrackFinderCheater::~ShowerTrackFinderCheater()
  {
  }

  int ShowerTrackFinderCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, 
						 art::Event& Event, 
						 reco::shower::ShowerElementHolder& ShowerEleHolder){

    const simb::MCParticle* trueParticle;

    //Get the hits from the shower:
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerTrackFinderCheater") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    if (ShowerEleHolder.CheckElement(fTrueParticleIntputLabel)){
      ShowerEleHolder.GetElement(fTrueParticleIntputLabel,trueParticle);
    } else {

      //Could store these in the shower element holder and just calculate once?
      std::map<int,const simb::MCParticle*> trueParticles = fTRACSCheatingAlg.GetTrueParticleMap();
      std::map<int,std::vector<int> > showersMothers = fTRACSCheatingAlg.GetTrueChain(trueParticles);


      //Get the clusters
      art::Handle<std::vector<recob::Cluster> > clusHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
        throw cet::exception("ShowerTrackFinderCheater") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
        return 1;
      }
      art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
      std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

      //Get the hit association
      art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

      std::vector<art::Ptr<recob::Hit> > showerHits;
      for(auto const& cluster: clusters){

        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
        showerHits.insert(showerHits.end(),hits.begin(),hits.end());
      }

      //Get the true particle from the shower
      std::pair<int,double> ShowerTrackInfo = fTRACSCheatingAlg.TrueParticleIDFromTrueChain(showersMothers,showerHits,2);

      if(ShowerTrackInfo.first==-99999) {
        mf::LogError("ShowerStartPosition") << "True Shower Not Found";
        return 1;
      }
      trueParticle = trueParticles[ShowerTrackInfo.first];
    }

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputTag)){
      mf::LogError("ShowerTrackFinderCheater") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputTag)){
      mf::LogError("ShowerTrackFinderCheater") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputTag,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputTag,ShowerDirection);

    art::Handle<std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hits;
    if(Event.getByLabel(fHitModuleLabel, hitHandle)){
      art::fill_ptr_vector(hits, hitHandle);
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::SpacePoint> fmsph(hitHandle, Event, fPFParticleModuleLabel);
    if(!fmsph.isValid()){
      throw cet::exception("ShowerTrackFinderCheater") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    std::vector<art::Ptr<recob::Hit> > trackHits;
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

    //Get the hits from the true particle
    for (auto hit : hits){
      int trueParticleID = fTRACSCheatingAlg.TrueParticleID(hit);
      if (trueParticleID == trueParticle->TrackId()){
        trackHits.push_back(hit);
        std::vector<art::Ptr<recob::SpacePoint> > sps = fmsph.at(hit.key());
        if (sps.size() == 1){
          trackSpacePoints.push_back(sps.front());
        }
      }
    }

    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOutputLabel);
    ShowerEleHolder.SetElement(trackSpacePoints,fInitialTrackSpacePointsOutputLabel);

    if (fDebugEVD){
      fTRACSCheatingAlg.CheatDebugEVD(trueParticle, Event, ShowerEleHolder, pfparticle);
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinderCheater)
