//############################################################################
//### Name:        ShowerDirectionCheater                                  ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/TRACSCheatingAlg.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//C++ Includes
#include <iostream>

//Root Includes
#include "TMath.h"
#include "TTree.h"

namespace ShowerRecoTools {

  class ShowerDirectionCheater:IShowerTool {

  public:

    ShowerDirectionCheater(const fhicl::ParameterSet& pset);

    ~ShowerDirectionCheater();

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;

  private:

    TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints_pfp, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre);
    double RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction);
    double CalculateRMS(std::vector<float> perps);

    //Algorithm functions
    shower::TRACSCheatingAlg fTRACSCheatingAlg;

    //Services
    art::ServiceHandle<art::TFileService> tfs;

    //fcl
    art::InputTag fPFParticleModuleLabel;
    float fNSegments; //Number of segement to split the shower into the perforam the RMSFlip.
    bool fRMSFlip;    //Flip the direction by considering the rms.
    bool fVertexFlip; //Flip the direction by considering the vertex position relative to the center position.

    //TTree Branch variables
    TTree* Tree;
    float vertexDotProduct;
    float rmsGradient;

    std::string fShowerStartPositionInputLabel;
    std::string fTrueParticleInputLabel;
    std::string fShowerDirectionOuputLabel;

  };


  ShowerDirectionCheater::ShowerDirectionCheater(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fTRACSCheatingAlg(pset.get<fhicl::ParameterSet>("TRACSCheatingAlg")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fNSegments(pset.get<float>("NSegments")),
    fRMSFlip(pset.get<bool>("RMSFlip")),
    fVertexFlip(pset.get<bool>("VertexFlip")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fTrueParticleInputLabel(pset.get<std::string>("TrueParticleInputLabel")),
    fShowerDirectionOuputLabel(pset.get<std::string>("ShowerDirectionOuputLabel"))
  {
    if (vertexDotProduct||rmsGradient){
      Tree = tfs->make<TTree>("DebugTreeDirCheater", "DebugTree from shower direction cheater");
      if (fVertexFlip) Tree->Branch("vertexDotProduct",&vertexDotProduct);
      if (fRMSFlip)    Tree->Branch("rmsGradient",&rmsGradient);
    }
  }

  ShowerDirectionCheater::~ShowerDirectionCheater()
  {
  }

  int ShowerDirectionCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, 
					       art::Event& Event, 
					       reco::shower::ShowerElementHolder& ShowerEleHolder){


    const simb::MCParticle* trueParticle;

    //Get the hits from the shower:
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerDirectionCheater") 
	<< "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    if (ShowerEleHolder.CheckElement(fTrueParticleInputLabel)){
      ShowerEleHolder.GetElement(fTrueParticleInputLabel,trueParticle);
    } else {

      //Could store these in the shower element holder and just calculate once?
      std::map<int,const simb::MCParticle*> trueParticles = fTRACSCheatingAlg.GetTrueParticleMap();
      std::map<int,std::vector<int> > showersMothers = fTRACSCheatingAlg.GetTrueChain(trueParticles);

      //Get the clusters
      art::Handle<std::vector<recob::Cluster> > clusHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
        throw cet::exception("ShowerDirectionCheater") 
	  << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
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
        mf::LogError("ShowerDirectionCheater") << "True shower not found, returning";
        return 1;
      }
      trueParticle = trueParticles[ShowerTrackInfo.first];
    }

    TVector3 trueDir = {trueParticle->Px(),trueParticle->Py(),trueParticle->Pz()};
    trueDir = trueDir.Unit(); // TODO: Can probably remove?

    TVector3 trueDirErr = {-999,-999,-999};
    ShowerEleHolder.SetElement(trueDir,trueDirErr,fShowerDirectionOuputLabel);

    if (fRMSFlip || fVertexFlip){
      //Get the SpacePoints and hits
      art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);

      if (!fmspp.isValid()){
        throw cet::exception("ShowerDirectionCheater") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
        return 1;
      }

      art::Handle<std::vector<recob::SpacePoint> > spHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
        throw cet::exception("ShowerDirectionCheater") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
        return 1;
      }
      art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
      if(!fmh.isValid()){
        throw cet::exception("ShowerDirectionCheater") << "Spacepoint and hit association not valid. Stopping.";
        return 1;
      }
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

      if (spacePoints.size()==0) return 1;

      //Get Shower Centre
      float TotalCharge;
      TVector3 ShowerCentre = IShowerTool::GetTRACSAlg().ShowerCentre(spacePoints, fmh, TotalCharge);

      //Check if we are pointing the correct direction or not, First try the start position
      if(ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel) && fVertexFlip){

        //Get the General direction as the vector between the start position and the centre
        TVector3 StartPositionVec = {-999, -999, -999};
        ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,StartPositionVec);

        TVector3 GeneralDir       = (ShowerCentre - StartPositionVec).Unit();

        //Dot product
        vertexDotProduct = trueDir.Dot(GeneralDir);

        //If the dotproduct is negative the Direction needs Flipping
        if(vertexDotProduct < 0){
          trueDir[0] = - trueDir[0];
          trueDir[1] = - trueDir[1];
          trueDir[2] = - trueDir[2];
        }

        //ShowerEleHolder.SetShowerDirection(trueDir);
      }

      if (fRMSFlip){
        //Otherwise Check against the RMS of the shower. Method adapated from EMShower Thanks Mike.
        rmsGradient = RMSShowerGradient(spacePoints,ShowerCentre,trueDir);
        if(rmsGradient < 0){

          trueDir[0] = - trueDir[0];
          trueDir[1] = - trueDir[1];
          trueDir[2] = - trueDir[2];
        }

        //ShowerEleHolder.SetShowerDirection(trueDir);
      }
      Tree->Fill();
    }

    return 0;
  }


  //Function to calculate the RMS at segements of the shower and calculate the gradient of this. If negative then the direction is pointing the opposite way to the correct one
  double ShowerDirectionCheater::RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction){

    //Order the spacepoints
    IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(sps,ShowerCentre,Direction);

    //Get the length of the shower.
    double minProj =IShowerTool::GetTRACSAlg().SpacePointProjection(sps[0],ShowerCentre,Direction);
    double maxProj =IShowerTool::GetTRACSAlg().SpacePointProjection(sps[sps.size()-1],ShowerCentre,Direction);

    double length = (maxProj-minProj);
    double segmentsize = length/fNSegments;

    std::map<int, std::vector<float> > len_segment_map;

    //Split the the spacepoints into segments.
    for(auto const& sp: sps){

      //Get the the projected length
      double len = IShowerTool::GetTRACSAlg().SpacePointProjection(sp,ShowerCentre,Direction);

      //Get the length to the projection
      double  len_perp = IShowerTool::GetTRACSAlg().SpacePointPerpendiular(sp,ShowerCentre,Direction,len);

      int sg_len = round(len/segmentsize);
      //TODO: look at this:
      //int sg_len = round(len/segmentsize+fNSegments/2); //Add to make positive
      len_segment_map[sg_len].push_back(len_perp);
    }

    int counter = 0;
    float sumx  = 0;
    float sumy  = 0;
    float sumx2 = 0;
    float sumxy = 0;

    //Get the rms of the segments and caclulate the gradient.
    for(auto const& segment: len_segment_map){
      if (segment.second.size()<2) continue;
      float RMS = CalculateRMS(segment.second);
      //Calculate the gradient using regression
      sumx  += segment.first;
      sumy  += RMS;
      sumx2 += segment.first * segment.first;
      sumxy += RMS * segment.first;
      ++counter;
    }

    double RMSgradient = (counter*sumxy - sumx*sumy)/(counter*sumx2 - sumx*sumx);
    return RMSgradient;

  }

  double ShowerDirectionCheater::CalculateRMS(std::vector<float> perps){
    int counter = 0;
    double sum  = 0;
    for (const auto &perp : perps){
      sum= perp*perp;
      ++counter;
    }
    double rms = TMath::Sqrt(sum/(counter-1));

    return rms;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerDirectionCheater)
