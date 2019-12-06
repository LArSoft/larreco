//############################################################################
//### Name:        ShowerPMATrackFinder                                    ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for creating the initial track of the shower       ###
//###              using PMA. Derived form EMShower credit Mike Wallbank.  ###
//###                                                                      ###
//############################################################################
#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/TRACSAlg.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

//C++ Includes
#include <iostream>
//#include <math.h>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools{

  class ShowerPMATrackFinder:IShowerTool {
  public:

    ShowerPMATrackFinder(const fhicl::ParameterSet& pset);

    //Generic Track Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;

  private:

    void InitialiseProducers() override;

    //Function to add the assoctions
    int AddAssociations(art::Event& Event,
			reco::shower::ShowerElementHolder& ShowerEleHolder) override;



    // Define algoritms and services.
    pma::ProjectionMatchingAlg        fProjectionMatchingAlg;
    art::ServiceHandle<geo::Geometry> fGeom;

    //fcl paramters
    float fMinTrajectoryPoints; //Minimum number of trajectory points returned from the fit to decide the track is good.
    std::string fInitialTrackLengthOutputLabel;
    std::string fInitialTrackOutputLabel;
    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fInitialTrackHitsInputLabel;
  };


  ShowerPMATrackFinder::ShowerPMATrackFinder(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")),
    fMinTrajectoryPoints(pset.get<float>("MinTrajectoryPoints")),
    fInitialTrackLengthOutputLabel(pset.get<std::string>("InitialTrackLengthOutputLabel")),
    fInitialTrackOutputLabel(pset.get<std::string>("InitialTrackOutputLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel")),
    fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel"))
  {
  }

  void ShowerPMATrackFinder::InitialiseProducers(){
    //Set up the recob::Track and the Assns so they can be put in the event
    InitialiseProduct<std::vector<recob::Track> >(fInitialTrackOutputLabel);
    InitialiseProduct<art::Assns<recob::Shower, recob::Track > >("ShowerTrackAssn");
    InitialiseProduct<art::Assns<recob::Track, recob::Hit > >("ShowerTrackHitAssn");

  }



  int ShowerPMATrackFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      mf::LogError("ShowerPMATrackFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
      mf::LogError("ShowerPMATrackFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fInitialTrackHitsInputLabel)){
      mf::LogError("ShowerPMATrackFinder") << "Initial track hits are not set, returning "<< std::endl;
      return 1;
    }


    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    std::vector<art::Ptr<recob::Hit> > InitialTrackHits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel,InitialTrackHits);

    //Get the hits in term of planes.
    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_trackhits;
    for(auto const& hit: InitialTrackHits){
      geo::WireID wire = hit->WireID();
      geo::PlaneID plane = wire.asPlaneID();
      plane_trackhits[plane].push_back(hit);
    }

    //Decide which plane to use
    int maxhits       = 1; // set to 1 as we require at least 2 track hits per plane
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition);
    geo::PlaneID maxplane; // Note this contains both tpc and cryostat information

    //Find the plane with the max number of hits in the tpc where the vertex is.
    for(auto const& plane : plane_trackhits){
      std::vector<art::Ptr<recob::Hit> >  trackhits = plane.second;
      geo::TPCID maxTPC = (plane.first).asTPCID();
      if( maxTPC == vtxTPC){
        if((int) trackhits.size() > maxhits ){
          maxplane = plane.first;
          maxhits  = trackhits.size();
        }
      }
    }

    if( maxhits == 1 || !maxplane){
      mf::LogError("ShowerPMATrackFinder") << "Max Plane not set " << std::endl;
      return 1;
    }

    int nextmaxhits  = 1;
    geo::PlaneID nextmaxplane;
    //Find the next largest plane.
    for(auto const& plane : plane_trackhits){
      //Check clusters are not in same plane
      if( (plane.first) == maxplane){continue;}
      //Need to make sure clusters are in same tpc
      geo::TPCID nextmaxTPC = (plane.first).asTPCID();
      if( nextmaxTPC == vtxTPC){
        std::vector<art::Ptr<recob::Hit> > trackhits = plane.second;
        if((int) trackhits.size() > nextmaxhits){
          nextmaxplane = plane.first;
          nextmaxhits  = trackhits.size();
        }
      }
    }


    if( nextmaxhits == 1 || !nextmaxplane){
      mf::LogError("ShowerPMATrackFinder") << "Next Max Plane not set " << std::endl;
      return 1;
    }

    std::vector<art::Ptr<recob::Hit> > maxPlaneHits = (plane_trackhits.at(maxplane));
    std::vector<art::Ptr<recob::Hit> > nextmaxPlaneHits = (plane_trackhits.at(nextmaxplane));

    if(maxPlaneHits.size() < 2 && nextmaxPlaneHits.size() < 2){
      mf::LogError("ShowerPMATrackFinder") << "Not enough hit to make track " << std::endl;
      return 0;
    }

    //Build the 3D track
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(maxPlaneHits, nextmaxPlaneHits, ShowerStartPosition);

    if(!pmatrack){
      mf::LogError("ShowerPMATrackFinder") << "Failed fit " << std::endl;
      return 1;
    }

    //Get the spacepoints
    std::vector<TVector3> spts;
    for (size_t i = 0; i<pmatrack->size(); ++i){
      if ((*pmatrack)[i]->IsEnabled()){
        TVector3 p3d = (*pmatrack)[i]->Point3D();
        spts.push_back(p3d);
      }
    }

    if(spts.size() < fMinTrajectoryPoints){
      mf::LogWarning("ShowerPMATrackFinder") << "Not Enough Spacepoints" << std::endl;
      return 1;
    }

    //Make the track. Adapted from PandoraShowerCreation
    recob::tracking::Positions_t xyz;
    recob::tracking::Momenta_t pxpypz;
    recob::TrackTrajectory::Flags_t flags;

    TVector3 spt1 = spts[0];

    //Loop over the trajectory points and start making the track
    for(unsigned int sp=0; sp<spts.size(); ++sp){

      TVector3 spt = spts[sp];

      if(sp < spts.size() -1){
        spt1 = spts[sp+1];
      }
      else{
        spt1 = spts[sp];
      }

      //Make the xyz for the trajectory point
      xyz.emplace_back(recob::tracking::Point_t(spt.X(), spt.Y(), spt.Z()));

      TVector3 dir = -(spt - spt1);

      //Make the direction at point ofr the trajectory point
      pxpypz.emplace_back(recob::tracking::Vector_t(dir.X(), dir.Y(), dir.Z()));

      if (std::fabs(spt.X()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
          std::fabs(spt.X()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
          std::fabs(spt.X()-util::kBogusF)<std::numeric_limits<float>::epsilon())
      {
        flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex, recob::TrajectoryPointFlagTraits::NoPoint));
      }
      else {
        flags.emplace_back(recob::TrajectoryPointFlags());
      }
      spt1 = spt;
    }

    std::vector<art::Ptr<recob::Hit> > TrackHits;
    for(auto const& trackhits_p: plane_trackhits){
      TrackHits.insert(TrackHits.end(),trackhits_p.second.begin(),trackhits_p.second.end());
    }

    //Actually make the thing.
    recob::Track track = recob::Track(recob::TrackTrajectory(std::move(xyz), std::move(pxpypz),std::move(flags), false),
				      util::kBogusI, util::kBogusF, util::kBogusI,
				      recob::tracking::SMatrixSym55(), recob::tracking::SMatrixSym55(), pfparticle.key());
    
    ShowerEleHolder.SetElement(track,fInitialTrackOutputLabel);

    TVector3 Start = {track.Start().X(), track.Start().Y(), track.Start().Z()};
    TVector3 End   = {track.End().X(), track.End().Y(),track.End().Z()};
    float tracklength = (Start-End).Mag();

    ShowerEleHolder.SetElement(tracklength,fInitialTrackLengthOutputLabel);

    return 0;
  }


  int ShowerPMATrackFinder::AddAssociations(art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    //Check the track has been set
    if(!ShowerEleHolder.CheckElement(fInitialTrackOutputLabel)){
      mf::LogError("ShowerPMATrackFinderAddAssn") << "Track not set so the assocation can not be made  "<< std::endl;
      return 1;
    }

    //Get the size of the ptr as it is.
    int trackptrsize = GetVectorPtrSize(fInitialTrackOutputLabel);

    const art::Ptr<recob::Track> trackptr = GetProducedElementPtr<recob::Track>(fInitialTrackOutputLabel, ShowerEleHolder,trackptrsize-1);
    const art::Ptr<recob::Shower> showerptr = GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::Track> >(showerptr,trackptr,"ShowerTrackAssn");

    std::vector<art::Ptr<recob::Hit> > TrackHits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel,TrackHits);

    for(auto const& TrackHit: TrackHits){
      AddSingle<art::Assns<recob::Track, recob::Hit> >(trackptr,TrackHit,"ShowerTrackHitAssn");
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPMATrackFinder)
