//############################################################################
//### Name:        ShowerSmartTrackTrajectoryPointDirection                ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              first trajectory of the initial track                   ###
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
#include "larreco/RecoAlg/TRACSAlg.h"

//C++ Includes
#include <iostream>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools {


  class ShowerSmartTrackTrajectoryPointDirection:IShowerTool {

  public:

    ShowerSmartTrackTrajectoryPointDirection(const fhicl::ParameterSet& pset);

    ~ShowerSmartTrackTrajectoryPointDirection();

    //Calculate the direction from the inital track.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;

  private:

    //fcl
    bool  fUsePandoraVertex;    //Direction from point defined as 
                                //(Position of traj point - Vertex) rather than 
                                //(Position of traj point - Track Start Point). 
    bool  fAllowDynamicSliding; //Rather than evualte the angle from the start use 
                                //the previous trajectory point position.
    bool  fUsePositionInfo;     //Don't use the DirectionAtPoint rather than 
                                //definition above.
                                //((Position of traj point + 1)-(Position of traj point)
    bool  fUseStartPos;         //Rather the using the angles between the directions
                                //from start position to the trajectory points
                                //use the angle between the the points themselves
    float fAngleCut;
    
    std::string fInitialTrackInputLabel;
    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionOuputLabel;
  };


  ShowerSmartTrackTrajectoryPointDirection::ShowerSmartTrackTrajectoryPointDirection(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fUsePandoraVertex(pset.get<bool>("UsePandoraVertex")),
    fAllowDynamicSliding(pset.get<bool>("AllowDynamicSliding")),
    fUsePositionInfo(pset.get<bool>("UsePositionInfo")),
    fUseStartPos(pset.get<bool>("UseStartPos")),
    fAngleCut(pset.get<float>("AngleCut")),
    fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionOuputLabel(pset.get<std::string>("ShowerDirectionOuputLabel"))
    
  {
  }

  ShowerSmartTrackTrajectoryPointDirection::~ShowerSmartTrackTrajectoryPointDirection()
  {
  }

  int ShowerSmartTrackTrajectoryPointDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Check the Track has been defined
    if(!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)){
      mf::LogError("ShowerSmartTrackTrajectoryPointDirection")
	<< "Initial track not set"<< std::endl;
      return 1;
    }
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement(fInitialTrackInputLabel,InitialTrack);

    //Smartly choose the which trajectory point to look at by ignoring the smush of hits at the vertex.
    if(InitialTrack.NumberTrajectoryPoints() == 1){
      mf::LogError("ShowerSmartTrackTrajectoryPointDirection") 
	<< "Not Enough trajectory points."<< std::endl;
      return 1;
    }

    //Trajectory point which the direction is calcualted for.
    int trajpoint = 0;
    geo::Vector_t Direction_vec;

    if(fUsePositionInfo){
      //Get the start position.
      geo::Point_t StartPosition;

      if(fUsePandoraVertex){
        //Check the Track has been defined
        if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
          mf::LogError("ShowerSmartTrackTrajectoryPointDirection")
	    << "Shower start position not set"<< std::endl;
          return 1;
        }
        TVector3 StartPosition_vec = {-999,-999,-999};
        ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,StartPosition_vec);
        StartPosition.SetCoordinates(StartPosition_vec.X(),StartPosition_vec.Y(),StartPosition_vec.Z());
      }
      else{
        StartPosition = InitialTrack.Start();
      }

      //Loop over the trajectory points and find two corresponding trajectory points where the angle between themselves (or themsleves and the start position) is less the fMinAngle.
      for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints()-2; ++traj){
        ++trajpoint;

        //ignore bogus info.
        auto trajflags = InitialTrack.FlagsAtPoint(trajpoint);
        if(trajflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}


        bool bail = false;

        //ignore bogus info.
        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}

        //find the next non bogus traj point.
        int nexttraj = traj+1;
        auto nextflags = InitialTrack.FlagsAtPoint(nexttraj);
        while(nextflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)){
          if(nexttraj == (int)InitialTrack.NumberTrajectoryPoints()-2){bail = true;break;}
          ++nexttraj;
          nextflags = InitialTrack.FlagsAtPoint(nexttraj);
        }

        //find the next next non bogus traj point.
        int nextnexttraj = nexttraj+1;
        auto nextnextflags = InitialTrack.FlagsAtPoint(nextnexttraj);
        while(nextnextflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)){
          if(nexttraj == (int)InitialTrack.NumberTrajectoryPoints()-1){bail = true;break;}
          ++nextnexttraj;
          nextnextflags = InitialTrack.FlagsAtPoint(nextnexttraj);
        }

        if(bail){
          mf::LogError("ShowerSmartTrackTrajectoryPointDirection") << "Trajectory point not set as rest of the traj points are bogus."<< std::endl;
          break;
        }

	//Get the directions.
        geo::Vector_t TrajPosition_vec      = InitialTrack.LocationAtPoint(traj) - StartPosition ;
        geo::Vector_t NextTrajPosition_vec;
        geo::Vector_t NextNextTrajPosition_vec;
        if(fUseStartPos){
          NextTrajPosition_vec  = InitialTrack.LocationAtPoint(nexttraj) - StartPosition;
          NextNextTrajPosition_vec  = InitialTrack.LocationAtPoint(nextnexttraj) - StartPosition;
        }
        else{
          NextTrajPosition_vec  = InitialTrack.LocationAtPoint(nexttraj)  - InitialTrack.LocationAtPoint(traj);
          NextNextTrajPosition_vec  = InitialTrack.LocationAtPoint(nextnexttraj) - InitialTrack.LocationAtPoint(traj+1);
        }

	//Get the directions.
        TVector3 TrajPosition = {TrajPosition_vec.X(), TrajPosition_vec.Y(),TrajPosition_vec.Z()};
        TVector3 NextTrajPosition = {NextTrajPosition_vec.X(), NextTrajPosition_vec.Y(),NextTrajPosition_vec.Z()};
        TVector3 NextNextTrajPosition = {NextNextTrajPosition_vec.X(), NextNextTrajPosition_vec.Y(),NextNextTrajPosition_vec.Z()};

        //Might still be bogus and we can't use the start point
        if(TrajPosition.Mag() == 0){continue;}
        if(NextTrajPosition.Mag() == 0){continue;}
        if(NextNextTrajPosition.Mag() == 0){continue;}

	//Check to see if the angle between the directions is small enough.
        if(TrajPosition.Angle(NextTrajPosition) < fAngleCut && TrajPosition.Angle(NextNextTrajPosition) <fAngleCut){break;}

	//Move the start position onwords.
        if(fAllowDynamicSliding){
          StartPosition = {InitialTrack.LocationAtPoint(traj).X(), InitialTrack.LocationAtPoint(traj).Y(), InitialTrack.LocationAtPoint(traj).Z()};
        }
      }
      
      geo::Point_t  TrajPosition   = InitialTrack.LocationAtPoint(trajpoint);
      Direction_vec  = (TrajPosition - StartPosition).Unit();
    }
    else{
      //Loop over the trajectory points and find two corresponding trajectory points where the angle between themselves (or themsleves and the start position) is less the fMinAngle.
      for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints()-2; ++traj){
        ++trajpoint;

        //ignore bogus info.
        auto trajflags = InitialTrack.FlagsAtPoint(trajpoint);
        if(trajflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}


        //ignore bogus info.
        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}

        bool bail = false;

        geo::Vector_t TrajDirection_vec;

	//Get the next non bogus trajectory points
        if(fUseStartPos){
          int prevtraj = 0;
          auto prevflags = InitialTrack.FlagsAtPoint(prevtraj);
          while(prevflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)){
            if(prevtraj == (int)InitialTrack.NumberTrajectoryPoints()-2){bail = true;break;}
            ++prevtraj;
            prevflags = InitialTrack.FlagsAtPoint(prevtraj);
          }
          TrajDirection_vec         = InitialTrack.DirectionAtPoint(prevtraj);
        }
        else if(fAllowDynamicSliding && traj!=0){
          int prevtraj = traj-1;
          auto prevflags = InitialTrack.FlagsAtPoint(prevtraj);
          while(prevflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)){
            if(prevtraj == 0){bail = true;break;}
            --prevtraj;
            prevflags = InitialTrack.FlagsAtPoint(prevtraj);
          }
          TrajDirection_vec         = InitialTrack.DirectionAtPoint(prevtraj);
        }
        else{
          TrajDirection_vec         = InitialTrack.DirectionAtPoint(traj);
        }

        //find the next non bogus traj point.
        int nexttraj = traj+1;
        auto nextflags = InitialTrack.FlagsAtPoint(nexttraj);
        while(nextflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)){
          if(nexttraj == (int)InitialTrack.NumberTrajectoryPoints()-2){bail = true;break;}
          ++nexttraj;
          nextflags = InitialTrack.FlagsAtPoint(nexttraj);
        }

        //find the next next non bogus traj point.
        int nextnexttraj = nexttraj+1;
        auto nextnextflags = InitialTrack.FlagsAtPoint(nextnexttraj);
        while(nextnextflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)){
          if(nexttraj == (int)InitialTrack.NumberTrajectoryPoints()-1){bail = true;break;}
          ++nextnexttraj;
          nextnextflags = InitialTrack.FlagsAtPoint(nextnexttraj);
        }

        if(bail){
          mf::LogError("ShowerSmartTrackTrajectoryPointDirection") 
	    << "Trajectory point not set as rest of the traj points are bogus."<< std::endl;
          break;
        }

	//Get the directions.
        geo::Vector_t NextTrajDirection_vec     = InitialTrack.DirectionAtPoint(nexttraj);
        geo::Vector_t NextNextTrajDirection_vec = InitialTrack.DirectionAtPoint(nextnexttraj);

        TVector3 TrajDirection = {TrajDirection_vec.X(), TrajDirection_vec.Y(),TrajDirection_vec.Z()};
        TVector3 NextTrajDirection = {NextTrajDirection_vec.X(), NextTrajDirection_vec.Y(),NextTrajDirection_vec.Z()};
        TVector3 NextNextTrajDirection = {NextNextTrajDirection_vec.X(), NextNextTrajDirection_vec.Y(),NextNextTrajDirection_vec.Z()};

        //Might still be bogus and we can't use the start point
        if(TrajDirection.Mag() == 0){continue;}
        if(NextTrajDirection.Mag() == 0){continue;}
        if(NextNextTrajDirection.Mag() == 0){continue;}

	//See if the angle is small enough.
        if(TrajDirection.Angle(NextTrajDirection) < fAngleCut && TrajDirection.Angle(NextNextTrajDirection) <fAngleCut){
	  break;
	}
      }
      Direction_vec  = InitialTrack.DirectionAtPoint(trajpoint).Unit();
    }

    if(trajpoint == (int) InitialTrack.NumberTrajectoryPoints() -3){
      mf::LogError("ShowerSmartTrackTrajectoryPointDirectio") << 
	"Trajectory point not set."<< std::endl;
      return 1;
    }

    //Set the direction.
    TVector3 Direction = {Direction_vec.X(), Direction_vec.Y(),Direction_vec.Z()};
    TVector3 DirectionErr = {-999,-999,-999};
    ShowerEleHolder.SetElement(Direction,DirectionErr,fShowerDirectionOuputLabel);
    return 0;
  }
}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerSmartTrackTrajectoryPointDirection)

