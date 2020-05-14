//############################################################################
//### Name:        ShowerSlidingStandardCalodEdx                           ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk)         ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. This        ###
//###              takes the sliding fit trajectory to make a 3D dEdx.     ###
//###              This module is best used with the sliding linear fit    ###
//###              and ShowerTrackTrajToSpacepoint                         ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TRACSAlg.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools {

  class ShowerSlidingStandardCalodEdx : IShowerTool {

  public:
    ShowerSlidingStandardCalodEdx(const fhicl::ParameterSet& pset);

    ~ShowerSlidingStandardCalodEdx();

    //Physics Function. Calculate the dEdx.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //Servcies and Algorithms
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg fCalorimetryAlg;
    detinfo::DetectorProperties const* fDetProp;

    //fcl parameters
    float fMinAngleToWire;  //Minimum angle between the wire direction and the shower
                            //direction for the spacepoint to be used. Default means
                            //the cut has no effect. In radians.
    float fShapingTime;     //Shaping time of the ASIC defualt so we don't cut on track
                            //going too much into the plane. In Microseconds
    float fMinDistCutOff;   //Distance in wires a hit has to be from the start position
                            //to be used
    float fMaxDist;         //Distance in wires a that a trajectory point can be from a
                            //spacepoint to match to it.
    float fdEdxTrackLength; //Max Distance a spacepoint can be away from the start of the
                            //track. In cm
    bool fUseMedian;        //Use the median value as the dEdx rather than the mean.
    bool fCutStartPosition; //Remove hits using MinDistCutOff from the vertex as well.
    art::InputTag fPFParticleModuleLabel;

    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackSpacePointsInputLabel;
    std::string fInitialTrackInputLabel;
    std::string fShowerdEdxOuputLabel;
    std::string fShowerBestPlaneOutputLabel;
    std::string fShowerdEdxVecOuputLabel;
  };

  ShowerSlidingStandardCalodEdx::ShowerSlidingStandardCalodEdx(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
    , fMinAngleToWire(pset.get<float>("MinAngleToWire"))
    , fShapingTime(pset.get<float>("ShapingTime"))
    , fMinDistCutOff(pset.get<float>("MinDistCutOff"))
    , fMaxDist(pset.get<float>("MaxDist"))
    , fdEdxTrackLength(pset.get<float>("dEdxTrackLength"))
    , fUseMedian(pset.get<bool>("UseMedian"))
    , fCutStartPosition(pset.get<bool>("CutStartPosition"))
    , fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel"))
    , fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel"))
    , fShowerdEdxOuputLabel(pset.get<std::string>("ShowerdEdxOuputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
    , fShowerdEdxVecOuputLabel(pset.get<std::string>("ShowerdEdxVecOuputLabel"))
  {}

  ShowerSlidingStandardCalodEdx::~ShowerSlidingStandardCalodEdx() {}

  int
  ShowerSlidingStandardCalodEdx::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    // Shower dEdx calculation
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      mf::LogError("ShowerSlidingStandardCalodEdx")
        << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)) {
      mf::LogError("ShowerSlidingStandardCalodEdx")
        << "Initial Track Spacepoints is not set returning" << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)) {
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Initial Track is not set" << std::endl;
      return 1;
    }

    //Get the initial track hits
    std::vector<art::Ptr<recob::SpacePoint>> tracksps;
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel, tracksps);

    if (tracksps.size() == 0) {
      mf::LogWarning("ShowerSlidingStandardCalodEdx")
        << "no spacepointsin the initial track" << std::endl;
      return 0;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint>> spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)) {
      throw cet::exception("ShowerSlidingStandardCalodEdx")
        << "Could not configure the spacepoint handle. Something is configured incorrectly. "
           "Stopping";
      return 1;
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit> fmsp(spHandle, Event, fPFParticleModuleLabel);
    if (!fmsp.isValid()) {
      throw cet::exception("ShowerSlidingStandardCalodEdx")
        << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    //Only consider hits in the same tpcs as the vertex.
    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition);

    //Get the initial track
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement(fInitialTrackInputLabel, InitialTrack);

    //Don't care that I could use a vector.
    std::map<int, std::vector<double>> dEdx_vec;
    std::map<int, std::vector<double>> dEdx_vecErr;
    std::map<int, int> num_hits;

    for (geo::PlaneID plane_id : fGeom->IteratePlaneIDs()) {
      dEdx_vec[plane_id.Plane] = {};
      dEdx_vecErr[plane_id.Plane] = {};
      num_hits[plane_id.Plane] = 0;
    }

    //Loop over the spacepoints
    for (auto const sp : tracksps) {

      //Get the associated hit
      std::vector<art::Ptr<recob::Hit>> hits = fmsp.at(sp.key());
      if (hits.size() == 0) {
        mf::LogWarning("ShowerSlidingStandardCalodEdx")
          << "no hit for the spacepoint. This suggest the find many is wrong." << std::endl;
        continue;
      }
      const art::Ptr<recob::Hit> hit = hits[0];
      double wirepitch = fGeom->WirePitch((geo::PlaneID)hit->WireID());

      //Only consider hits in the same tpc
      geo::PlaneID planeid = hit->WireID();
      geo::TPCID TPC = planeid.asTPCID();
      if (TPC != vtxTPC) { continue; }

      //Ignore spacepoints within a few wires of the vertex.
      double dist_from_start =
        (IShowerTool::GetTRACSAlg().SpacePointPosition(sp) - ShowerStartPosition).Mag();

      if (fCutStartPosition) {
        if (dist_from_start < fMinDistCutOff * wirepitch) { continue; }

        if (dist_from_start > fdEdxTrackLength) { continue; }
      }

      //Find the closest trajectory point of the track. These should be in order if the user has used ShowerTrackTrajToSpacepoint_tool but the sake of gernicness I'll get the cloest sp.
      unsigned int index = 999;
      double MinDist = 999;
      for (unsigned int traj = 0; traj < InitialTrack.NumberTrajectoryPoints(); ++traj) {

        geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(traj);
        TVector3 TrajPosition = {
          TrajPositionPoint.X(), TrajPositionPoint.Y(), TrajPositionPoint.Z()};

        //ignore bogus info.
        auto flags = InitialTrack.FlagsAtPoint(traj);
        if (flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)) { continue; }

        TVector3 pos = IShowerTool::GetTRACSAlg().SpacePointPosition(sp) - TrajPosition;

        if (pos.Mag() < MinDist && pos.Mag() < fMaxDist * wirepitch) {
          MinDist = pos.Mag();
          index = traj;
        }
      }

      //If there is no matching trajectory point then bail.
      if (index == 999) { continue; }

      geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(index);
      TVector3 TrajPosition = {TrajPositionPoint.X(), TrajPositionPoint.Y(), TrajPositionPoint.Z()};

      geo::Point_t TrajPositionStartPoint = InitialTrack.LocationAtPoint(0);
      TVector3 TrajPositionStart = {
        TrajPositionStartPoint.X(), TrajPositionStartPoint.Y(), TrajPositionStartPoint.Z()};

      //Ignore values with 0 mag from the start position
      if ((TrajPosition - TrajPositionStart).Mag() == 0) { continue; }
      if ((TrajPosition - ShowerStartPosition).Mag() == 0) { continue; }

      if ((TrajPosition - TrajPositionStart).Mag() < fMinDistCutOff * wirepitch) { continue; }

      if ((TrajPosition - TrajPositionStart).Mag() > fdEdxTrackLength) { continue; }

      //Get the direction of the trajectory point
      geo::Vector_t TrajDirection_vec = InitialTrack.DirectionAtPoint(index);
      TVector3 TrajDirection = {
        TrajDirection_vec.X(), TrajDirection_vec.Y(), TrajDirection_vec.Z()};

      //If the direction is in the same direction as the wires within some tolerance the hit finding struggles. Let remove these.
      TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection();

      if (TrajDirection.Angle(PlaneDirection) < fMinAngleToWire) {
        mf::LogWarning("ShowerSlidingStandardCalodEdx") << "remove from angle cut" << std::endl;
        continue;
      }

      //If the direction is too much into the wire plane then the shaping amplifer cuts the charge. Lets remove these events.
      double velocity = fDetProp->DriftVelocity(fDetProp->Efield(), fDetProp->Temperature());
      double distance_in_x = TrajDirection.X() * (wirepitch / TrajDirection.Dot(PlaneDirection));
      double time_taken = TMath::Abs(distance_in_x / velocity);

      //Shaping time doesn't seem to exist in a global place so add it as a fcl.
      if (fShapingTime < time_taken) {
        mf::LogWarning("ShowerSlidingStandardCalodEdx") << "move for shaping time" << std::endl;
        continue;
      }

      //If we still exist then we can be used in the calculation. Calculate the 3D pitch
      double trackpitch = (TrajDirection * (wirepitch / TrajDirection.Dot(PlaneDirection))).Mag();

      //Calculate the dQdx
      double dQdx = hit->Integral() / trackpitch;

      //Calculate the dEdx
      double dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, hit->PeakTime(), planeid.Plane);

      //Add the value to the dEdx
      dEdx_vec[planeid.Plane].push_back(dEdx);

      //Iterate the number of hits on the plane
      ++num_hits[planeid.Plane];
    }

    //Never have the stats to do a landau fit and get the most probable value. User decides if they want the median value or the mean.
    std::vector<double> dEdx_val;
    std::vector<double> dEdx_valErr;
    for (auto const& dEdx_plane : dEdx_vec) {

      if ((dEdx_plane.second).size() == 0) {
        dEdx_val.push_back(-999);
        dEdx_valErr.push_back(-999);
        continue;
      }

      if (fUseMedian) {
        dEdx_val.push_back(TMath::Median((dEdx_plane.second).size(), &(dEdx_plane.second)[0]));
      }
      else {
        //Else calculate the mean value.
        double dEdx_mean = 0;
        for (auto const& dEdx : dEdx_plane.second) {
          if (dEdx > 10 || dEdx < 0) { continue; }
          dEdx_mean += dEdx;
        }
        dEdx_val.push_back(dEdx_mean / (float)(dEdx_plane.second).size());
      }
    }

    //Work out which is the best plane from the most hits.
    int max_hits = -999;
    int best_plane = -999;
    for (auto const& num_hits_plane : num_hits) {
      if (num_hits_plane.second > max_hits) {
        best_plane = num_hits_plane.first;
        max_hits = num_hits_plane.second;
      }
    }

    //Need to sort out errors sensibly.
    ShowerEleHolder.SetElement(dEdx_val, dEdx_valErr, fShowerdEdxOuputLabel);
    ShowerEleHolder.SetElement(best_plane, fShowerBestPlaneOutputLabel);
    ShowerEleHolder.SetElement(dEdx_vec, fShowerdEdxVecOuputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerSlidingStandardCalodEdx)
