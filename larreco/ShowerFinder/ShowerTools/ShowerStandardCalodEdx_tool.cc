//############################################################################
//### Name:        ShowerStandardCalodEdx                                  ###
//### Author:      Ed Tyley                                                ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. Derived     ###
//###              from the EMShower_module.cc                             ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools {

  class ShowerStandardCalodEdx : IShowerTool {

  public:
    ShowerStandardCalodEdx(const fhicl::ParameterSet& pset);

    ~ShowerStandardCalodEdx();

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //Define the services and algorithms
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg fCalorimetryAlg;

    //fcl parameters.
    double fdEdxTrackLength; //Max length from a hit can be to the start point in cm.
    bool fMaxHitPlane;       //Set the best planes as the one with the most hits
    bool fMissFirstPoint;    //Do not use any hits from the first wire.
    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackHitsInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fShowerdEdxOutputLabel;
    std::string fShowerBestPlaneOutputLabel;
  };

  ShowerStandardCalodEdx::ShowerStandardCalodEdx(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fdEdxTrackLength(pset.get<float>("dEdxTrackLength"))
    , fMaxHitPlane(pset.get<bool>("MaxHitPlane"))
    , fMissFirstPoint(pset.get<bool>("MissFirstPoint"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    , fShowerdEdxOutputLabel(pset.get<std::string>("ShowerdEdxOutputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
  {}

  ShowerStandardCalodEdx::~ShowerStandardCalodEdx() {}

  int
  ShowerStandardCalodEdx::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                           art::Event& Event,
                                           reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    // Shower dEdx calculation
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      mf::LogError("ShowerStandardCalodEdx") << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackHitsInputLabel)) {
      mf::LogError("ShowerStandardCalodEdx") << "Initial Track Hits not set returning" << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
      mf::LogError("ShowerStandardCalodEdx") << "Shower Direction not set" << std::endl;
      return 1;
    }

    //Get the initial track hits
    std::vector<art::Ptr<recob::Hit>> trackhits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel, trackhits);

    if (trackhits.size() == 0) {
      mf::LogWarning("ShowerStandardCalodEdx") << "Not Hits in the initial track" << std::endl;
      return 0;
    }

    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    TVector3 showerDir = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, showerDir);

    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition);

    // Split the track hits per plane
    std::vector<double> dEdxVec;
    std::vector<std::vector<art::Ptr<recob::Hit>>> trackHits;
    unsigned int numPlanes = fGeom->Nplanes();
    trackHits.resize(numPlanes);

    // Loop over the track hits and split into planes
    for (unsigned int hitIt = 0; hitIt < trackhits.size(); ++hitIt) {
      art::Ptr<recob::Hit> hit = trackhits.at(hitIt);
      geo::PlaneID hitWire = hit->WireID();
      geo::TPCID TPC = hitWire.asTPCID();

      //only get hits from the same TPC as the vertex
      if (TPC == vtxTPC) { (trackHits.at(hitWire.Plane)).push_back(hit); }
    }

    int bestHitsPlane = 0;
    int bestPlaneHits = 0;
    int bestPlane = -999;
    double minPitch = 999;

    for (unsigned int plane = 0; plane < numPlanes; ++plane) {
      std::vector<art::Ptr<recob::Hit>> trackPlaneHits = trackHits.at(plane);

      if (trackPlaneHits.size()) {

        double dEdx = -999;
        double totQ = 0;
        double avgT = 0;
        double pitch = 0;

        //Calculate the pitch
        double wirepitch = fGeom->WirePitch(trackPlaneHits.at(0)->WireID().planeID());
        double angleToVert = fGeom->WireAngleToVertical(fGeom->Plane(plane).View(),
                                                        trackPlaneHits[0]->WireID().planeID()) -
                             0.5 * TMath::Pi();
        double cosgamma =
          std::abs(sin(angleToVert) * showerDir.Y() + cos(angleToVert) * showerDir.Z());

        pitch = wirepitch / cosgamma;

        if (pitch) { // Check the pitch is calculated correctly
          int nhits = 0;
          std::vector<float> vQ;

          //Get the first wire
          int w0 = trackPlaneHits.at(0)->WireID().Wire;

          for (auto const& hit : trackPlaneHits) {

            // Get the wire for each hit
            int w1 = hit->WireID().Wire;
            if (fMissFirstPoint && w0 == w1) { continue; }

            //Ignore hits that are too far away.
            if (std::abs((w1 - w0) * pitch) < fdEdxTrackLength) {
              vQ.push_back(hit->Integral());
              totQ += hit->Integral();
              avgT += hit->PeakTime();
              ++nhits;
            }
          }

          if (totQ) {
            // Check if the pitch is the smallest yet (best plane)
            if (pitch < minPitch) {
              minPitch = pitch;
              bestPlane = plane;
            }

            //Get the median and calculate the dEdx using the algorithm.
            double dQdx = TMath::Median(vQ.size(), &vQ[0]) / pitch;
            dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avgT / nhits, trackPlaneHits[0]->WireID().Plane);

            if (isinf(dEdx)) { dEdx = -999; };

            if (nhits > bestPlaneHits || ((nhits == bestPlaneHits) && (pitch < minPitch))) {
              bestHitsPlane = plane;
              bestPlaneHits = nhits;
            }
          }
          dEdxVec.push_back(dEdx);
        }
        else {
          throw cet::exception("ShowerStandardCalodEdx")
            << "pitch is 0. I can't think how it is 0? Stopping so I can tell you" << std::endl;
        }
      }
      else { // if not (trackPlaneHits.size())
        dEdxVec.push_back(-999);
      }
      trackPlaneHits.clear();
    } //end loop over planes

    //TODO
    std::vector<double> dEdxVecErr = {-999, -999, -999};

    ShowerEleHolder.SetElement(dEdxVec, dEdxVecErr, fShowerdEdxOutputLabel);

    //Set The best plane
    if (fMaxHitPlane) { bestPlane = bestHitsPlane; }

    if (bestPlane == -999) {
      throw cet::exception("ShowerStandardCalodEdx") << "No best plane set";
      return 1;
    }
    else {
      ShowerEleHolder.SetElement(bestPlane, fShowerBestPlaneOutputLabel);
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerStandardCalodEdx)
