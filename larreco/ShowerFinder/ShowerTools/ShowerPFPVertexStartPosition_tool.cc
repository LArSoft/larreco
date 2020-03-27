//############################################################################
//### Name:        ShowerPFPVertexStartPosition                            ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the start poistion                     ###
//###              methods.                                                ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/TRACSAlg.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes
#include "TMath.h"
#include "TVector3.h"

namespace ShowerRecoTools {

  class ShowerPFPVertexStartPosition : public IShowerTool {

  public:
    ShowerPFPVertexStartPosition(const fhicl::ParameterSet& pset);

    ~ShowerPFPVertexStartPosition();

    //Calculate the start position
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //fcl parameters
    art::InputTag fPFParticleModuleLabel;
    std::string fShowerStartPositionOutputLabel;
    std::string fShowerDirectionInputLabel;
  };

  ShowerPFPVertexStartPosition::ShowerPFPVertexStartPosition(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel"))
    , fShowerStartPositionOutputLabel(pset.get<std::string>("ShowerStartPositionOutputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  {}

  ShowerPFPVertexStartPosition::~ShowerPFPVertexStartPosition() {}

  int
  ShowerPFPVertexStartPosition::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                                 art::Event& Event,
                                                 reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //Get the vertices.
    art::Handle<std::vector<recob::Vertex>> vtxHandle;
    std::vector<art::Ptr<recob::Vertex>> vertices;
    if (Event.getByLabel(fPFParticleModuleLabel, vtxHandle))
      art::fill_ptr_vector(vertices, vtxHandle);
    else {
      throw cet::exception("ShowerPFPVertexStartPosition")
        << "Could not get the pandora vertices. Something is not configured correctly. Please give "
           "the correct pandora module label. Stopping";
      return 1;
    }

    //Get the spacepoints handle and the hit assoication
    art::Handle<std::vector<recob::SpacePoint>> spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)) {
      throw cet::exception("ShowerPFPVertexStartPosition")
        << "Could not configure the spacepoint handle. Something is configured incorrectly. "
           "Stopping";
      return 1;
    }
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if (!fmh.isValid()) {
      throw cet::exception("ShowerPFPVertexStartPosition")
        << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)) {
      throw cet::exception("ShowerPFPVertexStartPosition")
        << "Could not get the pandora pf particles. Something is not cofingured coreectly Please "
           "give the correct pandoa module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Vertex> fmv(pfpHandle, Event, fPFParticleModuleLabel);
    if (!fmv.isValid()) {
      throw cet::exception("ShowerPFPVertexStartPosition")
        << "Vertex and PF particle association is somehow not valid. Stopping";
      return 1;
    }
    std::vector<art::Ptr<recob::Vertex>> vtx_cand;
    try {
      vtx_cand = fmv.at(pfparticle.key());
    }
    catch (...) {
      mf::LogError("ShowerPFPVertexStartPosition") << "PFP-Vertex assan not set, returning";
      return 1;
    }
    //If there is more than one then fail becuase I don't think that this can be the case
    if (vtx_cand.size() != 1) {
      mf::LogError("ShowerPFPVertexStartPosition")
        << "Wrong number of vertices: " << vtx_cand.size() << ", returning";
      return 1;
    }

    //If there is only one vertex good news we just say that is the start of the shower.
    if (vtx_cand.size() == 1) {
      art::Ptr<recob::Vertex> StartPositionVertex = vtx_cand[0];
      double xyz[3] = {-999, -999, -999};
      StartPositionVertex->XYZ(xyz);
      TVector3 ShowerStartPosition = {xyz[0], xyz[1], xyz[2]};
      TVector3 ShowerStartPositionErr = {-999, -999, -999};
      ShowerEleHolder.SetElement(
        ShowerStartPosition, ShowerStartPositionErr, fShowerStartPositionOutputLabel);
      return 0;
    }

    //If we there have none then use the direction to find the neutrino vertex
    if (ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {

      TVector3 ShowerDirection = {-999, -999, -999};
      ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

      art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);

      if (!fmspp.isValid()) {
        throw cet::exception("ShowerPFPVertexStartPosition")
          << "Trying to get the spacepoints and failed. Something is not configured correctly. "
             "Stopping ";
        return 1;
      }

      //Get the spacepoints
      std::vector<art::Ptr<recob::SpacePoint>> spacePoints_pfp = fmspp.at(pfparticle.key());

      //Cannot continue if we have no spacepoints
      if (spacePoints_pfp.size() == 0) { return 0; }

      // Get the Shower Center
      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);
      TVector3 ShowerCentre =
        IShowerTool::GetTRACSAlg().ShowerCentre(clockData, detProp, spacePoints_pfp, fmh);

      //Order the Hits from the shower centre. The most negative will be the start position.
      IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(
        spacePoints_pfp, ShowerCentre, ShowerDirection);

      //Set the start position.
      TVector3 ShowerStartPosition =
        IShowerTool::GetTRACSAlg().SpacePointPosition(spacePoints_pfp[0]);

      TVector3 ShowerStartPositionErr = {-999, -999, -999};
      ShowerEleHolder.SetElement(
        ShowerStartPosition, ShowerStartPositionErr, fShowerStartPositionOutputLabel);

      return 0;
    }

    mf::LogWarning("ShowerPFPVertexStartPosition")
      << "Start Position has not been set yet. If you are not calculating the start position again "
         "then maybe you should stop";
    return 0;
  }

}
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPFPVertexStartPosition)
