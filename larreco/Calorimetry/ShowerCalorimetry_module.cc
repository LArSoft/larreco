////////////////////////////////////////////////////////////////////////
// Class:       ShowerCalorimetry
// Plugin Type: producer (art v3_02_06)
// File:        ShowerCalorimetry_module.cc
// Authors:     calcuttj@msu.edu
//              ahiguera@Central.UH.EDU
//              tjyang@fnal.gov
// Generated at Fri Jul 12 14:14:46 2019 by Jacob Calcutt using cetskelgen
// from cetlib version v3_07_02.
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include <TVector3.h>
#include <memory>

namespace calo {
  class ShowerCalorimetry;
}

class calo::ShowerCalorimetry : public art::EDProducer {
public:
  explicit ShowerCalorimetry(fhicl::ParameterSet const& p);

private:
  void produce(art::Event& e) override;

  art::InputTag fShowerTag;
  art::InputTag fSpacePointTag;
  bool fSCE;
  CalorimetryAlg caloAlg;
};

calo::ShowerCalorimetry::ShowerCalorimetry(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fShowerTag(p.get<art::InputTag>("ShowerTag"))
  , fSpacePointTag(p.get<art::InputTag>("SpacePointTag"))
  , fSCE(p.get<bool>("CorrectSCE"))
  , caloAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  produces<std::vector<anab::Calorimetry>>();
  produces<art::Assns<recob::Shower, anab::Calorimetry>>();
}

void
calo::ShowerCalorimetry::produce(art::Event& e)
{

  art::ServiceHandle<geo::Geometry> geom;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  //Make the container for the calo product to put onto the event.
  auto caloPtr = std::make_unique<std::vector<anab::Calorimetry>>();
  auto& caloVector = *caloPtr;

  //Make a container for the track<-->calo associations.
  //One entry per track, with entry equal to index in calorimetry collection of associated object.
  std::vector<size_t> assnShowerCaloVector;
  auto associationPtr = std::make_unique<art::Assns<recob::Shower, anab::Calorimetry>>();

  auto showerHandle = e.getValidHandle<std::vector<recob::Shower>>(fShowerTag);

  //Turn it into a vector of art pointers
  std::vector<art::Ptr<recob::Shower>> recoShowers;
  art::fill_ptr_vector(recoShowers, showerHandle);

  //Also get the hits from all the showers
  art::FindManyP<recob::Hit> findHitsFromShowers(showerHandle, e, fShowerTag);

  //Go through all of the reconstructed showers in the event
  for (size_t i = 0; i < recoShowers.size(); ++i) {
    auto& shower = recoShowers[i];

    int shower_index = shower.key();
    MF_LOG_INFO("ShowerCalorimetry") << "Getting Calorimetry info for " << shower_index << "\n";

    //This wil be used in the calorimetry object later
    float shower_length = shower->Length();
    //Get the hits from this shower
    auto const& hits = findHitsFromShowers.at(shower_index);

    art::FindManyP<recob::SpacePoint> spFromShowerHits(hits, e, fSpacePointTag);

    //Sort the hits by their plane
    //This vector stores the index of each hit on each plane
    std::vector<std::vector<size_t>> hit_indices_per_plane(geom->Nplanes());
    for (size_t j = 0; j < hits.size(); ++j) {
      hit_indices_per_plane[hits[j]->WireID().Plane].push_back(j);
    }

    //Go through each plane and make calorimetry objects
    for (size_t j = 0; j < geom->Nplanes(); ++j) {

      size_t hits_in_plane = hit_indices_per_plane[j].size();

      //Reserve vectors for each part of the calorimetry object
      std::vector<float> dEdx(hits_in_plane);
      std::vector<float> dQdx(hits_in_plane);
      std::vector<float> pitch(hits_in_plane);

      //residual range, xyz, and deadwire default for now
      std::vector<float> resRange(hits_in_plane, 0.);
      std::vector<TVector3> xyz(hits_in_plane, TVector3(0., 0., 0.));
      std::vector<float> deadwires(hits_in_plane, 0.);
      std::vector<size_t> hitIndex(hits_in_plane);

      geo::PlaneID planeID;

      float kineticEnergy = 0.;

      for (size_t k = 0; k < hits_in_plane; ++k) {

        size_t hit_index = hit_indices_per_plane[j][k];
        auto& theHit = hits[hit_index];
        if (!planeID.isValid) { planeID = theHit->WireID(); }
        hitIndex[k] = theHit.key();
        float wire_pitch = geom->WirePitch(theHit->View());

        float theHit_Xpos = -999.;
        float theHit_Ypos = -999.;
        float theHit_Zpos = -999.;

        auto const& sp = spFromShowerHits.at(hit_index);
        if (!sp.empty()) {
          //only use first space point
          theHit_Xpos = sp[0]->XYZ()[0];
          theHit_Ypos = sp[0]->XYZ()[1];
          theHit_Zpos = sp[0]->XYZ()[2];
        }
        else {
          MF_LOG_INFO("ShowerCalorimetry")
            << "no sp associated w/this hit ... we will skip this hit";
          continue;
        }

        TVector3 pos(theHit_Xpos, theHit_Ypos, theHit_Zpos);
        //correct pitch for hit direction
        float this_pitch = wire_pitch;
        float angleToVert = geom->WireAngleToVertical(theHit->View(), theHit->WireID());
        angleToVert -= 0.5 * ::util::pi<>();

        float cosgamma = std::abs(sin(angleToVert) * shower->Direction().Y() +
                                  cos(angleToVert) * shower->Direction().Z());
        if (cosgamma > 0) this_pitch /= cosgamma;
        if (this_pitch < wire_pitch) this_pitch = wire_pitch;
        pitch[k] = this_pitch;

        //Correct for SCE
        if (fSCE && sce->EnableCalSpatialSCE()) {

          geo::Vector_t posOffsets = {0., 0., 0.};
          geo::Vector_t dirOffsets = {0., 0., 0.};

          posOffsets = sce->GetCalPosOffsets(geo::Point_t(pos), theHit->WireID().TPC);

          //For now, use the shower direction from Pandora...a better idea?
          dirOffsets =
            sce->GetCalPosOffsets(geo::Point_t{pos.X() + this_pitch * shower->Direction().X(),
                                               pos.Y() + this_pitch * shower->Direction().Y(),
                                               pos.Z() + this_pitch * shower->Direction().Z()},
                                  theHit->WireID().TPC);

          TVector3 dir_corr = {
            this_pitch * shower->Direction().X() - dirOffsets.X() + posOffsets.X(),
            this_pitch * shower->Direction().Y() + dirOffsets.Y() - posOffsets.Y(),
            this_pitch * shower->Direction().Z() + dirOffsets.Z() - posOffsets.Z()};

          pitch[k] = dir_corr.Mag();
          //correct position for SCE
          theHit_Xpos -= posOffsets.X();
          theHit_Ypos += posOffsets.Y();
          theHit_Zpos += posOffsets.Z();
        }
        xyz[k].SetXYZ(theHit_Xpos, theHit_Ypos, theHit_Zpos);
        resRange[k] = std::hypot(theHit_Xpos - shower->ShowerStart().X(),
                                 theHit_Ypos - shower->ShowerStart().Y(),
                                 theHit_Zpos - shower->ShowerStart().Z());

        dQdx[k] = theHit->Integral() / pitch[k];

        // Just for now, use dQdx for dEdx
        // dEdx[k] = theHit->Integral() / this_pitch;
        dEdx[k] = caloAlg.dEdx_AREA(clockData, detProp, *theHit, pitch[k]),

        kineticEnergy += dEdx[k];
      }

      //Make a calo object in the vector
      caloVector.emplace_back(kineticEnergy,
                              dEdx,
                              dQdx,
                              resRange,
                              deadwires,
                              shower_length,
                              pitch,
                              recob::tracking::convertCollToPoint(xyz),
                              hitIndex,
                              planeID);

      //Place the shower index in the association object
      assnShowerCaloVector.emplace_back(shower_index);
    }
  }

  //Make the associations for ART
  for (size_t i = 0; i < assnShowerCaloVector.size(); i++) {
    if (assnShowerCaloVector[i] == std::numeric_limits<size_t>::max()) continue;

    art::Ptr<recob::Shower> shower_ptr(showerHandle, assnShowerCaloVector[i]);
    util::CreateAssn(e, caloVector, shower_ptr, *associationPtr, i);
  }

  //Finish up: Put the objects into the event
  e.put(std::move(caloPtr));
  e.put(std::move(associationPtr));
}

DEFINE_ART_MODULE(calo::ShowerCalorimetry)
