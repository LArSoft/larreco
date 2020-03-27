//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes

namespace ShowerRecoTools {

  class ShowerLinearEnergy : IShowerTool {
  public:
    explicit ShowerLinearEnergy(const fhicl::ParameterSet& pset);

    //Physics Function. Calculate the shower Energy.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerElementHolder) override;

  private:
    double CalculateEnergy(detinfo::DetectorClocksData const& clockData,
                           detinfo::DetectorPropertiesData const& detProp,
                           std::vector<art::Ptr<recob::Hit>> const& hits,
                           geo::View_t view);

    // fcl parameters
    double fUGradient;  // Gradient of the linear fit of total charge to total
                        // energy on the U plane.
    double fUIntercept; // Intercept of the linear fit of total charge to total
                        // energy on the U plane.
    double fVGradient;
    double fVIntercept;
    double fZGradient;
    double fZIntercept;
    double fXGradient;
    double fXIntercept;
    double fYGradient;
    double fYIntercept;
    double f3DGradient;
    double f3DIntercept;

    art::InputTag fPFParticleModuleLabel;

    std::string fShowerEnergyOutputLabel;

    //Services
    art::ServiceHandle<geo::Geometry> fGeom;
  };

  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fUGradient(pset.get<double>("UGradient"))
    , fUIntercept(pset.get<double>("UIntercept"))
    , fVGradient(pset.get<double>("VGradient"))
    , fVIntercept(pset.get<double>("VIntercept"))
    , fZGradient(pset.get<double>("ZGradient"))
    , fZIntercept(pset.get<double>("ZIntercept"))
    , fXGradient(pset.get<double>("XGradient"))
    , fXIntercept(pset.get<double>("XIntercept"))
    , fYGradient(pset.get<double>("YGradient"))
    , fYIntercept(pset.get<double>("YIntercept"))
    , f3DGradient(pset.get<double>("ThreeDGradient"))
    , f3DIntercept(pset.get<double>("ThreeDIntercept"))
    , fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel", ""))
    , fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel"))
  {}

  int
  ShowerLinearEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                       art::Event& Event,
                                       reco::shower::ShowerElementHolder& ShowerEleHolder)
  {
    //Holder for the final product
    std::vector<double> ShowerLinearEnergy;
    unsigned int numPlanes = fGeom->Nplanes();

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)) {
      throw cet::exception("ShowerLinearEnergy")
        << "Could not get the pandora pf particles. Something is not cofingured coreectly Please "
           "give the correct pandoa module label. Stopping";
    }

    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit>>> view_hits;

    //Get the clusters
    art::Handle<std::vector<recob::Cluster>> clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)) {
      throw cet::exception("ShowerLinearEnergy")
        << "Could not get the pandora clusters. Something is not cofingured coreectly Please give "
           "the correct pandoa module label. Stopping";
    }
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster>> clusters = fmc.at(pfparticle.key());

    //Get the hit association
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    //Loop over the clusters in the plane and get the hits
    for (auto const& cluster : clusters) {

      //Get the hits
      std::vector<art::Ptr<recob::Hit>> hits = fmhc.at(cluster.key());

      //Get the view.
      geo::View_t view = cluster->View();

      view_hits[view].insert(view_hits[view].end(), hits.begin(), hits.end());
    }

    std::map<unsigned int, double> view_energies;
    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    // Accounting for events crossing the cathode.
    for (auto const& [view, hits] : view_hits) {
      unsigned int viewNum = view;
      view_energies[viewNum] = CalculateEnergy(clockData, detProp, hits, view);
    }

    //TODO think of a better way of doing this
    for (unsigned int plane = 0; plane < numPlanes; ++plane) {
      double Energy;
      try {
        Energy = view_energies.at(plane);
        if (Energy < 0) {
          mf::LogWarning("ShowerLinearEnergy") << "Negative shower energy: " << Energy;
          Energy = -999;
        }
      }
      catch (...) {
        mf::LogWarning("ShowerLinearEnergy")
          << "No energy calculation for plane " << plane << std::endl;
        Energy = -999;
      }
      ShowerLinearEnergy.push_back(Energy);
    }

    if (ShowerLinearEnergy.size() == 0) {
      throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
    }

    //TODO
    std::vector<double> EnergyError = {-999, -999, -999};

    ShowerEleHolder.SetElement(ShowerLinearEnergy, EnergyError, fShowerEnergyOutputLabel);

    return 0;
  }

  //Function to calculate the energy of a shower in a plane. Using a linear map between charge and Energy.
  //Exactly the same method as the ShowerEnergyAlg.cxx. Thanks Mike.
  double
  ShowerLinearEnergy::CalculateEnergy(detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      std::vector<art::Ptr<recob::Hit>> const& hits,
                                      geo::View_t const view)
  {
    double totalCharge = 0, totalEnergy = 0;

    for (auto const& hit : hits) {
      totalCharge += (hit->Integral() * std::exp((sampling_rate(clockData) * hit->PeakTime()) /
                                                 (detProp.ElectronLifetime() * 1e3)));
    }

    switch (view) {
    case geo::kU: totalEnergy = (totalCharge * fUGradient) + fUIntercept; break;
    case geo::kV: totalEnergy = (totalCharge * fVGradient) + fVIntercept; break;
    case geo::kW: //same as geo::kZ
      totalEnergy = (totalCharge * fZGradient) + fZIntercept;
      break;
    case geo::kX: totalEnergy = (totalCharge * fXGradient) + fXIntercept; break;
    case geo::kY: totalEnergy = (totalCharge * fYGradient) + fYIntercept; break;
    case geo::k3D: totalEnergy = (totalCharge * f3DGradient) + f3DIntercept; break;
    default: throw cet::exception("ShowerLinearEnergy") << "View Not configured";
    }

    return totalEnergy;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)
