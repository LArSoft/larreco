////////////////////////////////////////////////////////////////////////
// Class: ShowerEnergyAlg
// File:  ShowerEnergyAlg.cxx
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Common/PtrVector.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/ArtDataHelper/ToElement.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"
#include "range/v3/numeric.hpp"
#include "range/v3/view.hpp"

#include <cmath>

using lar::to_element;

shower::ShowerEnergyAlg::ShowerEnergyAlg(fhicl::ParameterSet const& pset)
  : fLinearFunctions{{{pset.get<double>("UGradient"), pset.get<double>("UIntercept")},
                      {pset.get<double>("VGradient"), pset.get<double>("VIntercept")},
                      {pset.get<double>("ZGradient"), pset.get<double>("ZIntercept")}}}
{}

double
shower::ShowerEnergyAlg::ShowerEnergy(detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detprop,
                                      std::vector<art::Ptr<recob::Hit>> const& hits,
                                      geo::PlaneID::PlaneID_t const plane) const
{
  // Should we throw instead if the plane is not in the range [0,3)?
  if (plane >= fLinearFunctions.size()) { return 0.; }

  auto const coeff = sampling_rate(clockData) / (detprop.ElectronLifetime() * 1e3);

  auto in_plane = [plane](auto const& hit) { return hit.WireID().Plane == plane; };
  auto charge = [coeff](auto const& hit) {
    return hit.Integral() * std::exp(coeff * hit.PeakTime());
  };

  double const totalCharge =
    ranges::accumulate(hits | ranges::views::transform(to_element) | ranges::views::filter(in_plane) |
                         ranges::views::transform(charge),
                       0.);

  return fLinearFunctions[plane].energy_from(totalCharge);
}
