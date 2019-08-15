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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"

#include "TMath.h"

shower::ShowerEnergyAlg::ShowerEnergyAlg(fhicl::ParameterSet const& pset)
  : detprop(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
  fUGradient  = pset.get<double>("UGradient");
  fUIntercept = pset.get<double>("UIntercept");
  fVGradient  = pset.get<double>("VGradient");
  fVIntercept = pset.get<double>("VIntercept");
  fZGradient  = pset.get<double>("ZGradient");
  fZIntercept = pset.get<double>("ZIntercept");
}

double shower::ShowerEnergyAlg::ShowerEnergy(std::vector<art::Ptr<recob::Hit> > const& hits, int plane) {

  double totalCharge = 0, totalEnergy = 0;

  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit){
    if (int((*hit)->WireID().Plane)!=plane) continue;
    totalCharge += ( (*hit)->Integral() * TMath::Exp( (detprop->SamplingRate() * (*hit)->PeakTime()) / (detprop->ElectronLifetime()*1e3) ) );
  }

  switch (plane) {
  case 0:
    totalEnergy = (totalCharge * fUGradient) + fUIntercept;
    break;
  case 1:
    totalEnergy = (totalCharge * fVGradient) + fVIntercept;
    break;
  case 2:
    totalEnergy = (totalCharge * fZGradient) + fZIntercept;
    break;
  }

  return totalEnergy;

}
