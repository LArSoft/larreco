////////////////////////////////////////////////////////////////////////
// Class: ShowerEnergyAlg
// File:  ShowerEnergyAlg.cxx
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ShowerEnergyAlg.h"

shower::ShowerEnergyAlg::ShowerEnergyAlg(fhicl::ParameterSet const& pset) {
  fUGradient  = pset.get<double>("UGradient");
  fUIntercept = pset.get<double>("UIntercept");
  fVGradient  = pset.get<double>("VGradient");
  fVIntercept = pset.get<double>("VIntercept");
  fZGradient  = pset.get<double>("ZGradient");
  fZIntercept = pset.get<double>("ZIntercept");
}

double shower::ShowerEnergyAlg::ShowerEnergy(art::PtrVector<recob::Hit> const& hits, int plane) {

  /// Finds the total energy deposited by the shower in this view

  double totalCharge = 0, totalEnergy = 0;

  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    totalCharge += ( (*hit)->Integral() * TMath::Exp( (detprop->SamplingRate() * (*hit)->PeakTime()) / (larprop->ElectronLifetime()*1e3) ) );

  switch (plane) {
  case 0:
    totalEnergy = (double)(totalCharge - fUIntercept)/(double)fUGradient;
    break;
  case 1:
    totalEnergy = (double)(totalCharge - fVIntercept)/(double)fVGradient;
    break;
  case 2:
    totalEnergy = (double)(totalCharge - fZIntercept)/(double)fZGradient;
    break;
  }

  return totalEnergy;

}
