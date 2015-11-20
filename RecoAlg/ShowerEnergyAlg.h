////////////////////////////////////////////////////////////////////////
// Class: ShowerEnergyAlg
// File:  ShowerEnergyAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#ifndef ShowerEnergyAlg_hxx
#define ShowerEnergyAlg_hxx

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// larsoft
#include "DetectorInfo/DetectorProperties.h"
#include "RecoBase/Hit.h"

// ROOT
#include "TMath.h"

namespace shower {
  class ShowerEnergyAlg;
}

class shower::ShowerEnergyAlg {
 public:

  ShowerEnergyAlg(fhicl::ParameterSet const& pset);
  double ShowerEnergy(art::PtrVector<recob::Hit> const& hits, int plane);

 private:

  double fUGradient, fUIntercept;
  double fVGradient, fVIntercept;
  double fZGradient, fZIntercept;

  detinfo::DetectorProperties const* detprop = nullptr;

};

#endif
