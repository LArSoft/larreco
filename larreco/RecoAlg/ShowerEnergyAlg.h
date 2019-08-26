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
#include "canvas/Persistency/Common/Ptr.h"
namespace fhicl { class ParameterSet; }

// larsoft
#include "lardataobj/RecoBase/Hit.h"
namespace detinfo { class DetectorProperties; }

namespace shower {
  class ShowerEnergyAlg;
}

class shower::ShowerEnergyAlg {
 public:

  ShowerEnergyAlg(fhicl::ParameterSet const& pset);

  /// Finds the total energy deposited by the shower in this view
  double ShowerEnergy(std::vector<art::Ptr<recob::Hit> > const& hits, int plane) const;

 private:

  double fUGradient, fUIntercept;
  double fVGradient, fVIntercept;
  double fZGradient, fZIntercept;

  detinfo::DetectorProperties const* detprop = nullptr;

};

#endif
