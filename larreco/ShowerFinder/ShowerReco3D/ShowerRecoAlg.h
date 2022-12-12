#ifndef RECOTOOL_SHOWERRECOALG_H
#define RECOTOOL_SHOWERRECOALG_H

#include <vector>

#include "ShowerRecoAlgBase.h"
#include "larcorealg/Geometry/fwd.h"
#include "lardataobj/RecoBase/Shower.h"

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace showerreco {

  class ShowerRecoAlg : public ShowerRecoAlgBase {
  public:
    /// Function to decide if to use Area or Pulse Amplitude for calculations
    void SetUseArea(bool on) { fUseArea = on; }

    void Verbose(bool verbose) { fVerbosity = verbose; }

    void CaloAlgo(calo::CalorimetryAlg* alg) { fCaloAlg = alg; }

    /// Function to set whether to use E correction
    void setEcorrection(bool on) { _Ecorrection = on; }

    /// Function to reconstruct a shower
    recob::Shower RecoOneShower(geo::GeometryCore const& geom,
                                geo::ChannelMapAlg const& channelMapAlg,
                                detinfo::DetectorClocksData const& clockData,
                                detinfo::DetectorPropertiesData const& detProp,
                                std::vector<showerreco::ShowerCluster_t> const&);

  private:
    calo::CalorimetryAlg* fCaloAlg;
    bool _Ecorrection{true};
    bool fVerbosity{true};
    double fcalodEdxlength{1000};
    double fdEdxlength{2.4};
    bool fUseArea{true};
  };
}

#endif
/** @} */ // end of doxygen group
