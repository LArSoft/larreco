#include <string>
#include <optional>
#include <cmath>
#include <limits> // std::numeric_limits<>

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// ROOT includes
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

namespace calo {

class CaloChecker: public art::EDAnalyzer {
  public:
    explicit CaloChecker(fhicl::ParameterSet const& config);
    void analyze(const art::Event& evt) override;

  private:
    std::vector<std::string> fCaloLabels;
    std::string fTrackLabel;
}; 

} // end namespace calo


calo::CaloChecker::CaloChecker(fhicl::ParameterSet const& config):
  EDAnalyzer{config},
  fCaloLabels(config.get<std::vector<std::string>>("CaloLabels")),
  fTrackLabel(config.get<std::string>("TrackLabel"))
{
  assert(fCaloLabels.size() >= 2);
}

void calo::CaloChecker::analyze(const art::Event &evt) {
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackLabel, trackListHandle)) {
    art::fill_ptr_vector(tracklist, trackListHandle);
  }

  std::vector<art::FindManyP<anab::Calorimetry>> calos;
  for (const std::string &l: fCaloLabels) {
    calos.emplace_back(tracklist, evt, l);
  }

  float EPS = 1e-3;

  for (unsigned trk_i = 0; trk_i < tracklist.size(); trk_i++) {
    std::vector<art::Ptr<anab::Calorimetry>> base = calos[0].at(trk_i);

    std::cout << "New Track!\n";
    std::cout << "Base calo (" << fCaloLabels[0] << ") n calo: " << base.size() << std::endl;
    for (unsigned plane = 0; plane < base.size(); plane++) {
      std::cout << "N points on plane (" << plane << ") ID (" << base[plane]->PlaneID() << ") n points: " << base[plane]->dEdx().size() << std::endl;
    }

    for (unsigned i = 1; i < fCaloLabels.size(); i++) {
      bool equal = true;
      std::vector<art::Ptr<anab::Calorimetry>> othr = calos[i].at(trk_i);
      if (base.size() != othr.size()) {
        equal = false;
        std::cout << "Track: " << trk_i << " calo (" << fCaloLabels[0] << ") has (" << base.size() << "). Calo (" << fCaloLabels[i] << ") has size (" << othr.size() << ") mismatch." << std::endl;
      }

      for (unsigned plane = 0; plane < base.size(); plane++) {
        // check plane index
        if (base[plane]->PlaneID() != othr[plane]->PlaneID()) {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": ";
          std::cout << "Plane ID mismatch (" << base[plane]->PlaneID() << ") (" << othr[plane]->PlaneID() << ")" << std::endl;
        }

        // check the range
        if (abs(base[plane]->Range() - othr[plane]->Range()) > EPS) {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": ";
          std::cout << "Range mismatch (" << base[plane]->Range() << ") (" << othr[plane]->Range() << ")" << std::endl;
        }

        // check the kinetic energy
        if (abs(base[plane]->KineticEnergy() - othr[plane]->KineticEnergy()) > EPS) {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": ";
          std::cout << "KineticEnergy mismatch (" << base[plane]->KineticEnergy() << ") (" << othr[plane]->KineticEnergy() << ")" << std::endl;
        }

        // check dE dx
        const std::vector<float> &base_dedx = base[plane]->dEdx();        
        const std::vector<float> &othr_dedx = othr[plane]->dEdx();        

        if (base_dedx.size() == othr_dedx.size()) {
          for (unsigned i_dedx = 0; i_dedx < base_dedx.size(); i_dedx++) {
            if (abs(base_dedx[i_dedx] - othr_dedx[i_dedx]) > EPS) {
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "dEdx mismatch at index: " << i_dedx << " (" << base_dedx[i_dedx] << ") (" << othr_dedx[i_dedx] << ")" << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "dEdx size mismatch (" << base_dedx.size() << ") (" << othr_dedx.size() << ")" << std::endl;
        }

        // check dQdx
        const std::vector<float> &base_dqdx = base[plane]->dQdx();        
        const std::vector<float> &othr_dqdx = othr[plane]->dQdx();        

        if (base_dqdx.size() == othr_dqdx.size()) {
          for (unsigned i_dqdx = 0; i_dqdx < base_dqdx.size(); i_dqdx++) {
            if (abs(base_dqdx[i_dqdx] - othr_dqdx[i_dqdx]) > EPS) {
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "dQdx mismatch at index: " << i_dqdx << " (" << base_dqdx[i_dqdx] << ") (" << othr_dqdx[i_dqdx] << ")" << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "dQdx size mismatch (" << base_dqdx.size() << ") (" << othr_dqdx.size() << ")" << std::endl;
        }

        // check Residual Range
        const std::vector<float> &base_rr = base[plane]->ResidualRange();        
        const std::vector<float> &othr_rr = othr[plane]->ResidualRange();        

        if (base_rr.size() == othr_rr.size()) {
          for (unsigned i_rr = 0; i_rr < base_rr.size(); i_rr++) {
            if (abs(base_rr[i_rr] - othr_rr[i_rr]) > EPS) {
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "ResidualRange mismatch at index: " << i_rr << " (" << base_rr[i_rr] << ") (" << othr_rr[i_rr] << ")" << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "ResidualRange size mismatch (" << base_rr.size() << ") (" << othr_rr.size() << ")" << std::endl;
        }

        // check DeadWireRC
        const std::vector<float> &base_dwrr = base[plane]->DeadWireResRC();        
        const std::vector<float> &othr_dwrr = othr[plane]->DeadWireResRC();        

        if (base_dwrr.size() == othr_dwrr.size()) {
          for (unsigned i_dwrr = 0; i_dwrr < base_dwrr.size(); i_dwrr++) {
            if (abs(base_dwrr[i_dwrr] - othr_dwrr[i_dwrr]) > EPS) {
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "DeadWireResRC mismatch at index: " << i_dwrr << " (" << base_dwrr[i_dwrr] << ") (" << othr_dwrr[i_dwrr] << ")" << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "DeadWireResRC size mismatch (" << base_dwrr.size() << ") (" << othr_dwrr.size() << ")" << std::endl;
        }

        // check Track Pitch
        const std::vector<float> &base_tpv = base[plane]->TrkPitchVec();        
        const std::vector<float> &othr_tpv = othr[plane]->TrkPitchVec();        

        if (base_tpv.size() == othr_tpv.size()) {
          for (unsigned i_tpv = 0; i_tpv < base_tpv.size(); i_tpv++) {
            if (abs(base_tpv[i_tpv] - othr_tpv[i_tpv]) > EPS) {
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "TrkPitchVec mismatch at index: " << i_tpv << " (" << base_tpv[i_tpv] << ") (" << othr_tpv[i_tpv] << ")" << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "TrkPitchVec size mismatch (" << base_tpv.size() << ") (" << othr_tpv.size() << ")" << std::endl;
        }

        // check TP Indices
        const std::vector<size_t> &base_tpi = base[plane]->TpIndices();        
        const std::vector<size_t> &othr_tpi = othr[plane]->TpIndices();        

        if (base_tpi.size() == othr_tpi.size()) {
          for (unsigned i_tpi = 0; i_tpi < base_tpi.size(); i_tpi++) {
            if (base_tpi[i_tpi] != othr_tpi[i_tpi]) {
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "TpIndices mismatch at index: " << i_tpi << " (" << base_tpi[i_tpi] << ") (" << othr_tpi[i_tpi] << ")" << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "TpIndices size mismatch (" << base_tpi.size() << ") (" << othr_tpi.size() << ")" << std::endl;
        }

        // check XYZ
        const std::vector<geo::Point_t> &base_xyz = base[plane]->XYZ();        
        const std::vector<geo::Point_t> &othr_xyz = othr[plane]->XYZ();        

        if (base_xyz.size() == othr_xyz.size()) {
          for (unsigned i_xyz = 0; i_xyz < base_xyz.size(); i_xyz++) {
            if (abs(base_xyz[i_xyz].X() - othr_xyz[i_xyz].X()) > EPS || 
                abs(base_xyz[i_xyz].Y() - othr_xyz[i_xyz].Y() > EPS) || 
                abs(base_xyz[i_xyz].Z() - othr_xyz[i_xyz].Z()) > EPS) { 
              equal = false;
              std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": "; 
              std::cout << "XYZ mismatch at index: " << i_xyz; 
              std::cout << "X (" << base_xyz[i_xyz].X() << ") (" << othr_xyz[i_xyz].X() << ") ";
              std::cout << "Y (" << base_xyz[i_xyz].Y() << ") (" << othr_xyz[i_xyz].Y() << ") ";
              std::cout << "Z (" << base_xyz[i_xyz].Z() << ") (" << othr_xyz[i_xyz].Z() << ") " << std::endl;
            }
          }
        }
        else {
          equal = false;
          std::cout << "Track: " << trk_i << " calos (" << fCaloLabels[0] << ") (" << fCaloLabels[i] << ") plane " << plane << ": " << "XYZ size mismatch (" << base_xyz.size() << ") (" << othr_xyz.size() << ")" << std::endl;
        }

      }

      if (equal) {
        std::cout << "Track: " << trk_i << " calo (" << fCaloLabels[0] << ") is equal to calo (" << fCaloLabels[i] << ")" << std::endl;
      }

    }
  }

  
}

DEFINE_ART_MODULE(calo::CaloChecker)
