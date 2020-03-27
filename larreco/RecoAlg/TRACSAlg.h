#ifndef TRACSAlg_hxx
#define TRACSAlg_hxx

//Framework Includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/ShowerElementHolder.hh"

//C++ Includes
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

//Root Includes
#include "TCanvas.h"
#include "TH3F.h"
#include "TMath.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector.h"
#include "TVector3.h"

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace shower {
  class TRACSAlg;
}

class shower::TRACSAlg {
public:
  explicit TRACSAlg(const fhicl::ParameterSet& pset);

  void OrderShowerHits(detinfo::DetectorPropertiesData const& detProp,
                       std::vector<art::Ptr<recob::Hit>>& hits,
                       TVector3 const& ShowerDirection,
                       TVector3 const& ShowerPosition) const;

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint>>& showerspcs,
                              TVector3 const& vertex,
                              TVector3 const& direction) const;

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                              TVector3 const& vertex) const;

  TVector3 ShowerCentre(detinfo::DetectorClocksData const& clockData,
                        detinfo::DetectorPropertiesData const& detProp,
                        std::vector<art::Ptr<recob::SpacePoint>> const& showersps,
                        art::FindManyP<recob::Hit> const& fmh,
                        float& totalCharge) const;

  TVector3 ShowerCentre(detinfo::DetectorClocksData const& clockData,
                        detinfo::DetectorPropertiesData const& detProp,
                        std::vector<art::Ptr<recob::SpacePoint>> const& showerspcs,
                        art::FindManyP<recob::Hit> const& fmh) const;

  TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint>> const& showersps) const;

  TVector3 SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const;

  double DistanceBetweenSpacePoints(art::Ptr<recob::SpacePoint> const& sp_a,
                                    art::Ptr<recob::SpacePoint> const& sp_b) const;

  double TotalCorrectedCharge(detinfo::DetectorClocksData const& clockData,
                              detinfo::DetectorPropertiesData const& detProp,
                              std::vector<art::Ptr<recob::SpacePoint>> const& sps,
                              art::FindManyP<recob::Hit> const& fmh) const;

  double SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
                          art::FindManyP<recob::Hit> const& fmh) const;

  double SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
                        art::FindManyP<recob::Hit> const& fmh) const;

  TVector2 HitCoordinates(detinfo::DetectorPropertiesData const& detProp,
                          recob::Hit const& hit) const;

  double SpacePointProjection(art::Ptr<recob::SpacePoint> const& sp,
                              TVector3 const& vertex,
                              TVector3 const& direction) const;

  double SpacePointPerpendiular(art::Ptr<recob::SpacePoint> const& sp,
                                TVector3 const& vertex,
                                TVector3 const& direction,
                                double proj) const;

  void DebugEVD(art::Ptr<recob::PFParticle> const& pfparticle,
                art::Event const& Event,
                reco::shower::ShowerElementHolder& ShowerEleHolder,
                std::string evd_disp_name_append = "") const;

private:
  bool fUseCollectionOnly;
  art::InputTag fHitModuleLabel;
  art::InputTag fPFParticleModuleLabel;
  art::ServiceHandle<geo::Geometry const> fGeom;
  art::ServiceHandle<art::TFileService> tfs;

  std::string fInitialTrackInputLabel;
  std::string fShowerStartPositionInputLabel;
  std::string fShowerDirectionInputLabel;
  std::string fInitialTrackSpacePointsInputLabel;
};

#endif
