#ifndef TRACSAlg_hxx
#define TRACSAlg_hxx

//Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "art/Framework/Principal/Event.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/ShowerElementHolder.hh"

//C++ Includes
#include <iostream>
#include <vector>
#include <map>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TString.h"

namespace shower {
  class TRACSAlg;
}

class shower::TRACSAlg {
  public:
    TRACSAlg(const fhicl::ParameterSet& pset);

    void OrderShowerHits(std::vector<art::Ptr<recob::Hit> >& hits,
        TVector3 const& ShowerDirection,
        TVector3 const& ShowerPosition
        ) const;

    void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& showerspcs,
        TVector3 const& vertex, TVector3 const& direction) const;


    TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> > const& showersps,
        art::FindManyP<recob::Hit> const& fmh, float& totalCharge) const;


    TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> > const& showerspcs,
        art::FindManyP<recob::Hit> const& fmh) const;

    TVector3 SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const;

    double SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp, art::FindManyP<recob::Hit> const& fmh) const;

    double SpacePointTime(art::Ptr<recob::SpacePoint> const& sp, art::FindManyP<recob::Hit> const& fmh) const;

    TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit) const;


    double SpacePointProjection(art::Ptr<recob::SpacePoint> const& sp, TVector3 const& vertex,
        TVector3 const& direction) const;

    double SpacePointPerpendiular(art::Ptr<recob::SpacePoint> const& sp, TVector3 const& vertex,
        TVector3 const& direction, double proj) const;

    void DebugEVD(art::Ptr<recob::PFParticle> const& pfparticle,
        art::Event const& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder) const;

  private:

    bool fUseCollectionOnly;
    art::InputTag                           fHitModuleLabel;
    art::InputTag                           fPFParticleModuleLabel;
    detinfo::DetectorProperties const*      fDetProp = nullptr;
    art::ServiceHandle<geo::Geometry const> fGeom;
    art::ServiceHandle<art::TFileService>   tfs;

};

#endif
