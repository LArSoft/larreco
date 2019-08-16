#ifndef TRACSAlg_hxx
#define TRACSAlg_hxx

//Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
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
#include "larreco/ShowerFinder/ShowerElementHolder.hh"

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
        TVector3& ShowerDirection,
        TVector3& ShowerPosition
        );

    void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& showerspcs,
        TVector3& vertex, TVector3& direction);


    TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showersps,
        art::FindManyP<recob::Hit>& fmh, float& totalCharge);


    TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showerspcs,
        art::FindManyP<recob::Hit>& fmh);

    TVector3 SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp);

    double SpacePointCharge(art::Ptr<recob::SpacePoint> sp, art::FindManyP<recob::Hit>& fmh);

    double SpacePointTime(art::Ptr<recob::SpacePoint> sp, art::FindManyP<recob::Hit>& fmh);

    TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);


    double SpacePointProjection(const art::Ptr<recob::SpacePoint>&sp, TVector3& vertex,
        TVector3& direction);

    double SpacePointPerpendiular(const art::Ptr<recob::SpacePoint>&sp, TVector3& vertex,
        TVector3& direction, double proj);

    void DebugEVD(const art::Ptr<recob::PFParticle>& pfparticle,
        art::Event& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder);

  private:

    bool fUseCollectionOnly;
    art::InputTag                           fHitModuleLabel;
    art::InputTag                           fPFParticleModuleLabel;
    detinfo::DetectorProperties const*      fDetProp = nullptr;
    art::ServiceHandle<geo::Geometry const> fGeom;
    art::ServiceHandle<art::TFileService>   tfs;

};

#endif
