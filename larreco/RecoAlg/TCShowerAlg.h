////////////////////////////////////////////////////////////////////////                                         
// Class:       TCShower                                                                                         
// File:        TCShowerAlg.h                                                                               
//                                                                                                               
// Contact: roryfitz@umich.edu                                                                                   
//                                                                                                               
// module produces showers by selecting tracks surround by many                                                  
// showerLike trajectories as defined by trajcluster with negative                                               
// cluster IDs                                                                                                   
////////////////////////////////////////////////////////////////////////   

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TH1F.h"

#include <memory>

namespace shower {
  class TCShowerAlg {
  public:

    // shower parameters
    TVector3 shwDir;
    TVector3 dcosVtxErr;
    TVector3 shwvtx;
    TVector3 xyzErr;
    std::vector<double> totalEnergy;
    std::vector<double> totalEnergyErr;
    std::vector<double> dEdx;
    std::vector<double> dEdxErr;
    int bestplane;
    
    // vector of shower hits
    std::vector< art::Ptr<recob::Hit> > showerHits;

    TCShowerAlg(fhicl::ParameterSet const& pset);

    int makeShowers(std::vector<art::Ptr<recob::Track> > tracklist, std::vector<art::Ptr<recob::Cluster> > clusterlist, std::vector<art::Ptr<recob::Hit> > hitlist, art::FindManyP<recob::Hit> cls_fm, art::FindManyP<recob::Hit> trk_fm, art::FindManyP<recob::Track> hit_fm, art::FindManyP<recob::Cluster> hitcls_fm, art::FindManyP<anab::Calorimetry> fmcal, art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm);
    
  private: 

    calo::CalorimetryAlg fCalorimetryAlg;

    int goodHit(art::Ptr<recob::Hit>, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2);

    int goodHit(art::Ptr<recob::Hit>, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2, int& pull);

    bool addShowerHit(art::Ptr<recob::Hit> hit, std::vector< art::Ptr<recob::Hit> > showerhits);

  }; // class TCShowerAlg
  
} // namespace shower
