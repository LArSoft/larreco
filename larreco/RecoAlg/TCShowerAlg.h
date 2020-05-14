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

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
namespace fhicl {
  class ParameterSet;
}

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"

#include "TVector3.h"

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
    std::vector<art::Ptr<recob::Hit>> showerHits;

    TCShowerAlg(fhicl::ParameterSet const& pset);

    int makeShowers(std::vector<art::Ptr<recob::PFParticle>> pfplist,
                    std::vector<art::Ptr<recob::Vertex>> vertexlist,
                    std::vector<art::Ptr<recob::Cluster>> clusterlist,
                    std::vector<art::Ptr<recob::Hit>> hitlist,
                    art::FindManyP<recob::Hit> cls_fm,
                    art::FindManyP<recob::Cluster> clspfp_fm,
                    art::FindManyP<recob::Vertex> vtxpfp_fm,
                    art::FindManyP<recob::PFParticle> hit_fm,
                    art::FindManyP<recob::Cluster> hitcls_fm,
                    art::FindManyP<recob::Track> trkpfp_fm,
                    art::FindManyP<anab::Calorimetry> fmcal);

  private:
    calo::CalorimetryAlg fCalorimetryAlg;
    pma::ProjectionMatchingAlg fProjectionMatchingAlg;

    int goodHit(art::Ptr<recob::Hit>,
                double maxDist,
                double minDistVert,
                std::map<geo::PlaneID, double> trk_wire1,
                std::map<geo::PlaneID, double> trk_tick1,
                std::map<geo::PlaneID, double> trk_wire2,
                std::map<geo::PlaneID, double> trk_tick2) const;

    int goodHit(art::Ptr<recob::Hit>,
                double maxDist,
                double minDistVert,
                std::map<geo::PlaneID, double> trk_wire1,
                std::map<geo::PlaneID, double> trk_tick1,
                std::map<geo::PlaneID, double> trk_wire2,
                std::map<geo::PlaneID, double> trk_tick2,
                int& pull) const;

    bool addShowerHit(art::Ptr<recob::Hit> hit, std::vector<art::Ptr<recob::Hit>> showerhits) const;

  }; // class TCShowerAlg

} // namespace shower
