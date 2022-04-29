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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"

#include "TVector3.h"

#include <map>
#include <vector>

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

    int makeShowers(detinfo::DetectorClocksData const& dataClock,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector<art::Ptr<recob::PFParticle>> const& pfplist,
                    std::vector<art::Ptr<recob::Vertex>> const& vertexlist,
                    std::vector<art::Ptr<recob::Cluster>> const& clusterlist,
                    std::vector<art::Ptr<recob::Hit>> const& hitlist,
                    art::FindManyP<recob::Hit> const& cls_fm,
                    art::FindManyP<recob::Cluster> const& clspfp_fm,
                    art::FindManyP<recob::Vertex> const& vtxpfp_fm,
                    art::FindManyP<recob::PFParticle> const& hit_fm,
                    art::FindManyP<recob::Cluster> const& hitcls_fm,
                    art::FindManyP<recob::Track> const& trkpfp_fm,
                    art::FindManyP<anab::Calorimetry> const& fmcal);

  private:
    calo::CalorimetryAlg fCalorimetryAlg;
    pma::ProjectionMatchingAlg fProjectionMatchingAlg;

    int goodHit(detinfo::DetectorClocksData const& dataClock,
                detinfo::DetectorPropertiesData const& detProp,
                art::Ptr<recob::Hit> const&,
                double maxDist,
                double minDistVert,
                std::map<geo::PlaneID, double> const& trk_wire1,
                std::map<geo::PlaneID, double> const& trk_tick1,
                std::map<geo::PlaneID, double> const& trk_wire2,
                std::map<geo::PlaneID, double> const& trk_tick2) const;

    int goodHit(detinfo::DetectorClocksData const& dataClock,
                detinfo::DetectorPropertiesData const& detProp,
                art::Ptr<recob::Hit> const&,
                double maxDist,
                double minDistVert,
                std::map<geo::PlaneID, double> const& trk_wire1,
                std::map<geo::PlaneID, double> const& trk_tick1,
                std::map<geo::PlaneID, double> const& trk_wire2,
                std::map<geo::PlaneID, double> const& trk_tick2,
                int& pull) const;

    bool addShowerHit(art::Ptr<recob::Hit> hit, std::vector<art::Ptr<recob::Hit>> showerhits) const;

  }; // class TCShowerAlg

} // namespace shower
