#include "ShowerRecoAlgBase.h"

#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParams.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include "larreco/ShowerFinder/ShowerReco3D/ShowerRecoException.h"

namespace showerreco {

  void ShowerRecoAlgBase::Reset() { fInputClusters.clear(); }

  void ShowerRecoAlgBase::AppendInputClusters(
    const std::vector<::cluster::ClusterParamsAlg>& cpan_v)
  {
    std::vector<::showerreco::ShowerCluster_t> clusters;
    clusters.reserve(cpan_v.size());

    for (auto const& cpan : cpan_v) {

      clusters.push_back(::showerreco::ShowerCluster_t());

      (*clusters.rbegin()).start_point = cpan.GetParams().start_point;
      (*clusters.rbegin()).end_point = cpan.GetParams().end_point;
      (*clusters.rbegin()).angle_2d = cpan.GetParams().angle_2d;
      (*clusters.rbegin()).plane_id = cpan.Plane();
      (*clusters.rbegin()).hit_vector = cpan.GetHitVector();
    }

    fInputClusters.push_back(clusters);
  }

  std::vector<::recob::Shower> ShowerRecoAlgBase::Reconstruct(
    geo::GeometryCore const& geom,
    geo::WireReadoutGeom const& wireReadoutGeom,
    detinfo::DetectorClocksData const& clockData,
    detinfo::DetectorPropertiesData const& detProp)
  {
    ProcessInputClusters();

    std::vector<::recob::Shower> output;
    output.reserve(fInputClusters.size());

    for (auto const& clusters : fInputClusters)
      output.push_back(RecoOneShower(geom, wireReadoutGeom, clockData, detProp, clusters));

    return output;
  }

}
