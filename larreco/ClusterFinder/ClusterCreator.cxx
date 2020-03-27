/** ****************************************************************************
 * @file   ClusterCreator.cxx
 * @brief  Helper functions to create a cluster - implementation file
 * @date   January 21, 2015
 * @author petrillo@fnal.gov
 * @see    Cluster.h ClusterCreator.h
 *
 * ****************************************************************************/

// declaration header
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"

// LArSoft libraries
#include "lardataobj/RecoBase/Cluster.h"

namespace cluster {

  ClusterCreator::ClusterCreator(
    util::GeometryUtilities const& gser,
    ClusterParamsAlgBase& algo,
    float start_wire,
    float sigma_start_wire,
    float start_tick,
    float sigma_start_tick,
    float end_wire,
    float sigma_end_wire,
    float end_tick,
    float sigma_end_tick,
    recob::Cluster::ID_t ID,
    geo::View_t view,
    geo::PlaneID const& plane,
    recob::Cluster::SentryArgument_t sentry /* = recob::Cluster::Sentry */
    )
    : cluster(CreateCluster(gser,
                            algo,
                            start_wire,
                            sigma_start_wire,
                            start_tick,
                            sigma_start_tick,
                            end_wire,
                            sigma_end_wire,
                            end_tick,
                            sigma_end_tick,
                            ID,
                            view,
                            plane,
                            sentry))
  {} // ClusterCreator::ClusterCreator()

  //----------------------------------------------------------------------
  recob::Cluster
  ClusterCreator::CreateCluster(
    util::GeometryUtilities const& gser,
    ClusterParamsAlgBase& algo,
    float start_wire,
    float sigma_start_wire,
    float start_tick,
    float sigma_start_tick,
    float end_wire,
    float sigma_end_wire,
    float end_tick,
    float sigma_end_tick,
    recob::Cluster::ID_t ID,
    geo::View_t view,
    geo::PlaneID const& plane,
    recob::Cluster::SentryArgument_t sentry /* = recob::Cluster::Sentry */
  )
  {
    return recob::Cluster(start_wire,
                          sigma_start_wire,
                          start_tick,
                          sigma_start_tick,
                          algo.StartCharge(gser).value(),
                          algo.StartAngle().value(),
                          algo.StartOpeningAngle().value(),
                          end_wire,
                          sigma_end_wire,
                          end_tick,
                          sigma_end_tick,
                          algo.EndCharge(gser).value(),
                          algo.EndAngle().value(),
                          algo.EndOpeningAngle().value(),
                          algo.Integral().value(),
                          algo.IntegralStdDev().value(),
                          algo.SummedADC().value(),
                          algo.SummedADCStdDev().value(),
                          algo.NHits(),
                          algo.MultipleHitDensity(),
                          algo.Width(gser),
                          ID,
                          view,
                          plane,
                          sentry);
  } // ClusterCreator::CreateCluster()

} // namespace cluster
