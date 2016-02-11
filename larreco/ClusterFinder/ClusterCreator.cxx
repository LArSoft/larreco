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

// C/C++ standard library
#include <utility> // std::move()
#include <algorithm> // std::accumulate()

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"


/// Cluster reconstruction namespace
namespace cluster {
  
  //****************************************************************************
  //***  ClusterCreator
  //----------------------------------------------------------------------
  ClusterCreator::ClusterCreator(
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
    ):
    cluster(CreateCluster(
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
      sentry
      ))
  {} // ClusterCreator::ClusterCreator()
  
  
  //----------------------------------------------------------------------
  recob::Cluster ClusterCreator::CreateCluster(
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
  ) {
    return recob::Cluster(
      start_wire,
      sigma_start_wire,
      start_tick,
      sigma_start_tick,
      algo.StartCharge().value(),       // start_charge
      algo.StartAngle().value(),        // start_angle
      algo.StartOpeningAngle().value(), // start_opening
      end_wire,
      sigma_end_wire,
      end_tick,
      sigma_end_tick,
      algo.EndCharge().value(),         // end_charge
      algo.EndAngle().value(),          // end_angle
      algo.EndOpeningAngle().value(),   // end_opening
      algo.Integral().value(),          // integral
      algo.IntegralStdDev().value(),    // integral_stddev
      algo.SummedADC().value(),         // summedADC
      algo.SummedADCStdDev().value(),   // summedADC_stddev
      algo.NHits(),                     // n_hits
      algo.MultipleHitDensity(),           // multiple_hit_density
      algo.Width(),                     // width
      ID,
      view,
      plane,
      sentry
      );
  } // ClusterCreator::CreateCluster()
  
  
  //----------------------------------------------------------------------
} // namespace cluster
