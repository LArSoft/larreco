/** ****************************************************************************
 * @file   ClusterCreator.h
 * @brief  Helper functions to create a cluster
 * @date   January 21, 2015
 * @author petrillo@fnal.gov
 * @see    Cluster.h ClusterCreator.cxx
 *
 * ****************************************************************************/

#ifndef CLUSTERCREATOR_H
#define CLUSTERCREATOR_H

// C/C++ standard library
#include <string>
#include <vector>
#include <utility> // std::move()

// LArSoft libraries
#include "lardataobj/RecoBase/Cluster.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"


/// Cluster reconstruction namespace
namespace cluster {

  /** **************************************************************************
   * @brief Class managing the creation of a new recob::Cluster object
   *
   * A Creator is a class that creates a temporary data product, and at the
   * end it yields it to the caller for storage.
   * This last step should be by move construction, although a copy method is
   * also provided.
   *
   * An example of creating a Cluster object:
   * @todo Add the example!
   *
   *     cluster::ClusterCreator cluster(...);
   *     cluster.push_back(cluster.move()); // cluster content becomes invalid
   *
   * This is a one-step creation object: the cluster is constructed at the same
   * time the ClusterCreator is, and no facility is offered to modify the
   * constructed cluster, or to create another one.
   */
  class ClusterCreator {
    public:

      // destructor, copy and move constructor and assignment as default

      /**
       * @brief Constructor: computes some information from hit list
       * @param algo a configured and initialized algorithm set
       * @param start_wire wire coordinate of the start of the cluster
       * @param sigma_start_wire uncertainty on start_wire
       * @param start_tick tick coordinate of the start of the cluster
       * @param sigma_start_tick uncertainty on start_tick
       * @param end_wire wire coordinate of the end of the cluster
       * @param sigma_end_wire uncertainty on end_wire
       * @param end_tick tick coordinate of the end of the cluster
       * @param sigma_end_tick uncertainty on end_tick
       * @param ID cluster ID
       * @param view view for this cluster
       * @param plane location of the start of the cluster
       * @param sentry a check-point variable, optional
       *
       * The algorithm set algo will compute and return:
       * - `start_charge`: charge on the start wire
       * - `start_angle`: angle of the start of the cluster,
       *   in @f$ [ -\pi, \pi ] @f$
       * - `start_opening`: opening angle at the start of the cluster
       * - `sigma_end_tick`: uncertainty on end_tick
       * - `end_charge`: charge on the end wire
       * - `end_angle`: angle of the end of the cluster,
       *   in @f$ [ -\pi, \pi ] @f$
       * - `end_opening`: opening angle at the end of the cluster
       * - `integral`: total charge from fitted shape of hits
       * - `integral_stddev`: standard deviation of hit charge from fitted shape
       * - `summedADC`: total charge from signal ADC of hits
       * - `summedADC_stddev`: standard deviation of signal ADC of hits
       * - `n_hits`: number of hits in the cluster
       * - `multiple_hit_density`: wires covered by cluster divided by number of hits
       * - `width`: a measure of the cluster width
       *
       */
      ClusterCreator(
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
        recob::Cluster::SentryArgument_t sentry = recob::Cluster::Sentry
        );


      /**
       * @brief Prepares the constructed hit to be moved away
       * @return a right-value reference to the constructed hit
       *
       * Despite the name, no move happens in this function.
       * Move takes place in the caller code as proper; for example:
       *
       *     // be cluster a ClusterCreator instance:
       *     std::vector<recob::Cluster> Clusters;
       *     cluster.move();                        // nothing happens
       *     Clusters.push_back(cluster.move());    // here the move happens
       *     recob::Cluster single_cluster(cluster.move()); // wrong! cluster is empty now
       *
       */
      recob::Cluster&& move() { return std::move(cluster); }


      /**
       * @brief Returns the constructed wire
       * @return a constant reference to the constructed wire
       *
       * Despite the name, no copy happens in this function.
       * Copy takes place in the caller code as proper; for example:
       *
       *     // be cluster a ClusterCreator instance:
       *     std::vector<recob::Cluster> Clusters;
       *     cluster.copy();                        // nothing happens
       *     Clusters.push_back(cluster.copy());    // here a copy happens
       *     recob::Cluster single_cluster(cluster.copy()); // copied again
       *
       */
      recob::Cluster const& copy() const { return cluster; }

    protected:

      /// Local instance of the cluster being constructed
      recob::Cluster cluster;

      /**
       * @brief Creates a cluster from direct information and a hit list
       * @param algo a configured and initialized algorithm set
       * @param start_wire wire coordinate of the start of the cluster
       * @param sigma_start_wire uncertainty on start_wire
       * @param start_tick tick coordinate of the start of the cluster
       * @param sigma_start_tick uncertainty on start_tick
       * @param end_wire wire coordinate of the end of the cluster
       * @param sigma_end_wire uncertainty on end_wire
       * @param end_tick tick coordinate of the end of the cluster
       * @param sigma_end_tick uncertainty on end_tick
       * @param ID cluster ID
       * @param view view for this cluster
       * @param plane location of the start of the cluster
       * @param sentry a check-point variable, optional
       *
       * The algorithm set algo will compute and return:
       * - `start_charge`: charge on the start wire
       * - `start_angle`: angle of the start of the cluster,
       *   in @f$ [ -\pi, \pi ] @f$
       * - `start_opening`: opening angle at the start of the cluster
       * - `sigma_end_tick`: uncertainty on end_tick
       * - `end_charge`: charge on the end wire
       * - `end_angle`: angle of the end of the cluster,
       *   in @f$ [ -\pi, \pi ] @f$
       * - `end_opening`: opening angle at the end of the cluster
       * - `integral`: total charge from fitted shape of hits
       * - `integral_stddev`: standard deviation of hit charge from fitted shape
       * - `summedADC`: total charge from signal ADC of hits
       * - `summedADC_stddev`: standard deviation of signal ADC of hits
       * - `n_hits`: number of hits in the cluster
       * - `multiple_hit_density`: wires covered by cluster divided by number of hits
       * - `width`: a measure of the cluster width
       *
       */
      recob::Cluster CreateCluster(
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
        recob::Cluster::SentryArgument_t sentry = recob::Cluster::Sentry
        );

  }; // class ClusterCreator


} // namespace cluster

#endif // CLUSTERCREATOR_H
