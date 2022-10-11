/**
 *  @file   ClusterParamsBuilder.h
 *
 *  @brief  Interface definitions for selecting and building clusters after having run a clustering algorithm
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IClusterParamsBuilder_h
#define IClusterParamsBuilder_h

// Framework Includes
namespace fhicl {
  class ParameterSet;
}

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {
  /**
 *  @brief  ClusterParamsBuilder class definiton
 */
  class IClusterParametersBuilder {
  public:
    /**
     *  @brief  Destructor
     */
    virtual ~IClusterParametersBuilder() noexcept = default;

    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Given the results of running DBScan, format the clusters so that they can be
     *         easily transferred back to the larsoft world
     *
     *  @param hitPairClusterMap      map between view and a list of 3D hits
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *
     *                                The last two parameters are passed through to the IClusterParamsBuilder method
     */
    virtual void BuildClusterInfo(reco::ClusterParametersList& clusterParametersList) const = 0;

    /**
     *  @brief Fill the cluster parameters (expose to outside world for case of splitting/merging clusters)
     *
     *  @param ClusterParameters   The cluster parameters container to be modified
     *  @param Hit2DToClusterMap   Map to keep track of 2D hit to cluster association
     *  @param double              minimum fraction of unique hits
     *  @param double              maximum fraction of "lost" hits
     */

    virtual void FillClusterParams(reco::ClusterParameters&,
                                   reco::Hit2DToClusterMap&,
                                   double minUniqueFrac = 0.,
                                   double maxLostFrac = 1.) const = 0;
  };

} // namespace lar_cluster3d
#endif
