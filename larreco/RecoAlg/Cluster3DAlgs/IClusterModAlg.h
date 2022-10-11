/**
 *  @file   IClusterModAlg.h
 *
 *  @brief  This provides an art tool interface definition for 3D Cluster algorithms
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IClusterModAlg_h
#define IClusterModAlg_h

// Framework Includes
namespace fhicl {
  class ParameterSet;
}

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

//------------------------------------------------------------------------------------------------------------------------------------------
namespace art {
  class TFileDirectory;
}

namespace lar_cluster3d {
  /**
 *  @brief  IClusterModAlg interface class definiton
 */
  class IClusterModAlg {
  public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IClusterModAlg() noexcept = default;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Interface for initializing histograms if they are desired
     *         Note that the idea is to put hisgtograms in a subfolder
     *
     *  @param TFileDirectory - the folder to store the hists in
     */
    virtual void initializeHistograms(art::TFileDirectory&) = 0;

    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    virtual void ModifyClusters(reco::ClusterParametersList&) const = 0;

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    virtual float getTimeToExecute() const = 0;
  };

} // namespace lar_cluster3d
#endif
