/**
 *  @file   IClusterAlg.h
 * 
 *  @brief  This provides an art tool interface definition for 3D Cluster algorithms
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef IClusterAlg_h
#define IClusterAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief  IClusterAlg interface class definiton
 */
class IClusterAlg
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IClusterAlg() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    virtual void Cluster3DHits(reco::HitPairList&           hitPairList,
                       reco::ClusterParametersList& clusterParametersList) const = 0;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairListPtr        The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    virtual void Cluster3DHits(reco::HitPairListPtr&        hitPairList,
                               reco::ClusterParametersList& clusterParametersList) const = 0;

    /**
     *  @brief enumerate the possible values for time checking if monitoring timing
     */
    enum TimeValues {BUILDTHREEDHITS  = 0,
                     BUILDHITTOHITMAP = 1,
                     RUNDBSCAN        = 2,
                     BUILDCLUSTERINFO = 3,
                     PATHFINDING      = 4,
                     NUMTIMEVALUES
    };
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    virtual float getTimeToExecute(TimeValues index) const = 0;
    
};
    
} // namespace lar_cluster3d
#endif
