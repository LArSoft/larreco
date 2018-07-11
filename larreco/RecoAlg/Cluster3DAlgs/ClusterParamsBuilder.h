/**
 *  @file   ClusterParamsBuilder.h
 * 
 *  @brief  This algorithm will create and then cluster 3D hits using DBScan
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef ClusterParamsBuilder_h
#define ClusterParamsBuilder_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

// std includes
#include <vector>
#include <list>
#include <set>
#include <map>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief  ClusterParamsBuilder class definiton
 */
class ClusterParamsBuilder
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    ClusterParamsBuilder(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~ClusterParamsBuilder();

    void configure(const fhicl::ParameterSet&);
    
    /**
     *  @brief Given the results of running DBScan, format the clusters so that they can be
     *         easily transferred back to the larsoft world
     *
     *  @param hitPairClusterMap      map between view and a list of 3D hits
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *
     *                                The last two parameters are passed through to the FillClusterParams method
     */
    void BuildClusterInfo(reco::ClusterParametersList& clusterParametersList) const;
    
    /**
     *  @brief Fill the cluster parameters (expose to outside world for case of splitting/merging clusters)
     *
     *  @param ClusterParameters   The cluster parameters container to be modified
     *  @param Hit2DToClusterMap   Map to keep track of 2D hit to cluster association
     *  @param double              minimum fraction of unique hits
     *  @param double              maximum fraction of "lost" hits
     */
    
    void FillClusterParams(reco::ClusterParameters&,
                           reco::Hit2DToClusterMap&,
                           double minUniqueFrac = 0.,
                           double maxLostFrac=1.) const;
    
private:
    
    void removeUsedHitsFromMap(reco::ClusterParameters&, reco::HitPairListPtr&, reco::Hit2DToClusterMap&) const;
    
    /**
     *  @brief Data members to follow
     */
    size_t                 m_clusterMinHits;
    double                 m_clusterMinUniqueFraction;
    double                 m_clusterMaxLostFraction;
    
    PrincipalComponentsAlg m_pcaAlg;                // For running Principal Components Analysis
};
    
} // namespace lar_cluster3d
#endif
