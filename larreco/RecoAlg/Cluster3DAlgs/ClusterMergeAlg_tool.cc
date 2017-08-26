/**
 *  @file   ClusterMergeAlg.cxx
 * 
 *  @brief  Algorithm for comparing clusters and merging those that are consistent
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// Root includes
#include "TVector3.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {
    
class ClusterMergeAlg : virtual public IClusterModAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ClusterMergeAlg(const fhicl::ParameterSet&);
    
    /**
     *  @brief  Destructor
     */
    ~ClusterMergeAlg();
    
    void configure(fhicl::ParameterSet const &pset) override;
    
    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void ModifyClusters(reco::ClusterParametersList&) const override;
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute() const override {return m_timeToProcess;}
    
private:
    
    bool consistentClusters(const reco::PrincipalComponents&, const reco::PrincipalComponents&) const;
    
    void mergeClusters(reco::ClusterParameters&, reco::ClusterParameters&) const;
    
    /**
     *  @brief Data members to follow
     */
    bool                                 m_enableMonitoring;      ///<
    double                               m_minCosAxisAng;         ///<
    mutable float                        m_timeToProcess;         ///<
    
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    
    PrincipalComponentsAlg               m_pcaAlg;                // For running Principal Components Analysis
};

ClusterMergeAlg::ClusterMergeAlg(fhicl::ParameterSet const &pset) :
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergeAlg::~ClusterMergeAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergeAlg::configure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring = pset.get<bool>  ("EnableMonitoring",  true  );
    m_minCosAxisAng    = pset.get<double>("MinCosAxisAng",     0.95  );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    
    m_timeToProcess = 0.;
    
    return;
}
    
void ClusterMergeAlg::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
{
    /**
     *  @brief Top level interface for algorithm to consider pairs of clusters from the input
     *         list and determine if they are consistent with each other and, therefore, should
     *         be merged. This is done by looking at the PCA for each cluster and looking at the
     *         projection of the primary axis along the vector connecting their centers.
     */
    
    // Initial clustering is done, now trim the list and get output parameters
    cet::cpu_timer theClockBuildClusters;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockBuildClusters.start();
    
    reco::ClusterParametersList::iterator firstClusterItr = clusterParametersList.begin();
    
    while(firstClusterItr != clusterParametersList.end())
    {
        reco::ClusterParameters&              firstClusterParams = *firstClusterItr;
        reco::ClusterParametersList::iterator nextClusterItr     = firstClusterItr;
        
        // want the next one...
        nextClusterItr++;
    
        while(nextClusterItr != clusterParametersList.end())
        {
            reco::ClusterParameters& nextClusterParams = *nextClusterItr;
            
            // On any given loop through here it **might** be that the first cluster has been modified. So can't cache
            // the parameters, need the curret ones
            if (consistentClusters(firstClusterParams.getFullPCA(),nextClusterParams.getFullPCA()))
            {
                mergeClusters(firstClusterParams, nextClusterParams);
                
                // Now remove the "next" cluster
                nextClusterItr = clusterParametersList.erase(nextClusterItr);
            }
            else nextClusterItr++;
        }
        
        firstClusterItr++;
    }
    
    
    if (m_enableMonitoring)
    {
        theClockBuildClusters.stop();
        
        m_timeToProcess = theClockBuildClusters.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> Merge clusters done, found " << clusterParametersList.size() << " clusters" << std::endl;
    
    return;
}
    
bool ClusterMergeAlg::consistentClusters(const reco::PrincipalComponents& firstPCA, const reco::PrincipalComponents& nextPCA) const
{
    // Assume failure
    bool consistent(false);
    
    // Simplest check first, require some directionality agreement
    TVector3 firstAxis0(firstPCA.getEigenVectors()[0][0],firstPCA.getEigenVectors()[0][1],firstPCA.getEigenVectors()[0][2]);
    TVector3 nextAxis0(nextPCA.getEigenVectors()[0][0],nextPCA.getEigenVectors()[0][1],nextPCA.getEigenVectors()[0][2]);

    if (std::fabs(firstAxis0.Dot(nextAxis0)) > m_minCosAxisAng)
    {
        // Recover the position and direction of the first cluster
        TVector3 firstCenter(firstPCA.getAvePosition()[0],firstPCA.getAvePosition()[1],firstPCA.getAvePosition()[2]);
        TVector3 firstAxis1(firstPCA.getEigenVectors()[1][0],firstPCA.getEigenVectors()[1][1],firstPCA.getEigenVectors()[1][2]);
        TVector3 firstAxis2(firstPCA.getEigenVectors()[2][0],firstPCA.getEigenVectors()[2][1],firstPCA.getEigenVectors()[2][2]);
    
        // Get the eigen values
        TVector3 firstEigenVals(3. * sqrt(firstPCA.getEigenValues()[0]),
                                3. * sqrt(firstPCA.getEigenValues()[1]),
                                3. * sqrt(firstPCA.getEigenValues()[2]));
    
        // Recover the position of the next cluster
        TVector3 nextCenter(nextPCA.getAvePosition()[0],nextPCA.getAvePosition()[1],nextPCA.getAvePosition()[2]);
    
        // And form a vector between the two centers
        TVector3 firstPosToNextPos = nextCenter - firstCenter;

        // First check will be that the center of the second cluster lies within the cone of the projection of the first cluster
        // Start by getting the arc length to the doca of the second cluster center to the vector from the first
        double arcLenToDoca = firstPosToNextPos.Dot(firstAxis0);
    
        // Get coordinates of the point on this axis that is closest to the next cluster center
        TVector3 firstDocaPos = firstCenter + arcLenToDoca * firstAxis0;
    
        // Now get the vector from here to the next cluster center
        TVector3 firstDocaPosToNext = nextCenter - firstDocaPos;
    
        // Get projections of this vector onto the other two PCA axes for the first cluster
        double docaProj1 = firstAxis1.Dot(firstDocaPosToNext);
        double docaProj2 = firstAxis2.Dot(firstDocaPosToNext);
    
        // Get the ratio of the first eigen value to the arc length
        double scaleFactor = std::fabs(arcLenToDoca / firstEigenVals[0]);
        double cutValue    = 3. * scaleFactor * std::max(firstEigenVals[1],firstEigenVals[2]);
    
        // Is the next cluster position consistent with the first?
        if (std::fabs(docaProj1) < cutValue && std::fabs(docaProj2) < cutValue)
        {
            // Final check on colinearity
            // Get the distance from the first to the next
            double firstNextDistance = firstPosToNextPos.Mag();
        
            // Now make unit vector...
            firstPosToNextPos.SetMag(1.);
        
            // Get the projection of the first axis on this one
            double firstProj = firstEigenVals[0] * std::fabs(firstPosToNextPos.Dot(firstAxis0));
            
            // Get the projection of the next axis on this one
            double nextProj = nextPCA.getEigenValues()[0] * std::fabs(firstPosToNextPos.Dot(nextAxis0));
        
            // Require both projects are a significant fraction of distance
            if (firstProj + nextProj > 0.5 * firstNextDistance)
            {
                consistent = true;
            }
        }
    }
    
    return consistent;
}
    
void ClusterMergeAlg::mergeClusters(reco::ClusterParameters& firstClusterParams, reco::ClusterParameters& nextClusterParams) const
{
    // Merge the next cluster into the first one
    // Start by copying the hits and edges from the next cluster to the first one
    reco::HitPairListPtr& firstHitPairListPtr = firstClusterParams.getHitPairListPtr();
    reco::HitPairListPtr& nextHitPairListPtr  = nextClusterParams.getHitPairListPtr();
    
    for(const auto* hit : nextHitPairListPtr)
    {
        firstHitPairListPtr.push_back(hit);
        
        for(const auto* hit2D : hit->getHits())
            if (hit2D) firstClusterParams.UpdateParameters(hit2D);
    }
    
    reco::Hit3DToEdgeMap& firstEdgeMap = firstClusterParams.getHit3DToEdgeMap();
    reco::Hit3DToEdgeMap& nextEdgeMap  = nextClusterParams.getHit3DToEdgeMap();
    
    for(const auto& pair : nextEdgeMap) firstEdgeMap[pair.first] = pair.second;
    
    // Recalculate the PCA
    m_pcaAlg.PCAAnalysis_3D(firstClusterParams.getHitPairListPtr(), firstClusterParams.getFullPCA());
    
    // Must have a valid pca
    if (firstClusterParams.getFullPCA().getSvdOK())
    {
        // Set the skeleton PCA to make sure it has some value
        firstClusterParams.getSkeletonPCA() = firstClusterParams.getFullPCA();
    }
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(ClusterMergeAlg)
} // namespace lar_cluster3d
