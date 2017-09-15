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
    
    float closestApproach(const TVector3&, const TVector3&, const TVector3&, const TVector3&, TVector3&, TVector3&) const;
    
    /**
     *  @brief Data members to follow
     */
    bool                                 m_enableMonitoring;      ///<
    double                               m_minCosAxisAng;         ///< minimum Cos(angle) cut value
    double                               m_minEigenToProcess;     ///< Don't look anymore at clusters below this size
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
    m_enableMonitoring  = pset.get<bool>  ("EnableMonitoring",  true  );
    m_minCosAxisAng     = pset.get<double>("MinCosAxisAng",     0.975 );
    m_minEigenToProcess = pset.get<double>("MinEigenToProcess", 5.0   );
    
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
    
    // Resort by the largest first eigen value (so the "longest" cluster)
    clusterParametersList.sort([](auto& left, auto& right){return left.getFullPCA().getEigenValues()[0] > right.getFullPCA().getEigenValues()[0];});
    
    reco::ClusterParametersList::iterator firstClusterItr = clusterParametersList.begin();
    
    while(firstClusterItr != clusterParametersList.end())
    {
        reco::ClusterParameters&              firstClusterParams = *firstClusterItr;
        reco::ClusterParametersList::iterator nextClusterItr     = firstClusterItr;
        
        // Once you get down to the smallest clusters if they haven't already been absorbed there is no need to check them
        if (firstClusterParams.getFullPCA().getEigenValues()[0] < m_minEigenToProcess) break;
        
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
    
    // Recover the positions of the centers of the two clusters
    TVector3 firstCenter(firstPCA.getAvePosition()[0],firstPCA.getAvePosition()[1],firstPCA.getAvePosition()[2]);
    TVector3 nextCenter(nextPCA.getAvePosition()[0],nextPCA.getAvePosition()[1],nextPCA.getAvePosition()[2]);
    
    // And form a vector between the two centers
    TVector3 firstPosToNextPos = nextCenter - firstCenter;
    
    // Now get the first PCA's primary axis
    TVector3 firstAxis0(firstPCA.getEigenVectors()[0][0],firstPCA.getEigenVectors()[0][1],firstPCA.getEigenVectors()[0][2]);
    
    // Get the arclength to the point on the first axis which is the poca to the next point
    double arcLenToNextDoca = firstPosToNextPos.Dot(firstAxis0);
    
    // And recover that point
    TVector3 firstAxisPoca = firstCenter + arcLenToNextDoca * firstAxis0;
    
    // Get the eigen values
    TVector3 firstEigenVals(3. * sqrt(firstPCA.getEigenValues()[0]),
                            3. * sqrt(firstPCA.getEigenValues()[1]),
                            3. * sqrt(firstPCA.getEigenValues()[2]));
    
    // Get the ratio of the first eigen value to the arc length
    double scaleFactor = std::fabs(arcLenToNextDoca / firstEigenVals[0]);
    double cutValue    = 3. * scaleFactor * firstEigenVals[1];
    
    // Ok, the first selection is that the doca of the next point to the first PCA primary axis is within
    // the allowed maximum spread
    if ((nextCenter - firstAxisPoca).Mag() < cutValue && std::fabs(arcLenToNextDoca) < 6. * firstEigenVals[0])
    {
        // Now set up to look at the doca
        TVector3 firstPoca;
        TVector3 nextPoca;
        
        // Get next primary axis
        TVector3 nextAxis0(nextPCA.getEigenVectors()[0][0],nextPCA.getEigenVectors()[0][1],nextPCA.getEigenVectors()[0][2]);
        
        // Recover the doca of the two axes and their points of closest approach
        float lineDoca      = closestApproach(firstCenter, firstAxis0, nextCenter, nextAxis0, firstPoca, nextPoca);
        
        // Effectively the same cut as above but now with the point of closest approach
        float arcLenToPoca0 = (firstPoca - firstCenter).Mag();
        float docaCutValue  = std::fabs(arcLenToPoca0) * firstEigenVals[1] / firstEigenVals[0];
        
        if (lineDoca < 3. * docaCutValue)
        {
            // Final check, look at arc lengths to pocas
            double arcLenToPocaFirst = (firstPoca - firstCenter).Mag();
            double arcLenToPocaNext  = (nextPoca  - nextCenter).Mag();
            double totalArcLen       = std::fabs(arcLenToPocaFirst) + std::fabs(arcLenToPocaNext);
            
            if (2. * totalArcLen >firstPosToNextPos.Mag() && totalArcLen < firstPosToNextPos.Mag() / m_minCosAxisAng) consistent = true;
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
    
float ClusterMergeAlg::closestApproach(const TVector3& P0, const TVector3& u0,
                                       const TVector3& P1, const TVector3& u1,
                                       TVector3&       poca0,
                                       TVector3&       poca1) const
{
    // Technique is to compute the arclength to each point of closest approach
    TVector3 w0 = P0 - P1;
    float a(1.);
    float b(u0.Dot(u1));
    float c(1.);
    float d(u0.Dot(w0));
    float e(u1.Dot(w0));
    float den(a * c - b * b);
    
    float arcLen0 = (b * e - c * d) / den;
    float arcLen1 = (a * e - b * d) / den;
    
    poca0 = P0 + arcLen0 * u0;
    poca1 = P1 + arcLen1 * u1;
    
    return (poca0 - poca1).Mag();
}
    
DEFINE_ART_CLASS_TOOL(ClusterMergeAlg)
} // namespace lar_cluster3d
