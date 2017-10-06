/**
 *  @file   ClusterPathFinder_tool.cc
 * 
 *  @brief  art Tool for comparing clusters and merging those that are consistent
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
    
class ClusterPathFinder : virtual public IClusterModAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ClusterPathFinder(const fhicl::ParameterSet&);
    
    /**
     *  @brief  Destructor
     */
    ~ClusterPathFinder();
    
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
    
    /**
     *  @brief Use PCA to try to find path in cluster
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    void breakIntoTinyBits(reco::ClusterParameters&     cluster,
                           reco::ClusterParametersList& outputClusterList) const;

    float closestApproach(const TVector3&, const TVector3&, const TVector3&, const TVector3&, TVector3&, TVector3&) const;
    
    /**
     *  @brief Data members to follow
     */
    bool                                 m_enableMonitoring;      ///<
    size_t                               m_minTinyClusterSize;    ///< Minimum size for a "tiny" cluster
    mutable float                        m_timeToProcess;         ///<
    
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    
    PrincipalComponentsAlg               m_pcaAlg;                // For running Principal Components Analysis
};

ClusterPathFinder::ClusterPathFinder(fhicl::ParameterSet const &pset) :
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterPathFinder::~ClusterPathFinder()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterPathFinder::configure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring   = pset.get<bool>  ("EnableMonitoring",  true  );
    m_minTinyClusterSize = pset.get<size_t>("MinTinyClusterSize",40);

    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    
    m_timeToProcess = 0.;
    
    return;
}
    
void ClusterPathFinder::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
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
    
    int countClusters(0);

    // This is the loop over candidate 3D clusters
    // Note that it might be that the list of candidate clusters is modified by splitting
    // So we use the following construct to make sure we get all of them
    reco::ClusterParametersList::iterator clusterParametersListItr = clusterParametersList.begin();
    
    while(clusterParametersListItr != clusterParametersList.end())
    {
        // Dereference to get the cluster paramters
        reco::ClusterParameters& clusterParameters = *clusterParametersListItr;
        
        // Make sure our cluster has enough hits...
        if (clusterParameters.getHitPairListPtr().size() > m_minTinyClusterSize)
        {
            // get a new output list...
            reco::ClusterParametersList outputClusterList;
        
            // And break our cluster into smaller elements...
            breakIntoTinyBits(clusterParameters, outputClusterList);
        
            std::cout << "**> Broke Cluster " << countClusters++ << " into " << outputClusterList.size() << " sub clusters" << std::endl;
            
            // Add the daughters to the cluster
            clusterParameters.daughterList().insert(clusterParameters.daughterList().end(),outputClusterList.begin(),outputClusterList.end());
        }
    
        // Go to next cluster parameters object
        clusterParametersListItr++;
    }

    if (m_enableMonitoring)
    {
        theClockBuildClusters.stop();
        
        m_timeToProcess = theClockBuildClusters.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> Cluster Path finding done" << std::endl;
    
    return;
}
    
void ClusterPathFinder::breakIntoTinyBits(reco::ClusterParameters&     clusterToBreak,
                                          reco::ClusterParametersList& outputClusterList) const
{
    // This needs to be a recursive routine...
    // Idea is to take the input cluster and order 3D hits by arclength along PCA primary axis
    // If the cluster is above the minimum number of hits then we divide into two and call ourself
    // with the two halves. This means we form a new cluster with hits and PCA and then call ourself
    // If the cluster is below the minimum then we can't break any more, simply add this cluster to
    // the new output list.
    
    // Use the parameters of the incoming cluster to decide what to do
    // Recover the prime ingredients
    reco::PrincipalComponents& fullPCA     = clusterToBreak.getFullPCA();
    std::vector<double>        eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
        3. * std::sqrt(fullPCA.getEigenValues()[1]),
        3. * std::sqrt(fullPCA.getEigenValues()[2])};
    std::vector<TVector3>      pcaAxisVec  = {TVector3(fullPCA.getEigenVectors()[0][0],fullPCA.getEigenVectors()[0][1],fullPCA.getEigenVectors()[0][2]),
        TVector3(fullPCA.getEigenVectors()[1][0],fullPCA.getEigenVectors()[1][1],fullPCA.getEigenVectors()[1][2]),
        TVector3(fullPCA.getEigenVectors()[2][0],fullPCA.getEigenVectors()[2][1],fullPCA.getEigenVectors()[2][2])};
    
    // Look at projections in the planes of interest
    //    TVector3 zAxis(0.,0.,1.);
    //    TVector3 yAxis(0.,1.,0.);
    //    TVector3 widYZProj    = pcaAxisVec[1] - pcaAxisVec[1].Dot(zAxis) * zAxis;
    //    TVector3 heightXZProj = pcaAxisVec[2] - pcaAxisVec[2].Dot(yAxis) * yAxis;
    
    //    float magWidYZProj = eigenValVec[1] * widYZProj.Dot(pcaAxisVec[1]);
    //    float magWidXZProj = eigenValVec[2] * heightXZProj.Dot(pcaAxisVec[2]);
    TVector3 xAxis(1.,0.,0.);
    
    float    cosAng2ToX = std::fabs(pcaAxisVec[2].Dot(xAxis));
    
    // Conditions to split:
    // 1) ratio of secondary spread to maximum spread is more than 4
    // 2) the secondary spread must be "enough"
    // 3)
    //bool moreSplitting = (eigenValVec[1] / eigenValVec[0] > 0.05 && eigenValVec[1] > 0.1) &&  magWidYZProj > 0.1 && magWidXZProj > 0.01 && clusterToBreak.getHitPairListPtr().size() > m_minTinyClusterSize;
    bool moreSplitting = eigenValVec[1] / eigenValVec[0] > 0.02 && eigenValVec[1] > 0.05 && (cosAng2ToX < 0.95 || eigenValVec[2] > 0.05) && clusterToBreak.getHitPairListPtr().size() > m_minTinyClusterSize;
    
    //    std::cout << "===> Cluster with " << clusterToBreak.getHitPairListPtr().size() << " hits, eigen: " << eigenValVec[0] << "/" << eigenValVec[1] << "/" << eigenValVec[2] << ", magWidYZProj: " << magWidYZProj << ", magWidXZProj: " << magWidXZProj << " split: " << moreSplitting << std::endl;
    std::cout << "===> Cluster with " << clusterToBreak.getHitPairListPtr().size() << " hits, eigen: " << eigenValVec[0] << "/" << eigenValVec[1] << "/" << eigenValVec[2] << ", cosAng2ToX: " << cosAng2ToX << " split: " << moreSplitting << std::endl;
    
    // First question, are we done breaking?
    if (!moreSplitting)
    {
        // I think this is where we fill out the rest of the parameters?
        // Start by adding the 2D hits...
        // See if we can avoid duplicates by temporarily transferring to a set
        std::set<const reco::ClusterHit2D*> hitSet;
        
        // Loop through 3D hits to get a set of unique 2D hits
        for(const auto& hit3D : clusterToBreak.getHitPairListPtr())
        {
            for(const auto& hit2D : hit3D->getHits())
            {
                if (hit2D) hitSet.insert(hit2D);
            }
        }
        
        // Now add these to the new cluster
        for(const auto& hit2D : hitSet)
        {
            hit2D->setStatusBit(reco::ClusterHit2D::USED);
            clusterToBreak.UpdateParameters(hit2D);
        }
        
        std::cout << "        --> storing new subcluster of size " << clusterToBreak.getHitPairListPtr().size() << std::endl;
        
        outputClusterList.push_back(clusterToBreak);
    }
    // Otherwise, still chopping...
    else
    {
        // Recover the prime ingredients
        reco::HitPairListPtr& clusHitPairVector = clusterToBreak.getHitPairListPtr();
        
        // Calculate the docas and arc lengths
        m_pcaAlg.PCAAnalysis_calc3DDocas(clusHitPairVector, fullPCA);
        
        // Sort the hits along the PCA
        clusHitPairVector.sort([](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});
        
        // Break into two clusters
        reco::HitPairListPtr::iterator firstElem = clusHitPairVector.begin();
        reco::HitPairListPtr::iterator lastElem  = firstElem;
        
        std::advance(lastElem, clusHitPairVector.size() / 2);
        
        while(firstElem != lastElem)
        {
            reco::ClusterParameters clusterParams;
            
            clusterParams.getHitPairListPtr().resize(std::distance(firstElem,lastElem));
            
            std::copy(firstElem,lastElem,clusterParams.getHitPairListPtr().begin());
            
            // First stage of feature extraction runs here
            m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
            
            // Must have a valid pca
            if (clusterParams.getFullPCA().getSvdOK())
            {
                std::cout << " ******> breakIntoTinyBits has found a valid Full PCA" << std::endl;
                
                // Set the skeleton PCA to make sure it has some value
                clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
                
                breakIntoTinyBits(clusterParams, outputClusterList);
            }
            
            // Set up for next loop (if there is one)
            firstElem = lastElem;
            lastElem  = clusHitPairVector.end();
        }
    }
    
    return;
}

float ClusterPathFinder::closestApproach(const TVector3& P0, const TVector3& u0,
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
    
DEFINE_ART_CLASS_TOOL(ClusterPathFinder)
} // namespace lar_cluster3d
