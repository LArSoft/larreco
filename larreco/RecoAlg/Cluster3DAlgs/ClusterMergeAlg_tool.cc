/**
 *  @file   ClusterMergeAlg.cxx
 * 
 *  @brief  Algorithm for comparing clusters and merging those that are consistent
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
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
     *  @brief Interface for initializing histograms if they are desired
     *         Note that the idea is to put hisgtograms in a subfolder
     *
     *  @param TFileDirectory - the folder to store the hists in
     */
    void initializeHistograms(art::TFileDirectory&) override;

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
    
    bool mergeClusters(reco::ClusterParameters&, reco::ClusterParameters&) const;
    
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
    m_minEigenToProcess = pset.get<double>("MinEigenToProcess", 2.0   );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    
    m_timeToProcess = 0.;
    
    return;
}
    
void ClusterMergeAlg::initializeHistograms(art::TFileDirectory&)
{
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
    
    int clusCntr(0);
    
    while(firstClusterItr != clusterParametersList.end())
    {
        reco::ClusterParameters&              firstClusterParams = *firstClusterItr;
        reco::ClusterParametersList::iterator nextClusterItr     = firstClusterItr;
        
        // Once you get down to the smallest clusters if they haven't already been absorbed there is no need to check them
        if (firstClusterParams.getFullPCA().getEigenValues()[0] < m_minEigenToProcess) break;
        
        std::cout << "+++++++++++++++++++++++++++++++ Checking PCA for cluster # " << clusCntr++ << " +++++++++++++++++++++++++++" << std::endl;
        std::cout << "+++++++ eigen values: " << firstClusterParams.getFullPCA().getEigenValues()[0] << "/" << firstClusterParams.getFullPCA().getEigenValues()[1] << "/" << firstClusterParams.getFullPCA().getEigenValues()[2] << " +++++++++" << std::endl;
        
        // want the next one...
        nextClusterItr++;
    
        while(nextClusterItr != clusterParametersList.end())
        {
            reco::ClusterParameters& nextClusterParams = *nextClusterItr;
            
            // On any given loop through here it **might** be that the first cluster has been modified. So can't cache
            // the parameters, need the curret ones
            if (consistentClusters(firstClusterParams.getFullPCA(),nextClusterParams.getFullPCA()))
            {
                if (mergeClusters(firstClusterParams, nextClusterParams))
                {
                    // Now remove the "next" cluster
                    nextClusterItr = clusterParametersList.erase(nextClusterItr);
                
                    // Restart loop?
                    nextClusterItr = firstClusterItr;
                }
            }
            
            nextClusterItr++;
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
    
    // Two types of conditions to decide between to start:
    // 1) the cluster to merge lies along the trajectory of the primary cluster, in which case merge it
    //    - note that this should take care of relatively colinear trajectories
    // 2) the trajectories of the two clusters are consistent
    //    - their doca is small
    //    - the vectors to their centers are not inconsistent with the trajectory of the first,
    //    - etc.
    //
    // Initial set up to check if the merge candidate is within the trajectory of the primary cluster
    
    // Recover the positions of the centers of the two clusters
    TVector3 firstCenter(firstPCA.getAvePosition()[0],firstPCA.getAvePosition()[1],firstPCA.getAvePosition()[2]);
    TVector3 nextCenter(nextPCA.getAvePosition()[0],nextPCA.getAvePosition()[1],nextPCA.getAvePosition()[2]);
    
    // And form a vector between the two centers
    TVector3 firstPosToNextPos = nextCenter - firstCenter;
    
    // Now get the first PCA's primary axis and since we'll use them get all of them at once...
    TVector3 firstAxis0(firstPCA.getEigenVectors()[0][0],firstPCA.getEigenVectors()[0][1],firstPCA.getEigenVectors()[0][2]);
    TVector3 firstAxis1(firstPCA.getEigenVectors()[1][0],firstPCA.getEigenVectors()[1][1],firstPCA.getEigenVectors()[1][2]);
    TVector3 firstAxis2(firstPCA.getEigenVectors()[2][0],firstPCA.getEigenVectors()[2][1],firstPCA.getEigenVectors()[2][2]);
    
    // Adopt the convention that the cluster axis is in same direction as vector from first to next centers
    if (firstPosToNextPos.Dot(firstAxis0) < 0.) firstAxis0 = -firstAxis0;

    // Want the distance of closest approach of the next cluser's center to the primary, start by finding arc length
    double arcLenToNextDoca = firstPosToNextPos.Dot(firstAxis0);
    
    // Position on the first cluster's axis of doca to next center
    TVector3 firstAxisDocaPos = firstCenter + arcLenToNextDoca * firstAxis0;
    
    // Doca vector
    TVector3 nextDocaVec = nextCenter - firstAxisDocaPos;
    
    // Need the projection of the doca vector onto the two transverse axes of the first cluster
    double docaVecProj1 = std::fabs(firstAxis1.Dot(nextDocaVec));
    double docaVecProj2 = std::fabs(firstAxis2.Dot(nextDocaVec));
    
    // We will now compare these to the eigen values of the first cluster, so recover all of them
    TVector3 firstEigenVals(3. * sqrt(firstPCA.getEigenValues()[0]),
                            3. * sqrt(firstPCA.getEigenValues()[1]),
                            3. * sqrt(firstPCA.getEigenValues()[2]));
    
    // Use the angle between the vector between cluster centers and the first axis to moderate the selection cut
    double firstToNextDist  = firstPosToNextPos.Mag();
    double cosAngFTNtoAxis0 = arcLenToNextDoca / firstToNextDist;
    double docaVecProj1Cut  = std::max( 1., (1. + 2. * cosAngFTNtoAxis0) * firstEigenVals[1]);
    double docaVecProj2Cut  = std::max(0.5, (1. + 2. * cosAngFTNtoAxis0) * firstEigenVals[2]);
    
    std::cout << "   ==> Check in tube, doca: " << nextDocaVec.Mag() << ", proj: " << docaVecProj1 << "/" << docaVecProj2 << ", cut: " << docaVecProj1Cut << "/" << docaVecProj2Cut << ", eigen: " << firstEigenVals[0] << "/" << firstEigenVals[1] << "/" << firstEigenVals[2] << ", arcLenToNextDoca: " << arcLenToNextDoca << ", cos(ang): " << cosAngFTNtoAxis0 << std::endl;

    // Ok, the first selection is that the cluster to merge lies within an (elliptical) tube of the first cluster's axis
    if (docaVecProj1 < docaVecProj1Cut && docaVecProj2 < docaVecProj2Cut) consistent = true;
    
    // Otherwise we need to decide if the two clusters are consistent because they are "similar"...
    else
    {
        // Set up to find the distance of closeset approach of the two primary axes.
        // Results vectors
        TVector3 firstPoca;
        TVector3 nextPoca;
        
        // Get the primary axis for the next point
        TVector3 nextAxis0(nextPCA.getEigenVectors()[0][0],nextPCA.getEigenVectors()[0][1],nextPCA.getEigenVectors()[0][2]);
        
        // Convention on axis direction applied again
        if (firstPosToNextPos.Dot(nextAxis0) < 0.) nextAxis0 = -nextAxis0;
        
        // Recover the doca of the two axes and their points of closest approach
        float lineDoca = closestApproach(firstCenter, firstAxis0, nextCenter, nextAxis0, firstPoca, nextPoca);
        
        // Determine the arc lengths to the two pocas
        double arcLenToFirstPoca = (firstPoca - firstCenter).Dot(firstAxis0);  // is this faster than "Mag"?
        double arcLenToNextPoca  = (nextPoca  - nextCenter ).Dot(nextAxis0);
        
        std::cout << "       - arcLenToFirstPoca: " << arcLenToFirstPoca << ", arcLenToNextPoca: " << arcLenToNextPoca << ", first/Next dist: " << firstToNextDist << std::endl;
        
        // Require both of these to be less than length from first to next and to have the "right" sign where for the arc length
        // to the first axis poca this will be positive and for the next poca it will be negative
        // This prevents really long clusters that are not consistent from getting attached at their end points
        if (arcLenToFirstPoca >= 0. && arcLenToFirstPoca < firstToNextDist && arcLenToNextPoca <= 0. && arcLenToNextPoca > -firstToNextDist)
        {
            // Don't let clusters that are really far apart get joined together and really try to suppress clusters which
            // are not colinear but which have a small doca
            double nextEigenVal0     = std::max(1.,3. * sqrt(nextPCA.getEigenValues()[0]));
            double nextArcLenCut     = (1. + 5. * cosAngFTNtoAxis0) * nextEigenVal0;
            
            std::cout << "       - linedoca: " << lineDoca << ", cosAngFTNtoAxis0: " << cosAngFTNtoAxis0 << ", nextEigenVal0: " << nextEigenVal0 << ", nextArcLenCut: " << nextArcLenCut << std::endl;
            
            // Check the actual doca with a simple cut on the first eigen value, make sure "in range"
            if (lineDoca < firstEigenVals[1] && -arcLenToNextPoca < nextArcLenCut) consistent = true;
        }
    }
    
    return consistent;
    
    // hide all of this down here for now
/*
    // We need the next cluster's principal eigenvalue
    double nextEigenVal0 = 3. * std::sqrt(nextPCA.getEigenValues()[0]);
    
    // Now get position of endpoint of next cluster
    TVector3 nextEndPoint = nextCenter - nextEigenVal0 * nextAxis0;
    
    // Ok, now reset the firstPosToNextPos vector
    firstPosToNextPos = nextEndPoint - firstCenter;
    
    // Get the arclength to the point on the first axis which is the poca to the next point
    double arcLenToNextDoca = firstPosToNextPos.Dot(firstAxis0);
    
    // And recover that point
    TVector3 firstAxisPoca = firstCenter + arcLenToNextDoca * firstAxis0;
    
    // And then the vector for poca
    TVector3 nextToFirstDocaVec = nextEndPoint - firstAxisPoca;
    
    // Projections of the doca vec on the first and second eigen axes
    double docaVecProjMajor = std::fabs(nextToFirstDocaVec.Dot(firstAxis1));
    double docaVecProjMinor = std::fabs(nextToFirstDocaVec.Dot(firstAxis2));
    
    // Get the eigen values
    
    // Get the ratio of the first eigen value to the arc length
    double scaleFactor   = std::max(3.,arcLenToNextDoca / firstEigenVals[0]);
    double cutValueMajor = scaleFactor * firstEigenVals[1];
    double cutValueMinor = scaleFactor * firstEigenVals[2];
    
    
    // Ok, the first selection is that the doca of the next point to the first PCA primary axis is within
    // the allowed maximum spread
    if (docaVecProjMajor < cutValueMajor && docaVecProjMinor < cutValueMinor)
    {
        
        // Effectively the same cut as above but now with the point of closest approach
        float arcLenToPoca0 = (firstPoca - firstCenter).Dot(firstAxis0);
        float docaCutValue  = arcLenToPoca0 * firstEigenVals[1] / firstEigenVals[0];
        
        std::cout << "   --> arcLenToPoca0: " << arcLenToPoca0 << ", firstToNext: " << firstPosToNextPos.Mag() << ", docaCutValue: " << docaCutValue << ", lineDoca: " << lineDoca << std::endl;
        
        // Require that the doca of the two lines be within range
        if (lineDoca < 3. * docaCutValue)
        {
            // Final check, look at arc lengths to pocas
            double arcLenToPocaFirst     = (firstPoca - firstCenter).Dot(firstAxis0);
            double arcLenToPocaNext      = (nextPoca  - nextEndPoint).Dot(nextAxis0);
            double firstPosToNextPosDist = firstPosToNextPos.Mag();
            double totalArcLen           = arcLenToPocaFirst - arcLenToPocaNext;  // arcLenToPocaNext should be negative
            
            std::cout << "   >>> arcLenToPocaFirst: " << arcLenToPocaFirst << ", arcLenToPocaNext: " << arcLenToPocaNext << ", totalArcLen: " << totalArcLen << ", 2nd cut: " << firstPosToNextPos.Mag() / m_minCosAxisAng << std::endl;
            
            // It can be that that the lines are nearly parallel so the poca can be past the center of the two clusters
            // If that is the case then it should be that a simple angle check is sufficient
            if (arcLenToPocaFirst > firstPosToNextPosDist || arcLenToPocaNext > firstPosToNextPosDist)
            {
                if (firstPosToNextPos.Unit().Dot(firstAxis0) > m_minCosAxisAng) consistent = true;
            }
            // Otherwise the poca's are both "between" the cluster centers so require the sum of arc lengths to be within range
            else if (2. * totalArcLen >firstPosToNextPos.Mag() && totalArcLen < firstPosToNextPos.Mag() / m_minCosAxisAng) consistent = true;
        }
    }
*/
}
    
bool ClusterMergeAlg::mergeClusters(reco::ClusterParameters& firstClusterParams, reco::ClusterParameters& nextClusterParams) const
{
    bool merged(false);
    
    // Merge the next cluster into the first one
    // Do this by making a local copy of the input cluster parameters
    reco::ClusterParameters clusterParams = firstClusterParams;
    
    // Get the hits
    reco::HitPairListPtr& hitPairListPtr = clusterParams.getHitPairListPtr();
    
    for(const auto* hit : nextClusterParams.getHitPairListPtr())
    {
        hitPairListPtr.push_back(hit);
        
        for(const auto* hit2D : hit->getHits())
            if (hit2D) clusterParams.UpdateParameters(hit2D);
    }
    
    
    TVector3 origCenter(clusterParams.getFullPCA().getAvePosition()[0],clusterParams.getFullPCA().getAvePosition()[1],clusterParams.getFullPCA().getAvePosition()[2]);
    TVector3 origAxis0(clusterParams.getFullPCA().getEigenVectors()[0][0],clusterParams.getFullPCA().getEigenVectors()[0][1],clusterParams.getFullPCA().getEigenVectors()[0][2]);
    TVector3 origEigen(clusterParams.getFullPCA().getEigenValues()[0],clusterParams.getFullPCA().getEigenValues()[1],clusterParams.getFullPCA().getEigenValues()[2]);

    std::cout << "      **>> orig center: " << origCenter[0] << "/" << origCenter[1] << "/" << origCenter[2] << std::endl;
    std::cout << "           orig vector: " << origAxis0[0] << "/" << origAxis0[1] << "/" << origAxis0[2] << std::endl;
    std::cout << "           orig eigen:  " << origEigen[0] << "/" << origEigen[1] << "/" << origEigen[2] << std::endl;
    
    
    // Recalculate the PCA
    m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
    
    
    TVector3 newCenter(clusterParams.getFullPCA().getAvePosition()[0],clusterParams.getFullPCA().getAvePosition()[1],clusterParams.getFullPCA().getAvePosition()[2]);
    TVector3 newAxis0(clusterParams.getFullPCA().getEigenVectors()[0][0],clusterParams.getFullPCA().getEigenVectors()[0][1],clusterParams.getFullPCA().getEigenVectors()[0][2]);
    TVector3 newEigen(clusterParams.getFullPCA().getEigenValues()[0],clusterParams.getFullPCA().getEigenValues()[1],clusterParams.getFullPCA().getEigenValues()[2]);
    
    std::cout << "      >>>> new center:  " << newCenter[0] << "/" << newCenter[1] << "/" << newCenter[2] << std::endl;
    std::cout << "           new vector:  " << newAxis0[0] << "/" << newAxis0[1] << "/" << newAxis0[2] << std::endl;
    std::cout << "           new eigen:   " << newEigen[0] << "/" << newEigen[1] << "/" << newEigen[2] << ", cos orig to new: " << origAxis0.Dot(newAxis0) << std::endl;
    

    if (newEigen[0] > 0.8 * origEigen[0] && newEigen[1] * origEigen[0] < 5. * newEigen[0] * origEigen[1])
    {
        // Must have a valid pca
        if (clusterParams.getFullPCA().getSvdOK())
        {
            // Finish out the merging here
            reco::Hit3DToEdgeMap& firstEdgeMap = clusterParams.getHit3DToEdgeMap();
            reco::Hit3DToEdgeMap& nextEdgeMap  = nextClusterParams.getHit3DToEdgeMap();
            
            for(const auto& pair : nextEdgeMap) firstEdgeMap[pair.first] = pair.second;
            
            // Set the skeleton PCA to make sure it has some value
            clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
            
            // Copy back to the input
            firstClusterParams = clusterParams;
            
            merged = true;
        }
    }
    
    return merged;
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
