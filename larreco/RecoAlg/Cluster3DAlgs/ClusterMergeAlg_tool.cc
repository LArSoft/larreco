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
    float getTimeToExecute() const override {return fTimeToProcess;}
    
private:
    
    bool consistentClusters(const reco::PrincipalComponents&, const reco::PrincipalComponents&) const;
    
    bool linearClusters(const reco::PrincipalComponents&, const reco::PrincipalComponents&) const;

    bool mergeClusters(reco::ClusterParameters&, reco::ClusterParameters&) const;
    
    float closestApproach(const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, Eigen::Vector3f&, Eigen::Vector3f&, Eigen::Vector3f&) const;
    
    /**
     *  @brief Data members to follow
     */
    bool                                 fEnableMonitoring;       ///< If true then turn on monitoring (e.g. timing)
    float                                fAxisAngleScaleFactor;   ///< Addition scaling on the angle selection cut
    float                                fMinTransEigenVal;       ///< Set a mininum allowed value for the transverse eigen values
    float                                fMinEigenToProcess;      ///< Don't look anymore at clusters below this size
    mutable float                        fTimeToProcess;          ///< Keep track of how long it took to run this algorithm
    
    geo::Geometry*                       fGeometry;               //< pointer to the Geometry service
    
    PrincipalComponentsAlg               fPCAAlg;                 // For running Principal Components Analysis
};

ClusterMergeAlg::ClusterMergeAlg(fhicl::ParameterSet const &pset) :
    fPCAAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
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
    fEnableMonitoring     = pset.get<bool> ("EnableMonitoring",     true  );
    fAxisAngleScaleFactor = pset.get<float>("AxisAngleScaleFactor", 5.    );
    fMinTransEigenVal     = pset.get<float>("MinTransEigenVal",     0.09  );
    fMinEigenToProcess    = pset.get<float>("MinEigenToProcess",    2.0   );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    fGeometry = &*geometry;
    
    fTimeToProcess = 0.;
    
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
    if (fEnableMonitoring) theClockBuildClusters.start();
    
    // Resort by the largest first eigen value (so the "longest" cluster)
    clusterParametersList.sort([](auto& left, auto& right){return left.getFullPCA().getEigenValues()[0] > right.getFullPCA().getEigenValues()[0];});
    
    reco::ClusterParametersList::iterator firstClusterItr = clusterParametersList.begin();
    
    while(firstClusterItr != clusterParametersList.end())
    {
        reco::ClusterParameters&              firstClusterParams = *firstClusterItr;
        reco::ClusterParametersList::iterator nextClusterItr     = firstClusterItr;
        
        // Once you get down to the smallest clusters if they haven't already been absorbed there is no need to check them
        if (firstClusterParams.getFullPCA().getEigenValues()[0] < fMinEigenToProcess) break;
        
        // want the next one...
        nextClusterItr++;
    
        while(nextClusterItr != clusterParametersList.end())
        {
            reco::ClusterParameters& nextClusterParams = *nextClusterItr;
            
            // On any given loop through here it **might** be that the first cluster has been modified. So can't cache
            // the parameters, need the curret ones
            if (linearClusters(firstClusterParams.getFullPCA(),nextClusterParams.getFullPCA()))
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
    
    
    if (fEnableMonitoring)
    {
        theClockBuildClusters.stop();
        
        fTimeToProcess = theClockBuildClusters.accumulated_real_time();
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
    
    if (!consistent) return consistent;
    
    // Recover the positions of the centers of the two clusters
    const Eigen::Vector3f& firstCenter = firstPCA.getAvePosition();
    const Eigen::Vector3f& nextCenter  = nextPCA.getAvePosition();
    
    // And form a vector between the two centers
    Eigen::Vector3f firstPosToNextPosVec  = nextCenter - firstCenter;
//    Eigen::Vector3f firstPosToNextPosUnit = firstPosToNextPosVec.normalized();
    
    // Now get the first PCA's primary axis and since we'll use them get all of them at once...
    Eigen::Vector3f firstAxis0(firstPCA.getEigenVectors().row(0));
    Eigen::Vector3f firstAxis1(firstPCA.getEigenVectors().row(1));
    Eigen::Vector3f firstAxis2(firstPCA.getEigenVectors().row(2));

    // Want the distance of closest approach of the next cluser's center to the primary, start by finding arc length
    float arcLenToNextDoca = firstPosToNextPosVec.dot(firstAxis0);

    // Adopt the convention that the cluster axis is in same direction as vector from first to next centers
    // And preserve the overall orientation of the PCA by flipping all if we flip the first one
    if (arcLenToNextDoca < 0.)
    {
        firstAxis0       = -firstAxis0;
        firstAxis1       = -firstAxis1;
        firstAxis2       = -firstAxis2;
        arcLenToNextDoca = -arcLenToNextDoca;
    }
    
    // Position on the first cluster's axis of doca to next center
    Eigen::Vector3f firstAxisDocaPos = firstCenter + arcLenToNextDoca * firstAxis0;
    
    // Doca vector
    Eigen::Vector3f nextDocaVec = nextCenter - firstAxisDocaPos;
    
    // Need the projection of the doca vector onto the two transverse axes of the first cluster
    float docaVecProj1 = std::fabs(firstAxis1.dot(nextDocaVec));
    float docaVecProj2 = std::fabs(firstAxis2.dot(nextDocaVec));
    
    // Recover the eigen values of the first and second PCAs for selection cuts
    Eigen::Vector3f firstEigenVals(3. * sqrt(firstPCA.getEigenValues()[0]),
                                   3. * sqrt(firstPCA.getEigenValues()[1]),
                                   3. * sqrt(firstPCA.getEigenValues()[2]));
    
    Eigen::Vector3f nextEigenVals( 3. * sqrt(nextPCA.getEigenValues()[0]),
                                   3. * sqrt(nextPCA.getEigenValues()[1]),
                                   3. * sqrt(nextPCA.getEigenValues()[2]));
    
   // Use the angle between the vector between cluster centers and the first axis to moderate the selection cut
    float firstToNextDist  = firstPosToNextPosVec.norm();
    float cosAngFTNtoAxis0 = arcLenToNextDoca / firstToNextDist;  // by construction arcLenToNextDoca is projecton of total dist on first axis
    float docaVecProj1Cut  = std::min(50., std::max( 1., (1. + cosAngFTNtoAxis0) * firstEigenVals[1]));
    float docaVecProj2Cut  = std::min(40., std::max(0.5, (1. + cosAngFTNtoAxis0) * firstEigenVals[2]));

    // Ok, the first selection is that the cluster to merge lies within an (elliptical) tube of the first cluster's axis
    if (firstEigenVals(0) > 11000. * nextEigenVals(0) && docaVecProj1 < docaVecProj1Cut && docaVecProj2 < docaVecProj2Cut)
    {
        consistent = true;
    }
    // Otherwise we need to decide if the two clusters are consistent because they are "similar"...
    else
    {
        // Get the primary axis for the next point
        Eigen::Vector3f nextAxis0(nextPCA.getEigenVectors().row(0));
        
        // Convention on axis direction applied again
        if (firstPosToNextPosVec.dot(nextAxis0) < 0.) nextAxis0 = -nextAxis0;
        
        // It should be that the two primary vectors point in the same direction... so now we can select on angle
        if (firstAxis0.dot(nextAxis0) > 0.5)
        {
            // Set up to find the distance of closeset approach of the two primary axes.
            // Results vectors
            Eigen::Vector3f firstPoca;
            Eigen::Vector3f nextPoca;
            Eigen::Vector3f firstNextVec;
            
            // Recover the doca of the two axes and their points of closest approach along each axis
            float lineDoca = closestApproach(firstCenter, firstAxis0, nextCenter, nextAxis0, firstPoca, nextPoca, firstNextVec);
            
            // Same sort of logic as before, the distance of closest approach needs to be sensible
            float firstNextProj1 = std::abs(firstAxis1.dot(firstNextVec));
            float firstNextProj2 = std::abs(firstAxis2.dot(firstNextVec));
            
            if (firstNextProj1 < 5. * firstEigenVals(1) && firstNextProj2 < 5. * firstEigenVals(2))
            {
                // Determine the arc lengths to the two pocas
                // Note that we do it this way so the arc lengths will be signed (which is important!)
                float arcLenToFirstPoca = (firstPoca - firstCenter).dot(firstAxis0);
                float arcLenToNextPoca  = (nextPoca  - nextCenter ).dot(nextAxis0);
                
                // Require both of these to be less than length from first to next and to have the "right" sign where for the arc length
                // to the first axis poca this will be positive and for the next poca it will be negative
                // This prevents really long clusters that are not consistent from getting attached at their end points
                if (arcLenToFirstPoca >= 0. && arcLenToFirstPoca < firstToNextDist && arcLenToNextPoca <= 0. && arcLenToNextPoca > -firstToNextDist)
                {
                    // Don't let clusters that are really far apart get joined together and really try to suppress clusters which
                    // are not colinear but which have a small doca
                    float nextEigenVal0     = std::max(1.,3. * sqrt(nextPCA.getEigenValues()[0]));
                    float nextArcLenCut     = (1. + 5. * cosAngFTNtoAxis0) * nextEigenVal0;
                    
                    // Check the actual doca with a simple cut on the first eigen value, make sure "in range"
                    if (lineDoca < firstEigenVals[1] && -arcLenToNextPoca < nextArcLenCut) consistent = true;
                }
            }
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
    
bool ClusterMergeAlg::linearClusters(const reco::PrincipalComponents& firstPCA, const reco::PrincipalComponents& nextPCA) const
{
    // Assume failure
    bool consistent(false);
    
    // The goal here is to compare the two input PCA's and determine if they are effectively colinear and
    // within reasonable range to consider merging them. Note that a key assumption is that the first input
    // PCA is from the "bigger" cluster, the second is "smaller" and may have a less reliable PCA.
    
    // Recover the positions of the centers of the two clusters
    const Eigen::Vector3f& firstCenter = firstPCA.getAvePosition();
    const Eigen::Vector3f& nextCenter  = nextPCA.getAvePosition();
    
    // And form a vector between the two centers
    Eigen::Vector3f firstPosToNextPosVec  = nextCenter - firstCenter;
    Eigen::Vector3f firstPosToNextPosUnit = firstPosToNextPosVec.normalized();
    
    // Now get the first PCA's primary axis and since we'll use them get all of them at once...
    Eigen::Vector3f firstAxis0(firstPCA.getEigenVectors().row(0));
    Eigen::Vector3f firstAxis1(firstPCA.getEigenVectors().row(1));
    Eigen::Vector3f firstAxis2(firstPCA.getEigenVectors().row(2));
    
    // Want the distance of closest approach of the next cluser's center to the primary, start by finding arc length
    float arcLenToNextDoca = firstPosToNextPosVec.dot(firstAxis0);
    
    // Adopt the convention that the cluster axis is in same direction as vector from first to next centers
    // And preserve the overall orientation of the PCA by flipping all if we flip the first one
    if (arcLenToNextDoca < 0.)
    {
        firstAxis0       = -firstAxis0;
        firstAxis1       = -firstAxis1;
        firstAxis2       = -firstAxis2;
        arcLenToNextDoca = -arcLenToNextDoca;
    }
    
    // Recover the eigen values of the first and second PCAs for selection cuts
    Eigen::Vector3f firstEigenVals(2.0 * sqrt(firstPCA.getEigenValues()[0]),
                                   2.0 * sqrt(std::max(firstPCA.getEigenValues()[1],fMinTransEigenVal)),
                                   2.0 * sqrt(std::max(firstPCA.getEigenValues()[2],fMinTransEigenVal)));
    
    // Position on the first cluster's axis of doca to next center
    Eigen::Vector3f firstAxisDocaPos = firstCenter + arcLenToNextDoca * firstAxis0;
    
    // And now we can compute the doca vector...
    Eigen::Vector3f nextDocaVec = nextCenter - firstAxisDocaPos;
    
    // And with this we can now get the projections of this doca vector onto the other two first PCA axes
    float nextDocaVecProj1 = std::fabs(firstAxis1.dot(nextDocaVec));
    float nextDocaVecProj2 = std::fabs(firstAxis2.dot(nextDocaVec));
    
    // Handle the special case that the cluster we are comparing to is "embedded" in the current cluster
    // (and could end up not passing more detailed cuts below)
    if (arcLenToNextDoca < firstEigenVals[0] && nextDocaVecProj1 < firstEigenVals[1] && nextDocaVecProj2 < firstEigenVals[2])
    {
        consistent = true;
    }
    // Otherwise we need to investigate further
    else
    {
        // Get a scaling factor based on the ratio of the projected distance to the next cluster
        float angleScaleFactor = fAxisAngleScaleFactor * std::max(arcLenToNextDoca / firstEigenVals(0),float(1.));
        
        // The first selection is then that these projects are "within range"
        // Here "range" is determined by the eigen values of those two axes but note we need to have a cutoff...
        float rangeAxis1 = angleScaleFactor * firstEigenVals(1);
        float rangeAxis2 = angleScaleFactor * firstEigenVals(2);
        
        if (nextDocaVecProj1 < rangeAxis1 && nextDocaVecProj2 < rangeAxis2)
        {
            // Now we recover the doca of the two primary axes of the input clusters
            // Get the primary axis for the next point
            Eigen::Vector3f nextAxis0(nextPCA.getEigenVectors().row(0));
            Eigen::Vector3f nextAxis1(nextPCA.getEigenVectors().row(1));
            Eigen::Vector3f nextAxis2(nextPCA.getEigenVectors().row(2));
            
            // Get the arc length from next to first centers - note convention for
            // PCA orientation
            float arcLenToFirstDoca = -firstPosToNextPosVec.dot(nextAxis0);
            
            // Convention on axis direction applied again
            if (arcLenToFirstDoca > 0.)
            {
                nextAxis0         = -nextAxis0;
                nextAxis1         = -nextAxis1;
                nextAxis2         = -nextAxis2;
                arcLenToFirstDoca = -arcLenToFirstDoca;
            }
            
            // Get the eigen values again
            Eigen::Vector3f nextEigenVals(2.0 * sqrt(nextPCA.getEigenValues()[0]),
                                          2.0 * sqrt(std::max(nextPCA.getEigenValues()[1],fMinTransEigenVal)),
                                          2.0 * sqrt(std::max(nextPCA.getEigenValues()[2],fMinTransEigenVal)));
            
            angleScaleFactor = fAxisAngleScaleFactor * std::max(-arcLenToFirstDoca / nextEigenVals(0),float(1.));
            
            // Position on the first cluster's axis of doca to next center
            Eigen::Vector3f nextAxisDocaPos = nextCenter + arcLenToFirstDoca * nextAxis0;
            
            // And now we can compute the doca vector...
            Eigen::Vector3f firstDocaVec = firstCenter - nextAxisDocaPos;
            
            // And with this we can now get the projections of this doca vector onto the other two first PCA axes
            float firstDocaVecProj1 = std::fabs(nextAxis1.dot(firstDocaVec));
            float firstDocaVecProj2 = std::fabs(nextAxis2.dot(firstDocaVec));
            
            // The first selection is then that these projects are "within range"
            // Here "range" is determined by the eigen values of those two axes but note we need to have a cutoff...
            rangeAxis1 = angleScaleFactor * nextEigenVals(1);
            rangeAxis2 = angleScaleFactor * nextEigenVals(2) * 2.;
            
            if (firstDocaVecProj1 < rangeAxis1 && firstDocaVecProj2 < rangeAxis2)
            {
                // Closest approach calculaiton results vectors
                Eigen::Vector3f firstPoca;
                Eigen::Vector3f nextPoca;
                Eigen::Vector3f firstNextVec;
                
                // Recover the doca of the two axes and their points of closest approach along each axis
                float lineDoca = closestApproach(firstCenter, firstAxis0, nextCenter, nextAxis0, firstPoca, nextPoca, firstNextVec);
                
                // Its pointless to continue if the line doca is too large
                if (lineDoca < firstEigenVals(0) && firstAxis0.dot(nextAxis0) > 0.7)
                {
                    // We'll use the average of the two pocas for determining the bend
                    Eigen::Vector3f aveDocaPos = 0.5 * (firstPoca + nextPoca);
                    
                    // We're aiming to compute distance from the line joining the two centers to thie point
                    float arcLenToAveDoca = (aveDocaPos - firstCenter).dot(firstPosToNextPosUnit);
                    
                    Eigen::Vector3f aveDocaVec = aveDocaPos - firstCenter - arcLenToAveDoca * firstPosToNextPosUnit;
                    
                    float bendDist = aveDocaVec.norm();
                    
                    // Now get projections of the two pocas on the axis from first to second center
                    float projectionFirstPoca = (firstPoca - firstCenter).dot(firstPosToNextPosUnit);
                    float projectionNextPoca  = (nextPoca  - nextCenter).dot(firstPosToNextPosUnit);
                    
                    // Note by conventions chosen, the projection to the next poca will be negative.
                    // As well, we expect the first projection to be positive
                    if (bendDist < 20. && projectionFirstPoca > 0. && projectionNextPoca < 0.)
                    {
                        // Last up try to compute the gap between the two clusters. This will be taken as
                        // the projections of the eigenvectors of length given by their eigenvalues along
                        // the vector between the two centers
                        float projectionFirstClus  = firstEigenVals(0) * firstPosToNextPosUnit.dot(firstAxis0);
                        float projectionNextClus   = nextEigenVals(0) * firstPosToNextPosUnit.dot(nextAxis0);
                        float firstPosToNextPosLen = firstPosToNextPosVec.norm();
                        
                        float gap = firstPosToNextPosLen - projectionFirstClus - projectionNextClus;
                        
                        if (gap < firstEigenVals(0)) consistent = true;
                    }
                }
            }
        }
    }
    
    return consistent;
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
    
    
//    const Eigen::Vector3f& origCenter = clusterParams.getFullPCA().getAvePosition();
//    Eigen::Vector3f        origAxis0(clusterParams.getFullPCA().getEigenVectors().row(0));
//    Eigen::Vector3f        origEigen(clusterParams.getFullPCA().getEigenValues());

//    std::cout << "      **>> orig center: " << origCenter[0] << "/" << origCenter[1] << "/" << origCenter[2] << std::endl;
//    std::cout << "           orig vector: " << origAxis0[0] << "/" << origAxis0[1] << "/" << origAxis0[2] << std::endl;
//    std::cout << "           orig eigen:  " << origEigen[0] << "/" << origEigen[1] << "/" << origEigen[2] << std::endl;

    
    // Recalculate the PCA
    fPCAAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
    
//    const Eigen::Vector3f& newCenter = clusterParams.getFullPCA().getAvePosition();
//    Eigen::Vector3f        newAxis0(clusterParams.getFullPCA().getEigenVectors().row(0));
//    Eigen::Vector3f        newEigen(clusterParams.getFullPCA().getEigenValues());
    
//    std::cout << "      >>>> new center:  " << newCenter[0] << "/" << newCenter[1] << "/" << newCenter[2] << std::endl;
//    std::cout << "           new vector:  " << newAxis0[0] << "/" << newAxis0[1] << "/" << newAxis0[2] << std::endl;
//    std::cout << "           new eigen:   " << newEigen[0] << "/" << newEigen[1] << "/" << newEigen[2] << ", cos orig to new: " << origAxis0.Dot(newAxis0) << std::endl;
    

//    if (newEigen[0] > 0.8 * origEigen[0] && newEigen[1] * origEigen[0] < 5. * newEigen[0] * origEigen[1])
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
    
float ClusterMergeAlg::closestApproach(const Eigen::Vector3f& P0, const Eigen::Vector3f& u0,
                                       const Eigen::Vector3f& P1, const Eigen::Vector3f& u1,
                                       Eigen::Vector3f&       poca0,
                                       Eigen::Vector3f&       poca1,
                                       Eigen::Vector3f&       firstNextVec) const
{
    // Technique is to compute the arclength to each point of closest approach
    Eigen::Vector3f w0 = P0 - P1;
    float a(1.);
    float b(u0.dot(u1));
    float c(1.);
    float d(u0.dot(w0));
    float e(u1.dot(w0));
    float den(a * c - b * b);
    
    float arcLen0 = (b * e - c * d) / den;
    float arcLen1 = (a * e - b * d) / den;
    
    poca0 = P0 + arcLen0 * u0;
    poca1 = P1 + arcLen1 * u1;
    
    firstNextVec = poca1 - poca0;
    
    return firstNextVec.norm();
}
    
DEFINE_ART_CLASS_TOOL(ClusterMergeAlg)
} // namespace lar_cluster3d
